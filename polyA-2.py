import argparse
import logging
import pickle
import os
import sys
import re
# External modules below
import pysam
# Custom modules below
import mapping
import transcript
import pybedtools

parser = argparse.ArgumentParser(description='this program tries to find polya cleavage sites through short-read assembly.it is expected that contigs are aligned to contigs, and reads aligned to contigs. these 2 alignment steps can be performed by trans-abyss. the aligners used are gmap for contig-genome and bwa-sw for read-contig alignments. annotations files for ensembl, knowngenes, refseq, and aceview are downloaded from ucsc. est data(optional) are also downloaded from ucsc. the analysis can be composed of 2 phases: 1. contig-centric phase - cleavage sites per contig are captured 2. coordinate-centric phase - contigs capturing the same cleavage site are consolidated into 1 report where expression/evidence-related data are summed. customized filtering based on evidence data can be performed.')
parser.add_argument('c2g', metavar='<contig-to-genome>', help='the contig-to-genome alignment file in bam format')
parser.add_argument('contigs', metavar='<contigs>', help='the file containing all contigs in fasta format')
parser.add_argument('ref_genome', metavar='<reference_genome>', help='the reference genome to use. default is hg19.', default='hg19', choices=['hg18','hg19'])
parser.add_argument('annot', metavar='<annotations>', help='the annotations file to use with the reference in gtf format.')
parser.add_argument('r2c', metavar='<reads-to-contigs>', help='the contigs-to-genome alignment file')
parser.add_argument('out', metavar='<output-file>', help='the file to output results')
parser.add_argument('-r', '--output_reads', dest='output_reads', help='enable to output read sequences', action='store_true', default=False)
parser.add_argument('-t', '--trim', dest='trim_reads', help='trim bases of quality <= this value. default is 3.', type=int, default=3)
parser.add_argument('--use_tmp', help='use tmp space', action='store_true', default=False)
parser.add_argument('--no_link', help='do not find link pairs', action='store_true', default=False)
parser.add_argument('-e', '--overlap_est', help='overlap expressed sequence tags', action='store_true', default=False)
parser.add_argument('--min_at', help='minimum number of a|t bases in tail. default is 4.', type=int, default=4)
parser.add_argument('--max_diff', help='maximum rate of x to y bases allowed in tail. where x are non a|t bases and y are a|t bases. default: 1 5.', nargs=2, default=[1,5])
parser.add_argument('--max_diff_link', help='Maximum number of non A|T bases in entire link read. Default is 2.', type=int, default=2)
parser.add_argument('--min_bridge_size', help='Minimum size of bridge. Default is 1.', type=int, default=1)
parser.add_argument('-k', '--track', metavar=('[name]','[description]'), help='Name and description of BED graph track to output.', nargs=2)
parser.add_argument('--rgb', help='RGB value of BED graph. Default is 0,0,255', default='0,0,255')

args = parser.parse_args()
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger('polyA_logger')

fh = logging.FileHandler('info.log')
fh.setLevel(logging.INFO)
logger.addHandler(fh)

# Globally used variables

# Reference genome sequence (refseq)
if (args.ref_genome == 'hg19'):
    reference_genome = '/projects/btl/trans-abyss/public_releases/v1.4.8/annotations/hg19/hg19.fa'
logger.debug("Loading reference genome via pysam 0.8.1")
refseq = pysam.FastaFile(reference_genome)
logger.debug("Reference genome successfully loaded!")

# Contigs to genome alignments (aligns)
logger.debug("Checking contig to genome alignment is in BAM format. Splitting file by '.': %s", args.c2g.split('.'))
if (args.c2g.split('.')[-1] == 'bam'):
    logger.debug("Contig to genome alignment file is in BAM format, loading BAM file...")
    aligns = pysam.AlignmentFile(args.c2g, "rb")
    logger.debug("Contig to genome alignment file loaded successfully")
else:
    sys.exit("Unrecognized format for <aln>, must be BAM format file. Exiting.")

# Features (features)
features = pysam.TabixFile(args.annot, parser=pysam.asGTF())
#features = pysam.tabix_iterator(open(args.annot), pysam.asGTF())

# Reads to contigs alignment (r2c)
r2c = pysam.AlignmentFile(args.r2c, "rb")

def findTailContig(align, contig, min_len, max_nonAT):
    junction_buffer = 50
    # clipped represents [start, end]
    clipped = [0,0]
    if int(align.qstart) > 1:
        clipped[0] = 1
    if int(align.qend) < int(align.query_len):
        clipped[1] = 1
    if (clipped[0] == 1):
        last_matched = int(align.qstart)
        clipped_seq = contig.query[:int(align.qstart)-1]
        junction_seq = contig.query[:int(align.qstart)-1+junction_buffer]
    else:
        last_matched = int(align.qend)
        clipped_seq = contig.query[int(align.qend):]
        junction_seq = contig.query[int(align.qend)-junction_buffer]

def qposToTpos(align, qpos):
    blocks = align.blocks
    for i in xrange(len(blocks)):
        qb0 = blocks[i][0]
        qb1 = blocks[i][1]
        if ((qpos >= qb0) and (qpos <= qb1)) or ((qpos <= qb0) and (qpos >= qb1)):
            block = i
            break
    tpos = None
    if block is not None:
        if (align.is_reverse == False):
            tpos = blocks[block][0] + qpos - blocks[block][0]
        else:
            tpos = blocks[block][1] - (qpos - blocks[block][1])
    return tpos

def isPolyATail(seq, min_len, max_nonAT):
    result = True
    nonATs = len(seq)-seq.count('A')-seq.count('T')
    if (seq is None) or (seq == '') or (nonATs > min_len):
        return False
    if (float(max_nonAT)[0]/max_nonAT[1]) < (float(nonATs)/(len(seq)-nonATs)):
        return False

def getQueryLenFromCigar(cigar):
    lens = []
    # if op is not 'D'(deletion), 'N'(skipped region), 'P'(padding)
    query_lens = [int(op[1]) for op in cigar if op[0] != 2 and op[0] != 3 and op[0] != 6]
    return sum(query_lens)

def cigarToBlocks(cigar, tstart, strand):
    query_len = getQueryLenFromCigar(cigar)
    qstart = 1 if strand == '+' else query_len
    tblocks = []
    qblocks = []          
    for i in range(len(cigar)):
        op, length = cigar[i]
        # 'D' (deletion)
        if op == 2 and i == 0:
            return None, None
        # 'S' or 'H' (clips)
        if op == 4 or op == 5:
            if i == 0:
                qstart = qstart + length if strand == '+' else qstart - length
                continue
        tblock = None
        qblock = None
        if not tblocks and op != 0:
            return None, None
        # match
        if op == 0:
            tend = tstart + length - 1
            qend = qstart + length - 1 if strand == '+' else qstart - length + 1
            tblock = [tstart, tend]
            qblock = [qstart, qend]
        # intron ('N'), skipped reference or deletion in reference ('D')
        elif op == 2 or op == 3:
            #intron
            if op == 3 and length < 3:
                continue
            qend = qend
            tend = tstart + length - 1 
        # insertion ('I') to reference
        elif op == 1:
            tend = tend
            qend = qstart + length - 1 if strand == '+' else qstart - length + 1
        if tblock:
            tblocks.append(tblock)
        if qblock:
            qblocks.append(qblock)
        tstart = tend + 1
        qstart = qend + 1 if strand == '+' else qend - 1
    return tblocks, qblocks

def revComp(seq):
    ndic = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
    revcomp = ''
    for nuc in seq:
        revcomp += ndic[nuc]
    return revcomp[::-1]

def getQueryBlocks(align):
    query_blocks = []
    start = align.qstart
    for block in align.blocks:
        diff = block[1]-block[0]
        query_blocks.append([start,start+diff])
        start = start+diff
    return query_blocks

def inferStrand(align, overlapping_features):
    half = len(overlapping_features)/2
    i = 0
    strand = None
    for feat in overlapping_features:
        if (feat.feature == '3UTR'):
            if (i < half):
                strand = '+'
            else:
                strand = '-'
            return strand
        elif (feat.feature == '5UTR'):
            if (i < half):
                strand = '-'
            else:
                strand = '+'
            return strand
        else:
            return None
        i += 1

# Iterate contig to genome alignments
for align in aligns:
    # Check if has polyA tail
    # ...
    # Obtain reads spanning contig
    chrom = aligns.getrname(align.tid)
    print 'Looking at contig: {}-{}-{}'.format(align.qname, align.reference_start, align.reference_end)
    read_count = 0
    strand = '+' if align.is_reverse == True else '-'
    
    # Check both ends of contig for 3UTR
    

#   print 'blocks: {}'.format(align.blocks)
#   print 'tblocks: {}'.format(cigarToBlocks(align.cigar, align.reference_start+1, strand)[0])
#   print 'qblocks: {}'.format(cigarToBlocks(align.cigar, align.reference_start+1, strand)[1])
#   print 'mqblocks: {}'.format(getQueryBlocks(align))
    for overlapping_read in r2c.fetch(align.qname, align.qstart, align.qend):
        read_count += 1
    feats = features.fetch(chrom, align.reference_start, align.reference_end)
    feature_list = []
    for feat in feats:
        feature_list.append(feat)
    utr3 = None
    for feature in feature_list:
        if (feature.feature == '3UTR'):
            if (utr3 is None):
                utr3 = feature
            elif (abs(feature.start - align.reference_start) < abs(utr3.start - align.reference_start)):
                utr3 = feature
        #print '\tFeature: {}\t{}\t{}\t{}'.format(feature.transcript_id,feature.feature,feature.start,feature.end)
    print '\tutr3: {}\t{}\t{}\t{}'.format(utr3.transcript_id, utr3.feature, utr3.start, utr3.end)
    inf_strand = inferStrand(align, feature_list)
    if (inf_strand is None):
        print 'Strand information could not be found nor inferred!'
    if (inf_strand == '+'):
        cleavage_site = align.reference_start
        if (strand == '+'):
            clipped_pos = 'start'
            last_matched = align.qstart
        else:
            clipped_pos = 'end'
            last_matched = align.qend
    else:
        cleavage_site = align.reference_end
        if (strand == '+'):
            clipped_pos = 'end'
            last_matched = align.qend
        else:
            clipped_pos = 'start'
            last_matched = align.qstart
    print '-'*50

def cigarScore(cigar_string, start, end):
    cigar = re.findall('(\d+)(\D+)', cigar_string)
    pos = 0
    rval = []
    for i in xrange(len(cigar)):
        if (pos < start <= pos+int(cigar[i][0])):
            rval.append(i)
        if (pos < end <= pos+int(cigar[i][0])):
            rval.append(i)
        pos += int(cigar[i][0])
    return rval
