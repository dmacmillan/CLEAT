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
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('polyA_logger')

fh = logging.FileHandler('info.log')
fh.setLevel(logging.INFO)
logger.addHandler(fh)

# Globally used variables

# Filters (filters)
filters = {}
filters['min_at'] = args.min_at
filters['max_diff'] = args.max_diff
filters['max_diff_link'] = args.max_diff_link

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

def int_to_base_qual(qual_int, offset):
    """Converts integer to base quality base"""
    if offset == 64 or offset == 33:
        return chr(qual_int + offset)

# Poor quals (poor_quals)
poor_quals = ''
if args.trim_reads:
    for i in range(args.trim_reads + 1):
        poor_qual = int_to_base_qual(i, 33)
        if poor_qual is not None:
            poor_quals += poor_qual

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

def check_freq(seq):
    """Returns frequency of each base in given sequence"""
    freq = {}
    for nt in seq:
        if not freq.has_key(nt.upper()):
            freq[nt.upper()] = 0
        freq[nt.upper()] += 1
    return freq

def is_polyA_tail(seq, expected_base, min_len, max_nonAT_allowed):
        """Determines if sequence can be possible tail
        
        min_len = minimum number of expected base in seq
        max_nonAT_allowed(N,M) = N base(s) other than expected base allowed
                                 per stretch of M bases
        """
        result = True
    
        if seq is None or seq == '' or expected_base is None or not expected_base in seq:
            return False
            
        # minimum number of As or Ts
        freq = check_freq(seq)
        if freq[expected_base] < min_len:
            return False

        for i in range(0, len(seq), max_nonAT_allowed[1]):
            subseq = seq[i:i+10]
        
            freq = self.check_freq(subseq)
            non_expected_freq = 0
            for base, f in freq.iteritems():
                if base.upper() != expected_base.upper():
                    non_expected_freq += f
            
            if non_expected_freq > max_nonAT_allowed[0]:
                result = False
                
        return result

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
                strand = '-'
            else:
                strand = '+'
            return strand
        elif (feat.feature == '5UTR'):
            if (i < half):
                strand = '+'
            else:
                strand = '-'
            return strand
        else:
            return None
        i += 1

def find_link_pairs(align, target, clipped_pos, last_matched, cleavage_site, feature, homo_len=20, max_mismatch=0):
    """Finds reads pairs where one mate is mapped to contig and the other mate (mapped elsewhere or unmapped)
    is a potential polyA tail
    
    The mate that is mapped to the contig (anchor) should be pointing outwards to the edge of the contig.
    The mate that is potential polyA tail can be unmapped or mapped to another contig.  If the mate is unmapped,
    it's expected to be stored under the same contig.  Unfortunately this is not the behaviour of BWA-SW so all the 
    unmapped mates are ignored for BWA-SW alignments.
    
    A matching transcript is provided.  The only place this is used is in the check of whether the vicinity of the 
    cleavage site has a homopolyer run.  The transcript strand is used for deciding whether the upstream or downstream
    region of the cleavage site should be checked.  If a polyT or polyA is in the neighborhood (200bp), then the case
    won't be further explored.
    """
    print '@{}\t{}\t{}\t{}\t{}\t{}'.format(align.qname, target, clipped_pos, last_matched, cleavage_site, feature)
    # determine if 3'UTR is first or last query block
    anchor_read_strand = None
    if clipped_pos == 'start':
        anchor_read_strand = '-'
    elif clipped_pos == 'end':
        anchor_read_strand = '+'
    
    if anchor_read_strand is None:
        return []
    
    # check if genomic region has polyA - if so, no good
    genome_buffer = 200
    if feature['feature'].strand == '-':
        span = (int(cleavage_site) - genome_buffer, int(cleavage_site))
    else:
        span = (int(cleavage_site), int(cleavage_site) + genome_buffer)
    #genome_seq = self.refseq.GetSequence(align.target, span[0], span[1])
    genome_seq = refseq.fetch(target, span[0], span[1])
    if re.search('A{%s,}' % (homo_len), genome_seq, re.IGNORECASE) or re.search('T{%s,}' % (homo_len), genome_seq, re.IGNORECASE):
        sys.stdout.write('genome sequence has polyAT tract - no reliable link pairs can be retrieved %s %s %s:%s-%s\n' % 
                         (align.query, cleavage_site, align.target, span[0], span[1]))
        return []

    mate_loc = {}
    #for read in self.bam.bam.fetch(align.query):
    cigarskip, anchor_rs_plus, anchor_rs_minus, read_count = 0,0,0,0
    for read in r2c.fetch(align.qname):
        read_count += 1
        # skip when both mates mapped to same contig
        if not read.mate_is_unmapped and read.tid == read.rnext:
            continue
        
        # anchor read must be mapped entirely within contig
        if len(read.cigar) > 1:
            cigarskip += 1
            continue
        
        if anchor_read_strand == '+' and (read.is_reverse or read.pos + 1 > last_matched):
            anchor_rs_plus += 1
            continue
        if anchor_read_strand == '-' and (not read.is_reverse or read.pos + read.rlen < last_matched):
            anchor_rs_minus += 1
            continue
        
        if read.rnext >= 0:
            mate_contig = r2c.getrname(read.rnext)       
            if not mate_loc.has_key(mate_contig):
                mate_loc[mate_contig] = {}
            mate_loc[mate_contig][read.qname] = read
        else:
            print 'cannot find unmapped mate %s' % read.qname
            
    print 'Total reads {}'.format(read_count)
    print '\tSkipped {} due to cigar > 1'.format(cigarskip)
    print '\tSkipped {} due to anchor read strand == +'.format(anchor_rs_plus)
    print '\tSkipped {} due to anchor read strand == -'.format(anchor_rs_minus)
    print 'Done looking at reads'
    link_pairs = []
    f = open('./'+align.qname+'.contigs','w')
    for contig in mate_loc.keys():
        f.write(contig+'\n')
        #for read in self.bam.bam.fetch(contig):
        for read in r2c.fetch(contig):
            if not mate_loc[contig].has_key(read.qname):
                continue
            
            trimmed_seq = read.seq
            if args.trim_reads:
                trimmed_seq = trim_bases(read.seq, read.qual)

            if trimmed_seq:
                for base in ('A', 'T'):
                    if is_bridge_read_good(trimmed_seq, base, len(trimmed_seq) - max_mismatch, mismatch=[max_mismatch, len(trimmed_seq)]):
                        link_pairs.append([read, mate_loc[contig][read.qname], trimmed_seq])
                        break
    
    return link_pairs

def is_bridge_read_good(clipped_seq, base, min_len, mismatch):
        """Determines if clipped sequence is possible polyA tail
        
        If clipped_seq is composed of base only, then it is automatically 
        considered a potential bridge read regardless of length
        Otherwise will check the frequecy of 'the other bases' using is_polyA_tail()
        to determine whether it's acceptable
        """
        good = False
        if clipped_seq[0].upper() == base and len(re.sub(r'(.)\1+', r'\1', clipped_seq)) == 1:
            good = True
        elif is_polyA_tail(clipped_seq, base, min_len=min_len, max_nonAT_allowed=mismatch):
            good = True

        return good

def trim_bases(seq, qual, end=None):
        """Trim poor quality bases from read sequence"""
        if poor_quals is None or poor_quals == '':
            return seq
        
        match_end = match_start = None
        match_end = re.search(r'[%s]+$' % poor_quals, qual)            
        match_start = re.search(r'^[%s]+' % poor_quals, qual)
        
        if match_start and not match_end:
            if end is None or end == 'start':
                return seq[match_start.end():]
            
        elif match_end and not match_start:
            if end is None or end == 'end':
                return seq[:match_end.start()]
            
        elif match_start and match_end:
            if end == 'start':
                return seq[match_start.end():]
            
            elif end == 'end':
                return seq[:match_end.start()]
            
            else:
                if len(match_start.group()) > len(match_end.group()):
                    return seq[match_start.end():]
                elif len(match_end.group()) > len(match_start.group()):
                    return seq[:match_end.start()]
        
        return seq

def inferCleavageSite(align, inf_strand):
    if (inf_strand == '-'):
        cleavage_site = align.reference_start
        if (align.is_reverse == True):
            clipped_pos = 'start'
            last_matched = align.qstart
        else:
            clipped_pos = 'end'
            last_matched = align.qend
    else:
        cleavage_site = align.reference_end
        if (align.is_reverse == False):
            clipped_pos = 'end'
            last_matched = align.qend
        else:
            clipped_pos = 'start'
            last_matched = align.qstart
    return [cleavage_site, clipped_pos, last_matched]

def show_trimmed_bases(seq, trimmed_bases):
    """Shows (link) read sequence with trimmed bases"""
    match_start = re.search('^' + trimmed_bases, seq)
    match_end = re.search(trimmed_bases + '$', seq)
    
    if match_start or match_end:
        if match_start:
            match = match_start
        else:
            match = match_end
        return seq[:match.start()].lower() + seq[match.start():match.end()].upper() + seq[match.end():].lower()
    
    else:
        return seq.lower()

def output_link_pairs(align, reads):
    """Outputs link pairs in FASTA format"""
    out = ''
    for r in reads:
        trimmed_seq = show_trimmed_bases(r[0].seq, r[2])
        num_trimmed_bases = r[0].rlen - len(r[2])
        if r[1].is_reverse:
            anchor_direction = 'L'
        else:
            anchor_direction = 'R'
        if r[0].is_unmapped:
            mate_contig = 'unmapped'
        else:
            mate_contig = r2c.getrname(r[0].tid)
        out += '>%s %s %s %s %s trimmed:%s\n%s\n' % (r[0].qname, align.query, r[1].pos, anchor_direction, 
                                                     mate_contig, num_trimmed_bases, trimmed_seq) 
    return out

# Iterate contig to genome alignments
for align in aligns:
    print 'Looking at contig: {}: {}-{}'.format(align.qname,align.reference_start, align.reference_end)
    # Get chromosome
    chrom = aligns.getrname(align.tid)
    # Find all overlapping features in gtf annotation
    feats = features.fetch(chrom, align.reference_start, align.reference_end)
    feature_list = []
    for feat in feats:
        feature_list.append({'feature': feat})
    # Try to get strandedness of alignment
    if (align.is_reverse in [True, False]):
        strand = '-' if align.is_reverse == True else '+'
    inf_strand = inferStrand(align, [x['feature'] for x in feature_list])
    if (strand is None) and (inf_strand is None):
        print 'Strand information could not be found nor inferred! Skipping alignment!'
        continue
    cleavage_site, clipped_pos, last_matched = inferCleavageSite(align, inf_strand)
    #read_count = 0
    #for overlapping_read in r2c.fetch(align.qname, align.qstart, align.qend):
        #read_count += 1
    for feat in feature_list:
        if (feat['feature'].strand == '+'):
            annotated_end = feat['feature'].end
        else:
            annotated_end = feat['feature'].start
        distance_from_end = annotated_end - cleavage_site
        feat['distance_from_end'] = distance_from_end
    print '\tcleavage site: {}'.format(cleavage_site)
    utr3 = None
    for feature in feature_list: 
        if (feature['feature'].feature == '3UTR'):
            if (utr3 is None):
                utr3 = feature
            elif (abs(feature['feature'].start - align.reference_start) < abs(utr3['feature'].start - align.reference_start)):
                utr3 = feature
        #print '\tFeature: {}\t{}\t{}\t{}\t{}'.format(feature['feature'].transcript_id,feature['feature'].feature,feature['feature'].start,feature['feature'].end, feature['feature'].strand)
    print '\tutr3: {}\t{}\t{}\t{}\t{}'.format(utr3['feature'].transcript_id, utr3['feature'].feature, utr3['feature'].strand, utr3['feature'].start, utr3['feature'].end)
    print '\t\tdistance from end: {}'.format(utr3['distance_from_end'])
    lines_link = ''
    link_pairs = find_link_pairs(align, chrom, clipped_pos, last_matched, cleavage_site, utr3)
    if link_pairs and args.output_reads:
        lines_link += output_link_pairs(align, link_pairs)
    print '-'*50
