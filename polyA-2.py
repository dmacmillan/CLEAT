import argparse
import logging
import pickle
import os
import sys
import re
import subprocess
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

def ucsc_chroms(genome):
    """Extracts conversion of UCSC chromosome names
    eg. hg19"""
    package_dir = "/".join(os.path.abspath(__file__).split("/")[:-2])
    conversion_file = package_dir + '/annotations/' + genome + "/ucsc_chr.txt"

    conversions = {}
    if os.path.exists(conversion_file):
        for line in open(conversion_file, 'r'):
            chr_from, chr_to = line.rstrip('\n').split()
            conversions[chr_from] = chr_to
            
    return conversions

# chrom_proper
chrom_proper = ucsc_chroms(args.ref_genome)

def qpos_to_tpos(align, qpos):
    blocks = align.blocks
    for i in xrange(len(blocks)):
        qb0 = blocks[i][0]
        qb1 = blocks[i][1]
        if ((qpos >= qb0) and (qpos <= qb1)) or ((qpos <= qb0) and (qpos >= qb1)):
            block = i
            break
    tpos = None
    try:
        if (align.is_reverse == False):
            tpos = blocks[block][0] + qpos - blocks[block][0]
        else:
            tpos = blocks[block][1] - (qpos - blocks[block][1])
    except NameError:
        pass
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
        
            freq = check_freq(subseq)
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
                         (align.query_sequence, cleavage_site, align.target, span[0], span[1]))
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
    #print '\tSkipped {} due to cigar > 1'.format(cigarskip)
    #print '\tSkipped {} due to anchor read strand == +'.format(anchor_rs_plus)
    #print '\tSkipped {} due to anchor read strand == -'.format(anchor_rs_minus)
    #print 'Done looking at reads'
    link_pairs = []
    for contig in mate_loc.keys():
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
        out += '>%s %s %s %s %s trimmed:%s\n%s\n' % (r[0].qname, align.qname, r[1].pos, anchor_direction, 
                                                     mate_contig, num_trimmed_bases, trimmed_seq) 
    return out

def annotate_cleavage_site(align, target, cleavage_site, clipped_pos, base, min_txt_match_percent=0.6):
    """Finds transcript where proposed cleavage site makes most sense, and also fetches matching ESTs

    This method assesses whether the cleavage site makes sense with
    a) the clipped position in the contig (start or end)
    b) the alignment strand
    c) the clipped base
    It determines what the strand the transcript should be on based on a) and b).
    Then it checks to see if the transcript strand agrees with c)
    (As an alternative, it can use splice sites, if the contig is large enough, to determine the 
    transcript strand, and then see whether both the clipped position and clipped base make sense)

    If the above makes sense, it moves to find transcript(s) that overlaps the alignment, given the
    'expected' strand and that cleavage site must exist 3' to the transcript's 3'UTR start, 
    using get_txts_with_min_match().
    If an annotated transcript whose end matches the cleavage site, that transcript is chosen.
    If not and more than one transcripts overlap, the one that is closest to the annotated 'end' is chosen.
    
    If a candidate transcript is chosen, it will try to find ESTs which end exactly at the 
    given cleavage site. (This is only done if the code is asked to at the command prompt)
    
    If no trancript can be found given the clipped postion, clipped base, and alignment, it will
    return None
    """     
    result = None
    
    chrom = proper_chrom(target, chrom_proper=chrom_proper)
    
    # determine which transcript strand can the cleavage site come from
    if clipped_pos == 'start':
        if align.is_reverse == False:
            txt_strand = '-'
        else:
            txt_strand = '+'
    else:
        if align.is_reverse == False:
            txt_strand = '+'
        else:
            txt_strand = '-'
            
    if (txt_strand == '+' and base == 'A') or\
       (txt_strand == '-' and base == 'T'):     
        ests = []
        
#        txts_screened = get_txts_with_min_match(align, cleavage_site, txt_strand)      
#        if txts_screened:
#            ## discard transcripts if proposed cleavage site is before 3'UTR
#            #txts_screened = [t for t in txts_screened if t['within_utr']]
#            #if txts_screened:
#            # if a transcript whose transcript end is exactly the same as proposed cleavage site,
#            # discard other candidates
#            exact_txts = [t for t in txts_screened if t['identical']]
#            if exact_txts:
#                txts_screened = exact_txts
#            
#            # pick transcript closest to proposed cleavage site
#            txts_screened.sort(key=lambda t: abs(t['distance_from_end']))
#            closest_txt = txts_screened[0]['txt']
#            
#            # find ESTs that end in exactly same site
#            if args.est_overlap:
#                if closest_txt.strand == '+':
#                    ests = extract_est(chrom, int(closest_txt.txStart), cleavage_site)
#                else:
#                    ests = extract_est(chrom, cleavage_site, int(closest_txt.txEnd))
#            else:
#                ests = None
#                
#            result = {
#                'ests': ests,
#                'txt': closest_txt,
#                'novel': not txts_screened[0]['identical'],
#                'within_utr': txts_screened[0]['within_utr'],
#                'coord': '%s:%d' % (align.target, cleavage_site),
#                'cleavage_site': cleavage_site,
#                'from_end': abs(txts_screened[0]['distance_from_end']),
#                'which_end': clipped_pos,
#                'base': base
#            }
#        else:
#            print "%s : %s:%s clipped is not 3' to CDS end of all matching transcripts" % (align.qname, target, cleavage_site)
#            
#    else:
#        print '%s : %s:%s clipped base(%s) does not match transcript strand(%s)' % (align.qname, target, cleavage_site, base, txt_strand)
        
    return result

def find_polyA_cleavage(align, target):
        """Finds PolyA cleavage sites of a given aligned contig
        
        This method first checks if the given contig captures a polyA tail (find_tail_contig),
        and then tries to capture bridge reads (find_bridge_reads) given the above result.
        These methods may return multiple polyA tails, which will then be assessed one-by-one
        against the annotated gene models to determine if each tail is plausible or not 
        (annotate_cleavage_site).  The final result is kept in a dictionary where the key is
        either the 'start' or the 'end' of the contig.
        If plausible tails are observed from both the 'start' and 'end' of the contig, the
        contig is dismissed.        
        """
        print 'checking', align.qname

        min_len = 1
        if filters and filters.has_key('min_at'):
            min_len = filters['min_at']
        mismatch = [1, 1]
        if filters and filters.has_key('max_diff'):
            mismatch = filters['max_diff']
        tail = find_bridge_reads(align, target, min_len, mismatch, tail=find_tail_contig(align, target, min_len, mismatch))
        
        results = []
        for clipped_pos in tail.keys():
            for event in tail[clipped_pos]:
                # contig coordinate of cleavage site                
                last_matched, cleavage_site, base, tail_seq, num_tail_reads, bridge_reads, bridge_clipped_seq = event
                result = annotate_cleavage_site(align, target, cleavage_site, clipped_pos, base)
                if result:
                    if bridge_reads:
                        result['num_bridge_reads'] = len(bridge_reads)
                        result['bridge_reads'] = bridge_reads
                        result['bridge_clipped_seq'] = bridge_clipped_seq
                    result['base'] = base
                    result['last_matched'] = last_matched
                    result['clipped_pos'] = clipped_pos
                    result['tail_seq'] = tail_seq
                    result['num_tail_reads'] = num_tail_reads
                    results.append(result)
        return results

def find_extended_bridge_reads(align, target, reads_to_screen, min_len, mismatch, genome_buffer=1000):
    """Finds bridge reads where only ending portion represent polyA tail"""
    query_seqs = {}
    for reads in reads_to_screen.values():
        for read in reads:
            query_seqs[read.qname] = read.seq
    
    #entirely_mapped = self.align_transcript_seq(align, query_seqs, 'extended', self.get_full_blat_aln)
    
    clipped_reads = {'start':{}, 'end':{}}
    for clipped_pos, reads in reads_to_screen.iteritems():
        if reads:
            if (clipped_pos == 'start' and align.is_reverse == False) or\
               (clipped_pos == 'end' and align.is_reverse == True):
                target_coord = [int(align.reference_start) - genome_buffer, int(align.reference_end)]
            else:
                target_coord = [int(align.reference_start), int(align.reference_end) + genome_buffer]
            
            query_seqs = dict((read.qname, read.seq) for read in reads) #if not entirely_mapped.has_key(read.qname))
            partial_aligns = align_genome_seq(align, target, query_seqs, target_coord, 'extended-bridge-genome', get_partial_blat_aln)
            
            read_objs = dict((read.qname, read) for read in reads)
            
            for read_name, mapped_coord in partial_aligns.iteritems():
                if mapped_coord[0] == 0:
                    clipped_seq = read_objs[read_name].seq[mapped_coord[1]:]
                else:
                    clipped_seq = read_objs[read_name].seq[:mapped_coord[0]]
                    
                if filters is not None and filters.has_key('min_bridge_size') and len(clipped_seq) < filters['min_bridge_size']:
                    continue
                                        
                # reverse complement to be in agreement with reference instead of contig
                clipped_seq_genome = clipped_seq
                if align.is_reverse == True:
                    clipped_seq_genome = revComp(clipped_seq)
                    
                if mapped_coord[0] == 0:
                    last_matched = read_objs[read_name].pos + mapped_coord[1]
                    pos_genome = target_coord[0] + mapped_coord[3] - 1
                else:
                    last_matched = read_objs[read_name].pos - mapped_coord[0]
                    pos_genome = target_coord[0] + mapped_coord[2]
                                        
                for base in ('A', 'T'):
                    if is_bridge_read_good(clipped_seq, base, min_len, mismatch):
                        #print 'possible extended bridge reads', align.query, read_objs[read_name].qname, read_objs[read_name].seq, is_seed, clipped_seq_genome
                        if not clipped_reads[clipped_pos].has_key(last_matched):
                            clipped_reads[clipped_pos][last_matched] = {}
                        if not clipped_reads[clipped_pos][last_matched].has_key(base):
                            clipped_reads[clipped_pos][last_matched][base] = []
                        clipped_reads[clipped_pos][last_matched][base].append([read_objs[read_name], clipped_seq_genome, pos_genome])
                        
    return clipped_reads

def merge_clipped_reads(clipped_reads, extended_clipped_reads):
    """Merges clipped_reads and extended_clipped_reads"""
    for clipped_pos in extended_clipped_reads.keys():
        for pos in extended_clipped_reads[clipped_pos].keys():
            if not clipped_reads[clipped_pos].has_key(pos):
                clipped_reads[clipped_pos][pos] = extended_clipped_reads[clipped_pos][pos]
            else:
                for base in extended_clipped_reads[clipped_pos][pos]:
                    if not clipped_reads[clipped_pos][pos].has_key(base):
                        clipped_reads[clipped_pos][pos][base] = extended_clipped_reads[clipped_pos][pos][base]
                    else:
                        clipped_reads[clipped_pos][pos][base].extend(extended_clipped_reads[clipped_pos][pos][base])

def filter_vs_reference(align, target, clipped_reads):
    """Filter bridge reads against reference genome sequence
    
    The reads are filtered against both the genomic region and the overlapping
    transcripts.
    This is done to check if bridge reads (including the clipped portion)
    can be mapped from end to end to a single location in the genome
    The reasons the reads were clipped could be because the contig was 
    reconstructed short of a polyA|T region.  The clipped sequences may be 
    aligned over an exon junction to the next exon (in which case checking 
    against the genome sequence won't help), or the clipped portion can 
    actually be aligned but there may be mismatches that cause the reads to
    be clipped by BWA-SW but not BLAT.
    """     
    # buffer on ends of alignment to extract reference genome sequence
    genome_buffer = 500
    same_reference = {}
    
    #for cleanup
    tmp_files = []
    
    if args.use_tmp:
        path = '/tmp'
    else:
        path = os.path.dirname(os.path.abspath(args.out))
        
    # create genome sequence file for Blatting
    target_file = '%s/%s-genome.fa' % (path, align.qname)
    if os.path.exists(target_file):
        os.remove(target_file)
    out = open(target_file, 'w')
    genome_start, genome_end = int(align.reference_start) - genome_buffer, int(align.reference_end) + genome_buffer
    seq = refseq.fetch(target, genome_start, genome_end)
    out.write('>%s:%s-%s\n%s\n' % (target, genome_start, genome_end, seq))
    out.close()
    tmp_files.append(target_file)
    

    # create reads for Blatting
    reads_file = '%s/%s-bridge.fa' % (path, align.qname)
    if os.path.exists(reads_file):
        os.remove(reads_file)
    out = open(reads_file, 'w')
    tmp_files.append(reads_file)
    
    for clipped_pos in clipped_reads.keys():
        for pos in clipped_reads[clipped_pos].keys():               
            for base in clipped_reads[clipped_pos][pos]:        
                for read in clipped_reads[clipped_pos][pos][base]:
                    out.write('>%s\n%s\n' % (read[0].qname, read[0].seq))
    out.close()

    # align reads against genome
    fully_aligned_genome = None
    aln_file = '%s/%s-bridge-genome.psl' % (path, align.qname)
    if os.path.exists(aln_file):
        os.remove(aln_file)
    tmp_files.append(aln_file)
    try:
        subprocess.call(['blat', target_file, reads_file, aln_file])
    except CalledProcessError as err:
        sys.stderr.write('error running blat:%s' % err.cmd)
    else:
        if os.path.exists(aln_file):
            fully_aligned_genome = get_full_blat_aln(aln_file)
            
    # create transcript sequence file for Blatting
    target_file = '%s/%s-txt.fa' % (path, align.qname)
    if os.path.exists(target_file):
        os.remove(target_file)
    target_empty = extract_transcript_seq(align, target, target_file)
    tmp_files.append(target_file)
    
    # align reads against transcripts
    fully_aligned_txt = None
    aln_file = '%s/%s-bridge-txt.psl' % (path, align.qname)
    if os.path.exists(aln_file):
        os.remove(aln_file)
    tmp_files.append(aln_file)
    if not target_empty:                
        try:
            subprocess.call(['blat', target_file, reads_file, aln_file])
        except CalledProcessError as err:
            sys.stderr.write('error running blat:%s' % err.cmd)
        else:
            if os.path.exists(aln_file):
                fully_aligned_txt = get_full_blat_aln(aln_file)
                                    
    # filtering
    for clipped_pos in clipped_reads.keys():
        for pos in clipped_reads[clipped_pos].keys():               
            for base in clipped_reads[clipped_pos][pos]:        
                # check if clipped sequence is genomic/transcriptomic, otherwise skip
                is_genomic = is_transcriptome = bad_neighbor = False
                for read in clipped_reads[clipped_pos][pos][base]:
                    if fully_aligned_genome is not None and fully_aligned_genome.has_key(read[0].qname):
                        is_genomic = True
                        break
                    
                    if fully_aligned_txt is not None and fully_aligned_txt.has_key(read[0].qname):
                        is_transcriptome = True
                        break
                        
                    if read[2] is not None and in_homopolymer_neighbor(target, read[2], read[1], base):
                        bad_neighbor = True
                        break
                                        
                if is_genomic or is_transcriptome or bad_neighbor:
                    if is_genomic:
                        reason = 'genome seq'
                    elif is_transcriptome:
                        reason = 'transcriptome seq'
                    else:
                        reason = 'neighborhood'
                        
                    print 'failed %s check %s %s:%s-%s %s %s %s %s %s' % (reason, align.qname, target, align.reference_start, align.reference_end, pos, read[2], read[0].qname, read[0].seq, read[1])
                    
                    if not same_reference.has_key(clipped_pos):
                        same_reference[clipped_pos] = []
                    same_reference[clipped_pos].append([pos, base])
                        
                    
    # remove entries that are same as reference
    for clipped_pos in same_reference.keys():
        for (pos, base) in same_reference[clipped_pos]:
            if clipped_reads[clipped_pos].has_key(pos) and clipped_reads[clipped_pos][pos].has_key(base):
                del clipped_reads[clipped_pos][pos][base]
            
    # clean up temporary alignment files
    for ff in tmp_files:
        if os.path.exists(ff):
            #print 'cleanup', ff
            os.remove(ff)

def find_bridge_reads(align, target, min_len, mismatch, genome_buffer=1000, tail=None):
    """Finds bridge reads
    
    It first checks to see clipped reads if the entire clipped portion is A's or T's.
    Then it will check, through find_exteneded_bridge_reads() if the ending portion of 
    the clipped sequence is A's or T's.
    The 2 results will be merged together.
    """
    # used for check if read is mapped to the aligned portion of the contig
    query_bounds = sorted([int(align.qstart), int(align.qend)])
    
    # identify clipped reads that are potential pA/pT
    clipped_reads = {'start':{}, 'end':{}}
    second_round = {'start':[], 'end':[]}
    for read in r2c.fetch(align.qname):
    #for read in self.bam.bam.fetch(align.query):
        if not read.cigar or len(read.cigar) != 2:
            continue
                    
        if (read.cigar[0][0] == 4 or read.cigar[0][0] == 5) or\
           (read.cigar[-1][0] == 4 or read.cigar[-1][0] == 5):
            # clipped at start
            if read.cigar[0][0] == 4 or read.cigar[0][0] == 5:
                clipped_pos = 'start'
                last_matched = read.pos + 1
                # to see whether the clipped seq is a pA tail
                clipped_seq = read.seq[:read.cigar[0][1]]
                # for trimming of poor quality bases if so descired
                clipped_qual = read.qual[:read.cigar[0][1]]
                
            # clipped at end
            else:
                clipped_pos = 'end'
                last_matched = read.pos + read.alen
                # to see whether the clipped seq is a pA tail
                clipped_seq = read.seq[-1 * read.cigar[-1][1]:]
                # for trimming of poor quality bases if so descired
                clipped_qual = read.qual[-1 * read.cigar[-1][1]:]
                
            # if last_match is beyond the limit of the alignment, adjust last_matched
            if last_matched < query_bounds[0] or last_matched > query_bounds[1]:
                if last_matched < query_bounds[0]:
                    diff = query_bounds[0] - last_matched
                else:
                    diff = last_matched - query_bounds[1]
                if clipped_pos == 'start':
                    last_matched = last_matched + diff
                    clipped_seq = read.seq[:read.cigar[0][1] + diff]
                    clipped_qual = read.qual[:read.cigar[0][1] + diff]
                else:
                    last_matched = last_matched - diff
                    clipped_seq = read.seq[-1 * (read.cigar[-1][1] + diff):]
                    clipped_qual = read.qual[-1 * (read.cigar[-1][1] + diff):]
                                    
            # trim poor quality base if desired
            if args.trim_reads:
                clipped_seq = trim_bases(clipped_seq, clipped_qual, clipped_pos)
                
            if len(clipped_seq) < 1:
                continue
            if filters is not None and filters.has_key('min_bridge_size') and len(clipped_seq) < filters['min_bridge_size']:
                continue
            
            # reverse complement to be in agreement with reference instead of contig
            clipped_seq_genome = clipped_seq
            if align.is_reverse == True:
                clipped_seq_genome = revComp(clipped_seq)
            
            # check for possible tail (stretch of A's or T's)
            #pos_genome = align.qpos_to_tpos(last_matched)
            pos_genome = qpos_to_tpos(align, last_matched)
            picked = False
            for base in ('A', 'T'):
                if is_bridge_read_good(clipped_seq_genome, base, min_len, mismatch):
                    if not clipped_reads[clipped_pos].has_key(last_matched):
                        clipped_reads[clipped_pos][last_matched] = {}
                    if not clipped_reads[clipped_pos][last_matched].has_key(base):
                        clipped_reads[clipped_pos][last_matched][base] = []
                    clipped_reads[clipped_pos][last_matched][base].append([read, clipped_seq_genome, pos_genome])
                    picked = True
                    
            if not picked:
                second_round[clipped_pos].append(read)
                    
    extended_clipped_reads = find_extended_bridge_reads(align, target, second_round, min_len, mismatch)
    merge_clipped_reads(clipped_reads, extended_clipped_reads)
                
    # filter events against reference sequence
    filter_vs_reference(align, target, clipped_reads)
    
    # translate into results
    if tail is None:
        results = {}
        tail_search = None
    else:
        results = tail
        tail_search = {'start':{}, 'end':{}}
        for clipped_pos in tail.keys():
            for i in range(len(tail[clipped_pos])):
                cleavage_site = tail[clipped_pos][i][1]
                tail_search[clipped_pos][cleavage_site] = i
                            
    for clipped_pos in clipped_reads.keys():
        if not results.has_key(clipped_pos):
            results[clipped_pos] = []
        for pos in clipped_reads[clipped_pos].keys():                       
            for base in clipped_reads[clipped_pos][pos]:
                pos_genome = clipped_reads[clipped_pos][pos][base][0][2]
                
                if pos_genome is not None:
                    if tail_search is not None and tail_search[clipped_pos].has_key(pos_genome):
                        results_idx = tail_search[clipped_pos][pos_genome]
                        results[clipped_pos][results_idx][5] = [r[0] for r in clipped_reads[clipped_pos][pos][base]]
                        results[clipped_pos][results_idx][6] = [r[1] for r in clipped_reads[clipped_pos][pos][base]]
                    else:
                        results[clipped_pos].append([pos,
                                                     pos_genome, 
                                                     base,
                                                     None, 
                                                     None,
                                                     [r[0] for r in clipped_reads[clipped_pos][pos][base]], 
                                                     [r[1] for r in clipped_reads[clipped_pos][pos][base]]
                                                     ])      
    return results

def proper_chrom(chrom, genome=None, chrom_proper=None):
    """Returns proper chromosome name
    UCSC-format if available
    """
    if not chrom_proper and genome:
        chrom_proper = ucsc_chroms(genome)
    
    if chrom_proper:
        if chrom_proper.has_key(chrom):
            chrom = chrom_proper[chrom]
        elif chrom[:3].lower() == 'chr' and chrom_proper.has_key(chrom[3:]):
            chrom = chrom_proper[chrom[3:]]
    
    if not re.match('^(chr|scaffold)', chrom, re.IGNORECASE):
        chrom = 'chr' + chrom
            
    return chrom

def in_homopolymer_neighbor(chrom, pos, tail, base):
    """Checks if tail is juxtaposed against a homopolymer of the given base in genome
    
    The window size for checking is 5 or 2*length of the tail sequence (tail), whichever larger.
    The homopolyer run is constructed by making a string of the given base (base) with length 
    equal to the frequency of that base in the tail sequence.
    If the homopolyer is found immedicately before or after the given position (pos),
    then the result is True.
    """
    min_len = 5
    length = max(min_len, 2 * len(tail))
    neighbor_seq = refseq.fetch(chrom, pos-length, pos+length)
    #neighbor_seq = self.refseq.GetSequence(chrom, pos-length, pos+length)
    homo = ''
    freq = check_freq(tail)
    
    if freq.has_key(base):
        for i in range(length):
            homo += base.upper()
                
        m = re.search(homo, neighbor_seq, re.IGNORECASE)
        if m:               
            # if homopolymer is touching the middle base of neighbor_seq
            if m.start() <= len(neighbor_seq)/2 + 1 and m.end() - 1 >= len(neighbor_seq)/2 - 1:
                return True
                    
    return False

def get_num_tail_reads(align, last_matched):
    """Reports number of reads spanning cleavage site in contig"""
    num = 0
    for read in r2c.fetch(align.qname):
        if not read.cigar or len(read.cigar) != 1:
            continue
        
        if read.pos + 1 <= last_matched and read.pos + read.alen > last_matched:
            num += 1
    
    return num

def find_tail_contig(align, target, min_len, mismatch):
    """Finds contigs that have polyA tail reconstructed"""
    # size of stretch of contig sequence to be added to polyA tail
    # against trasncript sequences to see it polyA tail is genomic
    junction_buffer=50
    results = {}
    
    clipped = {'start':False, 'end':False}
    if int(align.qstart) > 1:
        clipped['start'] = True
        
    if int(align.qend) < int(align.query_length):
        clipped['end'] = True
        
    for clipped_pos in ('start', 'end'):
        if clipped[clipped_pos]:
            if clipped_pos == 'start':
                last_matched = int(align.qstart)
                clipped_seq = align.query_sequence[:int(align.qstart)-1]
                junction_seq = align.query_sequence[:int(align.qstart) - 1 + junction_buffer]
            else:
                last_matched = int(align.qend)
                clipped_seq = align.query_sequence[int(align.qend):]
                junction_seq = align.query_sequence[int(align.qend) - junction_buffer:]
                                    
            cleavage_site = qpos_to_tpos(align, last_matched)
            clipped_seq_genome = clipped_seq
            if align.is_reverse == True:
                clipped_seq_genome = revComp(clipped_seq)
            
            matched_transcript = in_homopolymer = False
            for base in ('A', 'T'):
                if matched_transcript or in_homopolymer:
                    continue
                                        
                perfect = is_polyA_tail(clipped_seq_genome, base, min_len=1, max_nonAT_allowed=[0, 1])
                #imperfect = self.is_polyA_tail(clipped_seq_genome, base, min_len=4, max_nonAT_allowed=[1, 4])
                imperfect = is_polyA_tail(clipped_seq_genome, base, min_len=min_len, max_nonAT_allowed=mismatch)
                
                if perfect or imperfect:
                    # don't need to do the following 2 checks if it's not a potential tail
                    if len(clipped_seq) == 1 and in_homopolymer_neighbor(target, cleavage_site, clipped_seq_genome, clipped_seq[0]):
                        print '%s : clipped seq in middle of homopolyer run (%s) %s' % (align.qname, clipped_pos, clipped_seq)
                        in_homopolymer = True
                        continue
                    
                    #entirely_mapped = align_transcript_seq(align, target, {align.qname:junction_seq}, 'junction', get_full_blat_aln)
                    #if entirely_mapped.has_key(align.qname):
                    #    print '%s : clipped seq is just transcript seq (%s) %s %d' % (align.qname, clipped_pos, clipped_seq, len(clipped_seq))
                    #    matched_transcript = True
                    #    continue
                    
                    # find reads corresponding to tail
                    num_tail_reads = get_num_tail_reads(align, last_matched)
                                                                    
                    if not results.has_key(clipped_pos):
                        results[clipped_pos] = []
                    results[clipped_pos].append([last_matched,
                                                 cleavage_site,
                                                 base,
                                                 clipped_seq_genome,
                                                 num_tail_reads,
                                                 None,
                                                 None,
                                                 ]
                                                )
            
                                        
    if not clipped['start'] and not clipped['end']:
        results = None
        
    return results

def align_transcript_seq(align, target, query_seqs, label, parse_fn):
    """Aligns(BLAT) query sequences to transcripts overlapping alignment
    
    Query sequences is given in a hash (query_seq) where
    key=query name, value=query sequence
    Targets are all transcript sequences overlapping 'tstart' and 'tend'
    of the alignment.
    'label' will be added in addition to query name to all the 
    temporary files
    The BLAT alignments will be processed by parse_fn() to return the results
    """
    result = None
    
    #for cleanup
    tmp_files = []
    
    if args.use_tmp:
        path = '/tmp'
    else:
        path = os.path.dirname(os.path.abspath(args.out))
    
    # create transcript sequence file for Blatting
    target_file = '%s/%s-target-%s.fa' % (path, align.qname, label)
    if os.path.exists(target_file):
        os.remove(target_file)
    target_empty = extract_transcript_seq(align, target, target_file)
    tmp_files.append(target_file)

    # create query for Blatting
    query_file = '%s/%s-query-%s.fa' % (path, align.qname, label)
    if os.path.exists(query_file):
        os.remove(query_file)
    out = open(query_file, 'w')
    for query, seq in query_seqs.iteritems():
        out.write('>%s\n%s\n' % (query, seq))
    out.close()
    tmp_files.append(query_file)

    # align query against target
    aln_file = '%s/%s-%s.psl' % (path, align.qname, label)
    if os.path.exists(aln_file):
        os.remove(aln_file)
    tmp_files.append(aln_file)
    try:
        subprocess.call(['blat', target_file, query_file, aln_file])
    except CalledProcessError as err:
        sys.stderr.write('error running blat:%s' % err.cmd)
    else:
        if os.path.exists(aln_file):
            result = parse_fn(aln_file)
            
    # clean up temporary alignment files
    for ff in tmp_files:
        if os.path.exists(ff):
            #print 'cleanup', ff
            os.remove(ff)
            
    return result

def get_full_blat_aln(aln_file):
    """Extracts full hits from BLAT aligments
    
    This is used for removing false-positive bridge reads where their entirety in
    sequence can be mapped to a single transcript
    """
    fully_aligned = {}
    for line in open(aln_file, 'r'):
        if not re.search('^\d', line):
            continue
        cols = line.rstrip('\n').split('\t')
        query, qsize, qstart, qend, target = cols[9:14]
        block_count = cols[17]
        #print 'aln', query, qsize, qstart, qend, block_count, target
        
        if int(qstart) == 0 and int(qsize) == int(qend) and int(block_count) == 1:
            fully_aligned[query] = True
        
    return fully_aligned

def extract_transcript_seq(align, target, out_file):
    """Extracts transcripts overlapping an alignment and outputs their
    sequences to the given file
    
    Returns empty=True if it fails to find overlapping transcripts
    """
    
    chrom = proper_chrom(target, chrom_proper=chrom_proper)
    feats = features.fetch(chrom, align.reference_start, align.reference_end)
    out = open(out_file, 'w')
    empty = True
    transcripts = {}
    for feature in feats:
        if (feature.transcript_id not in transcripts):
            transcripts[feature.transcript_id] = ''
        if feature.feature == 'exon':
            empty = False
            transcripts[feature.transcript_id] += refseq.fetch(chrom, feature.start, feature.end).upper()
    for tid in transcripts:
        print "transcript_seq: {}".format(transcripts[tid])
        out.write('>%s\n%s\n' % (tid, transcripts[tid]))
    out.close()
    
    return empty

def get_partial_blat_aln(aln_file):
    """Extracts single-block, partial hits from BLAT aligments
    
    This is for capturing the polyA tails of extended bridge reads/
    The clipped portion of extended bridge reads contains both genomic sequence
    and the polyA tail. By alignment these reads to the genmomic region, the 
    polyA tail should be unaligned whereas the genomic portion would align.
    The return variable is a dictionary, where:
    key = query(read) name
    value = [qstart, qend, tstart, tend]
    """

    partially_aligned = {}
    for line in open(aln_file, 'r'):
        if not re.search('^\d', line):
            continue
        cols = line.rstrip('\n').split('\t')
        query, qsize, qstart, qend, target, tsize, tstart, tend = cols[9:17]
        block_count = cols[17]
        
        if int(block_count) == 1 and \
            ((int(qstart) == 0 and int(qend) < int(qsize)) or\
             (int(qstart) > 0 and int(qsize) == int(qend))):
            partially_aligned[query] = [int(qstart), int(qend), int(tstart), int(tend)]
        
    return partially_aligned

def align_genome_seq(align, target, query_seqs, coord, label, parse_fn):
    """Aligns(BLAT) query sequences to genomic sequence of given coordinate
    
    Query sequences is given in a hash (query_seq) where
    key=query name, value=query sequence
    Target is genomic sequence between coord[0] and coord[1]
    'label' will be added in addition to query name to all the 
    temporary files
    The BLAT alignments will be processed by parse_fn() to return the results
    """
    result = None
    
    #target_seq = self.refseq.GetSequence(align.target, coord[0], coord[1])
    target_seq = refseq.fetch(target, coord[0], coord[1])
    
    #for cleanup
    tmp_files = []
    
    if args.use_tmp:
        path = '/tmp'
    else:
        path = os.path.dirname(os.path.abspath(args.out))
    
    # create genome reference file for Blatting
    target_file = '%s/%s-genome-%s.fa' % (path, align.qname, label)
    if os.path.exists(target_file):
        os.remove(target_file)
    out = open(target_file, 'w')
    out.write('>%s:%d-%d\n%s\n' % (target, coord[0], coord[1], target_seq))
    out.close()
    tmp_files.append(target_file)

    # create query for Blatting
    query_file = '%s/%s-query-%s.fa' % (path, align.qname, label)
    if os.path.exists(query_file):
        os.remove(query_file)
    out = open(query_file, 'w')
    for query, seq in query_seqs.iteritems():
        out.write('>%s\n%s\n' % (query, seq))
    out.close()
    tmp_files.append(query_file)

    # align query against target
    aln_file = '%s/%s-%s.psl' % (path, align.qname, label)
    if os.path.exists(aln_file):
        os.remove(aln_file)
    tmp_files.append(aln_file)
    try:
        subprocess.call(['blat', target_file, query_file, aln_file])
    except CalledProcessError as err:
        sys.stderr.write('error running blat:%s' % err.cmd)
    else:
        if os.path.exists(aln_file):
            result = parse_fn(aln_file)
            
    # clean up temporary alignment files
    for ff in tmp_files:
        if os.path.exists(ff):
            #print 'cleanup', ff
            os.remove(ff)
            
    return result

def output_result(self, align, result, link_pairs=[]):
    """Outputs main results in tab-delimited format
    
    The output fields are specified in class variable 'output_fields'
    """
    data = {}
    data['contig'] = align.query
    
    data['transcript'] = result['txt'].name
    data['transcript_strand'] = result['txt'].strand
    data['chromosome'] = result['txt'].chrom
    if result['txt'].alias:
        data['gene'] = result['txt'].alias
        
    coding_type = result['txt'].coding_type()
    if coding_type == 'CODING':
        data['coding'] = 'yes'
    elif coding_type == 'NONCODING':
        data['coding'] = 'no'
    else:
        data['coding'] = 'unknown'
        
    if result['within_utr']:
        data['within_UTR'] = 'yes'
    else:
        data['within_UTR'] = 'no'
        
    data['chrom'] = result['txt'].chrom
    data['cleavage_site'] = result['cleavage_site']
    data['distance_from_annotated_site'] = result['from_end']
    
    if result.has_key('tail_seq'):
        if result['tail_seq'] is None:
            data['length_of_tail_in_contig'] = 0
            data['number_of_tail_reads'] = 0
        else:
            data['length_of_tail_in_contig'] = len(result['tail_seq'])
            data['number_of_tail_reads'] = result['num_tail_reads']
    else:
        data['length_of_tail_in_contig'] = data['number_of_tail_reads'] = '-'
        
    if result.has_key('bridge_reads') and result['bridge_reads']:
        data['number_of_bridge_reads'] = len(result['bridge_reads'])
        data['max_bridge_read_tail_length'] = max([len(s) for s in result['bridge_clipped_seq']])
        data['bridge_read_identities'] = ','.join([read.qname for read in result['bridge_reads']])
    else:
        data['number_of_bridge_reads'] = data['max_bridge_read_tail_length'] = 0
        data['bridge_read_identities'] = '-'
        
    if link_pairs is not None:
        if link_pairs:
            data['number_of_link_pairs'] = len(link_pairs)
            data['link_pair_identities'] = ','.join(r[0].qname for r in link_pairs)
            data['max_link_pair_length'] = max([len(r[-1]) for r in link_pairs])
        else:
            data['number_of_link_pairs'] = data['max_link_pair_length'] = 0
            
    elif not self.no_link:
        data['number_of_link_pairs'] = data['max_link_pair_length'] = 0
            
    if result['ests'] is not None:
        data['ESTs'] = '%s' % (len(result['ests']))
    else:
        data['ESTs'] = '-'
            
    data['tail+bridge_reads'] = 0
    if data['number_of_tail_reads'] != '-':
        data['tail+bridge_reads'] += data['number_of_tail_reads']
    if data['number_of_bridge_reads'] != '-':
        data['tail+bridge_reads'] += data['number_of_bridge_reads']
        
    cols = []
    for field in self.output_fields:
        if data.has_key(field):
            cols.append(str(data[field]))
        else:
            cols.append('-')
            
    result = '%s\n' % '\t'.join(cols)
    return result


# Iterate contig to genome alignments
for align in aligns:
    print 'Looking at contig: {}: {}-{}'.format(align.qname,align.reference_start, align.reference_end)
    # Get chromosome
    chrom = aligns.getrname(align.tid)
    #if (align.qname == 'k52.R224263'):
    #    extract_transcript_seq(align, chrom, './'+align.qname+'.look')
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
    for feat in feature_list:
        if (feat['feature'].strand == '+'):
            annotated_end = feat['feature'].end
        else:
            annotated_end = feat['feature'].start
        distance_from_end = annotated_end - cleavage_site
        feat['distance_from_end'] = distance_from_end
    #print '\tcleavage site: {}'.format(cleavage_site)
    utr3 = None
    for feature in feature_list: 
        if (feature['feature'].feature == '3UTR'):
            if (utr3 is None):
                utr3 = feature
            elif (abs(feature['feature'].start - align.reference_start) < abs(utr3['feature'].start - align.reference_start)):
                utr3 = feature
        #print '\tFeature: {}\t{}\t{}\t{}\t{}'.format(feature['feature'].transcript_id,feature['feature'].feature,feature['feature'].start,feature['feature'].end, feature['feature'].strand)
    #print '\tutr3: {}\t{}\t{}\t{}\t{}'.format(utr3['feature'].transcript_id, utr3['feature'].feature, utr3['feature'].strand, utr3['feature'].start, utr3['feature'].end)
    #print '\t\tdistance from end: {}'.format(utr3['distance_from_end'])
    lines_link = ''
    link_pairs = find_link_pairs(align, chrom, clipped_pos, last_matched, cleavage_site, utr3)
    if link_pairs and args.output_reads:
        lines_link += output_link_pairs(align, link_pairs)
    results = find_polyA_cleavage(align,chrom)
    print '-'*50
    if results:
        for result in results:
            lines_result += output_result(align, result, link_pairs=link_pairs)
            
            if options.output_reads and result.has_key('bridge_reads') and result['bridge_reads']:
                lines_bridge += output_bridge_reads(align, result)
                
    elif result_link is not None:
        lines_result += output_result(align, result_link, link_pairs=link_pairs)
              
# close output streams
pf.output(lines_result, lines_bridge, lines_link)
