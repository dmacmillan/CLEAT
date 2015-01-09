import argparse
import logging
import os
import sys
import re
import subprocess
# External modules below
import pysam

parser = argparse.ArgumentParser(description='this program tries to find polya cleavage sites through short-read assembly.it is expected that contigs are aligned to contigs, and reads aligned to contigs. these 2 alignment steps can be performed by trans-abyss. the aligners used are gmap for contig-genome and bwa-sw for read-contig alignments. annotations files for ensembl, knowngenes, refseq, and aceview are downloaded from ucsc. est data(optional) are also downloaded from ucsc. the analysis can be composed of 2 phases: 1. contig-centric phase - cleavage sites per contig are captured 2. coordinate-centric phase - contigs capturing the same cleavage site are consolidated into 1 report where expression/evidence-related data are summed. customized filtering based on evidence data can be performed.')
parser.add_argument('c2g', metavar='<contig-to-genome>', help='The contig-to-genome alignment file in bam format')
parser.add_argument('contigs', metavar='<contigs>', help='The file containing all contigs in fasta format')
parser.add_argument('ref_genome', metavar='<reference_genome>', help='The reference genome to use. default is hg19.', default='hg19', choices=['hg18','hg19'])
parser.add_argument('annot', metavar='<annotations>', help='The annotations file to use with the reference in gtf format.')
parser.add_argument('r2c', metavar='<reads-to-contigs>', help='The contigs-to-genome alignment file')
parser.add_argument('out', metavar='<output-file>', help='The file to output results')
parser.add_argument('-r', '--output_reads', dest='output_reads', help='Enable to output read sequences', action='store_true', default=False)
parser.add_argument('-t', '--trim', dest='trim_reads', help='Trim bases of quality <= this value. default is 3.', type=int, default=3)
parser.add_argument('-ss', '--strand_specific', action='store_true', help='Enable if reads are strand specific.')
parser.add_argument('--use_tmp', help='Use tmp space', action='store_true', default=False)
parser.add_argument('--no_link', help='Do not find link pairs', action='store_true', default=False)
parser.add_argument('-e', '--overlap_est', help='Overlap expressed sequence tags', action='store_true', default=False)
parser.add_argument('--min_at', help='Minimum number of a|t bases in tail. default is 4.', type=int, default=4)
parser.add_argument('--max_diff', help='Maximum rate of x to y bases allowed in tail. where x are non a|t bases and y are a|t bases. default: 1 5.', nargs=2, default=[1,5])
parser.add_argument('--max_diff_link', help='Maximum number of non A|T bases in entire link read. Default is 2.', type=int, default=2)
parser.add_argument('--min_bridge_size', help='Minimum size of bridge. Default is 1.', type=int, default=1)
parser.add_argument('-k', '--track', metavar=('[name]','[description]'), help='Name and description of BED graph track to output.', nargs=2)
parser.add_argument('--rgb', help='RGB value of BED graph. Default is 0,0,255', default='0,0,255')
parser.add_argument('-c', help='Specify a contig to look at.')
parser.add_argument('-x', '--combine', dest='combine', action='store_true', help='Path of discovery phase results for combining.')

args = parser.parse_args()
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('polyA_logger')

#fh = logging.FileHandler('info.log')
#fh.setLevel(logging.INFO)
#logger.addHandler(fh)

# Globally used variables
bridge_total = 0
link_pairs_total = 0
basic_total = 0


# Reference genome sequence (refseq)
if (args.ref_genome == 'hg19'):
    reference_genome = '/projects/btl/trans-abyss/public_releases/v1.4.8/annotations/hg19/hg19.fa'
logger.debug("Loading reference genome via pysam 0.8.1")
refseq = pysam.FastaFile(reference_genome)
logger.debug("Reference genome successfully loaded!")

# Contigs to genome alignments (aligns)
logger.debug("Checking contig to genome alignment is in BAM format. Splitting file by '.': %s", args.c2g.split('.'))
if (args.c2g.split('.')[-1].lower() == 'bam'):
    logger.debug("Contig to genome alignment file is in BAM format, loading BAM file...")
    aligns = pysam.AlignmentFile(args.c2g, "rb")
    logger.debug("Contig to genome alignment file loaded successfully")
elif (args.c2g.split('.')[-1].lower() == 'sam'):
    aligns = pysam.AlignmentFile(args.c2g, "r")
else:
    sys.exit("Unrecognized format for <aln>, must be BAM format file. Exiting.")

# Contig sequences
contigs = pysam.FastaFile(args.contigs)

# Features (features)
features = pysam.TabixFile(args.annot, parser=pysam.asGTF())

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

output_fields=['gene','transcript','transcript_strand','coding','contig','chromosome','cleavage_site','within_UTR','distance_from_annotated_site','ESTs','length_of_tail_in_contig','number_of_tail_reads','number_of_bridge_reads','max_bridge_read_tail_length','bridge_read_identities','tail+bridge_reads','number_of_link_pairs','max_link_pair_length','link_pair_identities']
if args.c:
    args.out += '.'+args.c
out_result = open(args.out, 'a')
out_result.write('%s\n' % '\t'.join(output_fields))
prefix = os.path.splitext(args.out)[0]
out_bridge_reads_file = prefix + '-reads.fa'
out_link_pairs_file = prefix + '-link.fa'
out_bridge_reads = open(out_bridge_reads_file, 'a')
out_link_pairs = open(out_link_pairs_file, 'a')

# Filters (filters)
filters = {}
filters['min_at'] = args.min_at
filters['max_diff'] = args.max_diff
filters['max_diff_link'] = args.max_diff_link
if args.combine:
    filters['min_bridge_size'] = args.min_bridge_size

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

def group_and_filter(lines_result, out_file, filters=None, make_track=None, rgb='0,0,0'):
    global output_fields
    """Consolidates contig-centric results into coordinate-centric results
    
    path = directory of tab-delimited results files
    """
    groups = {}
    
    lines_result = [x.split('\t') for x in lines_result.splitlines()]
    for line in lines_result:
        if line[0] == output_fields[0]:
            continue
        
        chrom, cleavage_site = line[5], line[6]
        if not groups.has_key(chrom):
            groups[chrom] = {}
        if not groups[chrom].has_key(cleavage_site):
            groups[chrom][cleavage_site] = []
        groups[chrom][cleavage_site].append(line)
            
    stats = {'gene': {},
             'transcript': {'coding': {}, 'noncoding':{}},
             'screened_out': 0,
             'cleavage_site': 0,
             'with_tail': 0,
             'tail_and_bridge_and_link': 0,
             'tail_and_bridge': 0,
             'tail_and_link': 0,
             'bridge_and_link': 0,
             'just_tail': 0,
             'just_bridge': 0,
             'just_link': 0,
             }
            
    out = open(out_file, 'w')
    out.write('%s\n' % '\t'.join(output_fields))
    
    track_plus = []
    track_minus = []
    for chrom in sorted(groups.keys(), cmp=compare_chr):
        for cleavage in sorted(groups[chrom].keys()):
            # Skip if there is no tail, bridge reads, and link pairs
            if (int(result[10]) == 0) and (int(result[12]) == 0) and (int(result[16]) == 0):
                continue
            results = groups[chrom][cleavage]
            # if more than one contig reports same cleavage site, add up the support numbers
            if len(results) > 1:
                result = merge_results(results)
            else:
                result = results[0]
                                    
            if filters is not None and filters.has_key('min_bridge_size') and int(result[13]) < filters['min_bridge_size']:
                continue
                
            out.write('%s\n' % '\t'.join(result))
            if make_track is not None:
                if result[2] == '+':
                    track_plus.append(show_expression(result))
                else:
                    track_minus.append(show_expression(result))
            
            # stats
            update_stats(stats, result)
                
    out.close()
        
    # prefix used for track and stats files
    prefix = os.path.splitext(out_file)[0]
    # output track
    if make_track is not None:
        track_file = prefix + '.bg'
        cls.output_track(track_file, make_track[0], make_track[1], rgb=rgb, plus=track_plus, minus=track_minus)
        
    # output stats file
    stats_file = prefix + '.stats'
    cls.output_stats(stats, stats_file)

def output_track(out_file, name, desc, rgb='0,0,0', plus=[], minus=[]):
    """Outputs bedgraph, separate tracks for plus and minus strands"""
    out = open(out_file, 'w')
    # plus strand
    out.write('%s\n' % prepare_track_header(name + ' (+)', desc + ' plus strand', rgb))
    for line in plus:
        out.write('%s\n' % line)
    # minus strand
    out.write('%s\n' % prepare_track_header(name + '(-)', desc + ' minus strand', rgb))
    for line in minus:
        out.write('%s\n' % line)
    out.close()

def merge_results(results):        
    """Merges results from different contigs of same cleavage site into single result"""
    # join fields: contig, bridge_name, link_name
    join = [4, 14, 18]              
    # add fields: num_tail_reads, num_bridge_reads, tail+bridge, num_link_pairs
    add = [11, 12, 15, 16]
    # max fields: tail_len, bridge_len, link_len
    biggest = [10, 13, 17]
    
    merged = []
    for i in range(len(results[0])):
        if i in add:
            data = [int(r[i]) for r in results if r[i].isdigit()]
            if data:
                merged.append(str(sum(data)))
            else:
                merged.append('-')              
        elif i in biggest:
            data = [int(r[i]) for r in results if r[i].isdigit()]
            if data:
                merged.append(str(max(data)))
            else:
                merged.append('-')              
        elif i in join:
            data = [r[i] for r in results if r[i] != '-']
            
            # if there is data other than '-', then list all items
            if data:
                merged.append(','.join([r[i] for r in results if r[i] != '-']))
            else:
                merged.append('-')                  
        else:
            merged.append(results[0][i])

    return merged

def update_stats(stats, result):
    """Updates summary stats with result
    
    stats = dictionary of final stats results
    result = list of values of each output line
    """
    has_tail = has_bridge = has_link = False
    if result[10].isdigit() and int(result[10]) > 0:
        has_tail = True
    if result[12].isdigit() and int(result[12]) > 0:
        has_bridge = True
    if result[16].isdigit() and int(result[16]) > 0:
        has_link = True
        
    if not (has_tail or has_bridge or has_link):
        return
        
    stats['cleavage_site'] += 1
    
    if not stats['gene'].has_key(result[0]):
        stats['gene'][result[0]] = 0
    stats['gene'][result[0]] += 1
    
    if result[3] == 'yes':
        stats['transcript']['coding'][result[1]] = True
    elif result[3] == 'no':
        stats['transcript']['noncoding'][result[1]] = True
    else:
        stats['transcript']['unknown'][result[1]] = True
                        
    if has_tail and has_bridge and has_link:
        stats['tail_and_bridge_and_link'] += 1
    elif has_tail and has_bridge:
        stats['tail_and_bridge'] += 1
    elif has_bridge and has_link:
        stats['bridge_and_link'] += 1
    elif has_tail and has_link:
        stats['tail_and_link'] += 1
    elif has_tail:
        stats['just_tail'] += 1
    elif has_bridge:
        stats['just_bridge'] += 1
    elif has_link:
        stats['just_link'] += 1

def output_stats(stats, out_file):
    """Outputs stats data into output file"""
    out = open(out_file, 'w')
    out.write('total cleavage sites: %d\n' % stats['cleavage_site'])
    out.write('cleavage sites with tail, bridge, link support: %d\n' % stats['tail_and_bridge_and_link'])
    out.write('cleavage sites with tail, bridge support: %d\n' % stats['tail_and_bridge'])
    out.write('cleavage sites with tail, link support: %d\n' % stats['tail_and_link'])
    out.write('cleavage sites with bridge, link support: %d\n' % stats['bridge_and_link'])
    out.write('cleavage sites with only tail support: %d\n' % stats['just_tail'])
    out.write('cleavage sites with only bridge support: %d\n' % stats['just_bridge'])
    out.write('cleavage sites with only link support: %d\n' % stats['just_link'])
    out.write('total genes: %d\n' % len(stats['gene'].keys()))
    out.write('average cleavage sites per gene: %.1f\n' % (float(stats['cleavage_site'])/len(stats['gene'].keys())))
    out.write('total transcripts: %d\n' % (len(stats['transcript']['coding'].keys()) + len(stats['transcript']['noncoding'].keys())))
    out.write('total coding transcripts: %d\n' % len(stats['transcript']['coding'].keys()))
    out.write('total noncoding transcripts: %d\n' % len(stats['transcript']['noncoding'].keys()))
    out.close()

def prepare_track_header(name, desc, rgb):
    """Creates header for track"""
    return 'track type=bedGraph name="%s" description="%s" visibility=full color=%s' % (name, desc, rgb)

def show_expression(result):
    """Creates bed-graph line depicting expression of cleavage site"""
    return '%s %s %s %s' % (result[5], int(result[6]) - 1, result[6], result[15])

# chrom_proper
chrom_proper = ucsc_chroms(args.ref_genome)

def qpos_to_tpos(a, qpos):
    blocks = align.blocks
    if align.is_reverse:
        strand = '-'
    else:
        strand = '+'
    for i in xrange(len(a['qblocks'])):
        qb0 = a['qblocks'][i][0]
        qb1 = a['qblocks'][i][1]
        if ((qpos >= qb0) and (qpos <= qb1)) or ((qpos <= qb0) and (qpos >= qb1)):
            block = i
            break
    tpos = None
    try:
        if (a['strand'] == '+'):
            tpos = blocks[block][0]+1 + qpos - a['qblocks'][block][0] # added the +1
        else:
            tpos = blocks[block][1] - (qpos - a['qblocks'][block][1])
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
        if (feat.feature == '3UTR') or (feat.feature == 'stop_codon'):
            if (i < half):
                strand = '-'
            else:
                strand = '+'
            return strand
        elif (feat.feature == '5UTR') or (feat.feature == 'start_codon'):
            if (i < half):
                strand = '+'
            else:
                strand = '-'
            return strand
        else:
            return None
        i += 1

def find_link_pairs(a, utr3, homo_len=20, max_mismatch=0):
    global link_pairs_total
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
    if utr3['clipped_pos'] == 'start':
        anchor_read_strand = '-'
    elif utr3['clipped_pos'] == 'end':
        anchor_read_strand = '+'
    
    if anchor_read_strand is None:
        return []
    
    # check if genomic region has polyA - if so, no good
    genome_buffer = 200
    if utr3['feature'].strand == '-':
        span = (int(utr3['cleavage_site']) - genome_buffer, int(utr3['cleavage_site']))
    else:
        span = (int(utr3['cleavage_site']), int(utr3['cleavage_site']) + genome_buffer)
    #genome_seq = self.refseq.GetSequence(align.target, span[0], span[1])
    genome_seq = refseq.fetch(a['target'], span[0]-1, span[1])
    #print 'genome_seq: {}'.format(genome_seq)
    if re.search('A{%s,}' % (homo_len), genome_seq, re.IGNORECASE) or re.search('T{%s,}' % (homo_len), genome_seq, re.IGNORECASE):
        sys.stdout.write('genome sequence has polyAT tract - no reliable link pairs can be retrieved %s %s %s:%s-%s\n' % 
                         (a['contig_seq'], utr3['cleavage_site'], a['target'], span[0], span[1]))
        return []

    mate_loc = {}
    #for read in self.bam.bam.fetch(align.query):
    cigarskip, anchor_rs_plus, anchor_rs_minus, read_count = 0,0,0,0
    for read in r2c.fetch(a['align'].qname):
        read_count += 1
        # skip when both mates mapped to same contig
        if not read.mate_is_unmapped and read.tid == read.rnext:
            continue
        
        # anchor read must be mapped entirely within contig
        if len(read.cigar) > 1:
            cigarskip += 1
            continue
        
        if anchor_read_strand == '+' and (read.is_reverse or read.pos + 1 > utr3['last_matched']):
            anchor_rs_plus += 1
            continue
        if anchor_read_strand == '-' and (not read.is_reverse or read.pos + read.rlen < utr3['last_matched']):
            anchor_rs_minus += 1
            continue
        
        if read.rnext >= 0:
            mate_contig = r2c.getrname(read.rnext)       
            if not mate_loc.has_key(mate_contig):
                link_pairs_total += 1
                mate_loc[mate_contig] = {}
            mate_loc[mate_contig][read.qname] = read
        else:
            pass
#            print 'cannot find unmapped mate %s' % read.qname
            
#    print 'Total reads {}'.format(read_count)
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
        #if good:
            #print '{}-{}'.format(clipped_seq, good)
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

def annotate_cleavage_site(a, feature_list, cleavage_site, clipped_pos, base, min_txt_match_percent=0.6):
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
    
    chrom = proper_chrom(a['target'], chrom_proper=chrom_proper)
    
    # determine which transcript strand can the cleavage site come from
    if clipped_pos == 'start':
        if a['strand'] == '+':
            txt_strand = '-'
        else:
            txt_strand = '+'
    else:
        if a['strand'] == '+':
            txt_strand = '+'
        else:
            txt_strand = '-'
            
    if (txt_strand == '+' and base == 'A') or\
       (txt_strand == '-' and base == 'T'):     
        ests = []
        

        txts_screened = feature_list
        exact_txts = None
        if txts_screened:
            exact_txts = [x for x in txts_screened if x['identical'] == True]
        if exact_txts:
            txts_screened = exact_txts
        
        txts_screened.sort(key=lambda t: abs(t['distance_from_end']))
        closest_txt = txts_screened[0]

        result = {
                'ests': ests,
                'txt': closest_txt,
                'novel': not txts_screened[0]['identical'],
                'within_utr': closest_txt['within_utr'],
                'coord': '%s:%d' % (a['target'], cleavage_site),
                'cleavage_site': cleavage_site,
                'from_end': abs(closest_txt['distance_from_end']),
                'which_end': clipped_pos,
                'base': base
            }

    return result

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
        
    #return result

def find_polyA_cleavage(a, feature_list):
    print
    print 'find_polyA_cleavage'
    print '~'*50
    global filters
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

    min_len = 1
    if filters and filters.has_key('min_at'):
        min_len = filters['min_at']
    mismatch = [1, 1]
    if filters and filters.has_key('max_diff'):
        mismatch = filters['max_diff']
    tail = find_bridge_reads(a, min_len, mismatch, tail=find_tail_contig(a, min_len, mismatch))
    results = []
    for clipped_pos in tail.keys():
        for event in tail[clipped_pos]:
            # contig coordinate of cleavage site                
            last_matched, cleavage_site, base, tail_seq, num_tail_reads, bridge_reads, bridge_clipped_seq = event
            result = annotate_cleavage_site(a, feature_list, cleavage_site, clipped_pos, base)
            #print '*'*25
            #print 'result: {}'.format(result)
            #print '*'*25
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

def find_extended_bridge_reads(a, reads_to_screen, min_len, mismatch, genome_buffer=1000):
    global bridge_total
    """Finds bridge reads where only ending portion represent polyA tail"""
    query_seqs = {}
    for reads in reads_to_screen.values():
        for read in reads:
            query_seqs[read.qname] = read.seq
    
    #entirely_mapped = self.align_transcript_seq(align, query_seqs, 'extended', self.get_full_blat_aln)
    
    clipped_reads = {'start':{}, 'end':{}}
    for clipped_pos, reads in reads_to_screen.iteritems():
        if reads:
            if (clipped_pos == 'start' and a['strand'] == '+') or\
               (clipped_pos == 'end' and a['strand'] == '-'):
                target_coord = [max(0, int(a['align'].reference_start)+1 - genome_buffer), int(a['align'].reference_end)]
            else:
                target_coord = [int(a['align'].reference_start)+1, int(a['align'].reference_end) + genome_buffer]
            
            query_seqs = dict((read.qname, read.seq) for read in reads) #if not entirely_mapped.has_key(read.qname))
            partial_aligns = align_genome_seq(a, query_seqs, target_coord, 'extended-bridge-genome', get_partial_blat_aln)
            
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
                if a['strand'] == '-':
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
                            bridge_total += 1
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
    genome_start, genome_end = max(0, int(align.reference_start) - genome_buffer), int(align.reference_end) + genome_buffer
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
        FNULL = open(os.devnull, 'w')
        task = subprocess.Popen(['blat', target_file, reads_file, aln_file], stdout=FNULL)
        task.communicate()
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
            FNULL = open(os.devnull, 'w')
            task = subprocess.Popen(['blat', target_file, reads_file, aln_file], stdout=FNULL)
            task.communicate()
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
                        
#                    print 'failed %s check %s %s:%s-%s %s %s %s %s %s' % (reason, align.qname, target, align.reference_start, align.reference_end, pos, read[2], read[0].qname, read[0].seq, read[1])
                    
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

def find_bridge_reads(a, min_len, mismatch, genome_buffer=1000, tail=None):
    print
    print 'find_bridge_reads'
    print '~'*50
    global bridge_total
    global filters
    #print "Running find_bridge_reads"
    """Finds bridge reads
    
    It first checks to see clipped reads if the entire clipped portion is A's or T's.
    Then it will check, through find_exteneded_bridge_reads() if the ending portion of 
    the clipped sequence is A's or T's.
    The 2 results will be merged together.
    """
    # used for check if read is mapped to the aligned portion of the contig
    query_bounds = sorted([int(a['qstart']), int(a['qend'])])
    
    # identify clipped reads that are potential pA/pT
    clipped_reads = {'start':{}, 'end':{}}
    second_round = {'start':[], 'end':[]}
    #for read in self.bam.bam.fetch(align.query):
    for read in r2c.fetch(a['align'].query_name):
        if not read.cigar or len(read.cigar) != 2:
            continue
        if (read.cigar[0][0] == 4 or read.cigar[0][0] == 5) or\
           (read.cigar[-1][0] == 4 or read.cigar[-1][0] == 5):
            # clipped at start
            if read.cigar[0][0] == 4 or read.cigar[0][0] == 5:
                clipped_pos = 'start'
                #print read
                #test = raw_input('Press any key to continue')
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
            if filters is not None and filters.has_key('min_bridge_size') and (len(clipped_seq) < filters['min_bridge_size']):
                continue
            
            # reverse complement to be in agreement with reference instead of contig
            clipped_seq_genome = clipped_seq
            if a['strand'] == '-':
                clipped_seq_genome = revComp(clipped_seq)
            
            # check for possible tail (stretch of A's or T's)
            #pos_genome = align.qpos_to_tpos(last_matched)
            pos_genome = qpos_to_tpos(a, last_matched)
                
            picked = False
            for base in ('A', 'T'):
                if is_bridge_read_good(clipped_seq_genome, base, min_len, mismatch):
                    if not clipped_reads[clipped_pos].has_key(last_matched):
                        bridge_total += 1
                        clipped_reads[clipped_pos][last_matched] = {}
                    if not clipped_reads[clipped_pos][last_matched].has_key(base):
                        clipped_reads[clipped_pos][last_matched][base] = []
                    clipped_reads[clipped_pos][last_matched][base].append([read, clipped_seq_genome, pos_genome])
                    picked = True
                    
            if not picked:
                second_round[clipped_pos].append(read)
                    
    extended_clipped_reads = find_extended_bridge_reads(a, second_round, min_len, mismatch)
    merge_clipped_reads(clipped_reads, extended_clipped_reads)

    # Make sure only one end is used for cleavage sites (pick side with more reads)
    #if bool(clipped_reads['start']) and bool(clipped_reads['end']):
    #    if (len(clipped_reads['end']) > len(clipped_reads['start'])):
    #        clipped_reads['start'] = {}

    # filter events against reference sequence
    #filter_vs_reference(align, target, clipped_reads)
    
    print
    print 'clipped_reads:'
    print '@'*50
    print clipped_reads
    

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
                            
    for clipped_pos in clipped_reads:
        if clipped_pos not in results:
            results[clipped_pos] = []
        for pos in clipped_reads[clipped_pos]:                       
            for base in clipped_reads[clipped_pos][pos]:
                pos_genome = clipped_reads[clipped_pos][pos][base][0][2]
                #print 'pos_genome: {}'.format(pos_genome)
                
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
    print 'Done find_bridge_reads: {}'.format(results)
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

def find_tail_contig(a, min_len, mismatch):
    print
    print 'find_tail_contig'
    print '~'*50
    global basic_total
    """Finds contigs that have polyA tail reconstructed"""
    # size of stretch of contig sequence to be added to polyA tail
    # against trasncript sequences to see it polyA tail is genomic
    junction_buffer=50
    results = {}
    
    clipped = {'start':False, 'end':False}
    if int(a['qstart']) > 1:
        clipped['start'] = True
    
    if int(a['qend']) < int(a['align'].infer_query_length(True)):
        clipped['end'] = True
        
    print 'clipped: {}'.format(clipped)
    for clipped_pos in ('start', 'end'):
        if clipped[clipped_pos]:
            if clipped_pos == 'start':
                last_matched = int(a['qstart'])
                clipped_seq = a['contig_seq'][:int(a['qstart'])-1] # used to have -1
                junction_seq = a['contig_seq'][:int(a['qstart']) -1 + junction_buffer] # used to have -1
            else:
                last_matched = int(a['qend'])
                clipped_seq = a['contig_seq'][int(a['qend']):]
                junction_seq = a['contig_seq'][int(a['qend']) - junction_buffer:]
                                    
            print 'clipped_seq: {}'.format(clipped_seq)
            print 'junction_seq: {}'.format(junction_seq)
            cleavage_site = qpos_to_tpos(a, last_matched)
            print 'cleavage_site: {}'.format(cleavage_site)
            clipped_seq_genome = clipped_seq
            if a['strand'] == '-':
                clipped_seq_genome = revComp(clipped_seq)
            
            matched_transcript = in_homopolymer = False
            for base in ('A', 'T'):
                if matched_transcript or in_homopolymer:
                    continue
                                        
                perfect = is_polyA_tail(clipped_seq_genome, base, min_len=1, max_nonAT_allowed=[0, 1])
                #imperfect = self.is_polyA_tail(clipped_seq_genome, base, min_len=4, max_nonAT_allowed=[1, 4])
                imperfect = is_polyA_tail(clipped_seq_genome, base, min_len=min_len, max_nonAT_allowed=mismatch)

                print 'perfect: {}'.format(perfect)
                print 'imperfect: {}'.format(imperfect)
                
                if perfect or imperfect:
                    basic_total += 1
                    # don't need to do the following 2 checks if it's not a potential tail
                    if len(clipped_seq) == 1 and in_homopolymer_neighbor(a['target'], cleavage_site, clipped_seq_genome, clipped_seq[0]):
#                        print '%s : clipped seq in middle of homopolyer run (%s) %s' % (align.qname, clipped_pos, clipped_seq)
                        in_homopolymer = True
                        continue
                    
                    #entirely_mapped = align_transcript_seq(align, target, {align.qname:junction_seq}, 'junction', get_full_blat_aln)
                    #if entirely_mapped.has_key(align.qname):
                    #    print '%s : clipped seq is just transcript seq (%s) %s %d' % (align.qname, clipped_pos, clipped_seq, len(clipped_seq))
                    #    matched_transcript = True
                    #    continue
                    
                    # find reads corresponding to tail
                    num_tail_reads = get_num_tail_reads(a['align'], last_matched)
                                                                    
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
        
    print 'Done find_tail_contig: {}'.format(results)
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
        FNULL = open(os.devnull, 'w')
        task = subprocess.Popen(['blat', target_file, query_file, aln_file], stdout=FNULL)
        task.communicate()
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
        #print "transcript_seq: {}".format(transcripts[tid])
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

def align_genome_seq(a, query_seqs, coord, label, parse_fn):
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
    target_seq = refseq.fetch(a['target'], coord[0], coord[1])
    
    #for cleanup
    tmp_files = []
    
    if args.use_tmp:
        path = '/tmp'
    else:
        path = os.path.dirname(os.path.abspath(args.out))
    
    # create genome reference file for Blatting
    target_file = '%s/%s-genome-%s.fa' % (path, a['align'].qname, label)
    if os.path.exists(target_file):
        os.remove(target_file)
    out = open(target_file, 'w')
    out.write('>%s:%d-%d\n%s\n' % (a['target'], coord[0], coord[1], target_seq))
    out.close()
    tmp_files.append(target_file)

    # create query for Blatting
    query_file = '%s/%s-query-%s.fa' % (path, a['align'].qname, label)
    if os.path.exists(query_file):
        os.remove(query_file)
    out = open(query_file, 'w')
    for query, seq in query_seqs.iteritems():
        out.write('>%s\n%s\n' % (query, seq))
    out.close()
    tmp_files.append(query_file)

    # align query against target
    aln_file = '%s/%s-%s.psl' % (path, a['align'].qname, label)
    if os.path.exists(aln_file):
        os.remove(aln_file)
    tmp_files.append(aln_file)
    try:
        FNULL = open(os.devnull, 'w')
        task = subprocess.Popen(['blat', target_file, query_file, aln_file], stdout=FNULL)
        task.communicate()
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

def output_result(align, result, chrom, link_pairs=[]):
    """Outputs main results in tab-delimited format
    
    The output fields are specified in class variable 'output_fields'
    """
    data = {}
    data['contig'] = align.qname
    
    data['transcript'] = result['txt']['feature'].transcript_id
    data['transcript_strand'] = result['txt']['feature'].strand
    data['chromosome'] = chrom
    if result['txt']['feature'].gene_id:
        data['gene'] = result['txt']['feature'].gene_id
        
    coding_type = result['txt']['coding_type']
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
        
    data['chrom'] = chrom
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
            
    elif not args.no_link:
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
    for field in output_fields:
        if data.has_key(field):
            cols.append(str(data[field]))
        else:
            cols.append('-')
            
    result = '%s\n' % '\t'.join(cols)
    return result

def output(report_lines, bridge_lines=None, link_lines=None):
    """Writes output lines to files"""      
    if out_result is not None and report_lines:
        for line in report_lines:
            out_result.write(line)
    out_result.close()
            
    if bridge_lines and out_bridge_reads is not None:
        for line in bridge_lines:
            out_bridge_reads.write(line)
    if out_bridge_reads is not None:
        out_bridge_reads.close()
        
    if link_lines and out_link_pairs is not None:
        for line in link_lines:
            out_link_pairs.write(line)
    if out_link_pairs is not None:
        out_link_pairs.close()

def fetchUtr3(alignd, features):
    utr3 = None
    for f in features:
        f['identical'] = False
        f['within_utr'] = False
        if (alignd['strand'] == '+'):
            f['clipped_pos'] = 'start'
            # Feature is exactly the 3UTR
            if (f['feature'].end == alignd['align'].reference_end) or (f['feature'].feature == '3UTR'):
                f['identical'] = f['within_utr'] = True
                f['distance_from_end'] = 0
                utr3 = f
            # Feature overlaps the 3UTR
            elif ((f['feature'].start+1) <= alignd['align'].reference_end < f['feature'].end):
                f['within_utr'] = True
                utr3 = f
            f['distance_from_end'] = f['feature'].end - alignd['align'].reference_end
        elif (alignd['strand'] == '-'):
            f['clipped_pos'] = 'end'
            # Feature is exactly the 3UTR
            if ((f['feature'].start+1) == alignd['align'].reference_start) or (f['feature'].feature == '3UTR'):
                f['identical'] = f['within_utr'] = True
                f['distance_from_end'] = 0
                utr3 = f
            # Feature overlaps the 3UTR
            elif ((f['feature'].start+1) < alignd['align'].reference_start <= f['feature'].end):
                f['within_utr'] = True
                utr3 = f
            f['distance_from_end'] = f['feature'].start+1 - alignd['align'].reference_start
    return utr3

lines_result = lines_bridge = lines_link = ''
i=0
for align in aligns:
    if args.c:
        if align.query_name != args.c:
            continue
    if (i > 10):
        continue
    i+=1
    a = {'align': align}
    a['contig_seq'] = contigs.fetch(align.query_name)
    if (align.is_reverse in [True, False]):
        a['strand'] = '-' if align.is_reverse else '+'
    else:
        a['strand'] = None
    print '*'*50
    print 'Looking at contig: {}-{}-{}'.format(align.qname,align.reference_start, align.reference_end)
    print '*'*50
    if (align.reference_start == None) or (align.reference_end == None):
        continue
    # Get target/chromosome
    a['target'] = aligns.getrname(align.tid)
    # Find all overlapping features in gtf annotation
    try:
        feats = features.fetch(a['target'], align.reference_start, align.reference_end)
    except ValueError:
        continue
    # Store relevant features in a list
    feature_list = []
    coding_type = False
    # Try to get strandedness of alignment
    for feat in feats:
        r = {'feature': feat}
        if (feat.strand != a['strand']):
            continue
        if (feat.feature != 'intron'):
            r['coding_type'] = True
        if (feat.feature == 'start_codon'):
            a['cstart'] = feat.start
        elif (feat.feature == 'end_codon'):
            a['cend'] = feat.end
        feature_list.append(r)
    #for f in feature_list:
        #print '\tFeature: gene:{}\ttid:{}\t{}\t{}\t{}\t{}'.format(f['feature'].gene_id,f['feature'].transcript_id,f['feature'].feature,f['feature'].start,f['feature'].end, f['feature'].strand)
    # Get the cleavage sites, clipped positions, and last matched for each feature
    a['inf_strand'] = inferStrand(align, [x['feature'] for x in feature_list])
    if (not a['strand']):
        a['strand'] = a['inf_strand']
    if (a['strand'] is None) and (a['inf_strand'] is None):
        continue
    # Get query blocks
    a['qblocks'] = cigarToBlocks(align.cigar, align.reference_start, a['strand'])[1]
    if a['qblocks'] == None:
        continue
    a['qstart'] = min(a['qblocks'][0][0], a['qblocks'][0][1], a['qblocks'][-1][0], a['qblocks'][-1][1])
    a['qend'] = max(a['qblocks'][0][0], a['qblocks'][0][1], a['qblocks'][-1][0], a['qblocks'][-1][1])
    utr3 = fetchUtr3(a, feature_list)
    result_link = link_pairs = None
    # Find feature closest to 3UTR
    if utr3:
        if (utr3['feature'].strand == '+'):
            utr3['cleavage_site'] = a['align'].reference_end
            utr3['clipped_pos'] = 'end'
            utr3['last_matched'] = a['qend']
        else:
            utr3['cleavage_site'] = a['align'].reference_start+1
            utr3['clipped_pos'] = 'start'
            utr3['last_matched'] = a['qstart']
        #print '\tutr3: gene:{}\ttid:{}\t{}\t{}\t{}\t{}'.format(utr3['feature'].gene_id,utr3['feature'].transcript_id,utr3['feature'].feature,utr3['feature'].start,utr3['feature'].end, utr3['feature'].strand)
        #print '\t{}'.format(utr3)
        utr3['chrom'] = a['target']
        result_link = {'clipped_pos': utr3['clipped_pos'], 'txt': utr3, 'cleavage_site': utr3['cleavage_site'], 'from_end': utr3['distance_from_end'], 'tail_seq': None, 'within_utr': utr3['within_utr'], 'novel': True, 'coord': '{}:{}'.format(a['target'], utr3['cleavage_site']), 'last_matched': utr3['last_matched'], 'ests': []}
    if result_link:
        max_mismatch = 0
        if filters is not None and filters.has_key('max_diff_link'):
            max_mismatch = filters['max_diff_link']
        link_pairs = find_link_pairs(a, utr3)
        if link_pairs and args.output_reads:
            lines_link += output_link_pairs(align, link_pairs)
    if (len(feature_list) == 0):
        continue
    print 'contig\ttarget\ttstrand\ttstart\ttend\tqstrand\tqstart\tqend\tblocks\tqblocks\tcigar\tseq'
    print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(align.qname, \
                                                                a['target'], \
                                                                None, \
                                                                align.reference_start, \
                                                                align.reference_end, \
                                                                a['strand'], \
                                                                a['qstart'], \
                                                                a['qend'], \
                                                                align.blocks, \
                                                                a['qblocks'], \
                                                                a['align'].cigar, \
                                                                a['contig_seq'])
    results = find_polyA_cleavage(a, feature_list)
    if results:
        for result in results:
            lines_result += output_result(a['align'], result, a['target'], link_pairs=link_pairs)
            
            if args.output_reads and result.has_key('bridge_reads') and result['bridge_reads']:
                lines_bridge += output_bridge_reads(align, result)
                
    elif result_link is not None:
        lines_result += output_result(align, result_link, a['target'], link_pairs=link_pairs)
              
# close output streams
if args.combine:
    group_and_filter(lines_result, args.out, filters=filters, make_track=args.track, rgb=args.rgb)
else:
    output(lines_result, lines_bridge, lines_link)

with open('./stuff', 'w') as f:
    f.write('basic: {}\n'.format(basic_total))
    f.write('bridge reads: {}\n'.format(bridge_total))
    f.write('link pairs: {}'.format(link_pairs_total))
