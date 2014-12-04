"""
This program tries to find polyA cleavage sites through short-read assembly.
It is expected that contigs are aligned to contigs, and reads aligned to contigs.
These 2 alignment steps can be performed by trans-ABySS.  The aligners used are
GMAP for contig-genome and BWA-SW for read-contig alignments.
Annotations files for Ensembl, KnownGenes, Refseq, and Aceview are downloaded
from UCSC. EST data(optional) are also downloaded from UCSC.
The analysis can be composed of 2 phases:
1. contig-centric phase - cleavage sites per contig are captured
2. coordinate-centric phase - contigs capturing the same cleavage site
   are consolidated into 1 report where expression/evidence-related data are summed.
   Customized filtering based on evidence data can be performed. 
   
Created by Readman Chiu on 2013-04-17.

v0.02 2013-04-25 by Readman Chiu: 
- line 625 of find_bridge_reads(), should use 'clipped_seq_genome' instead of 'clipped_seq'
  to check for possibility as polyA tail. Because the 'base' that is extracted from
  'clipped_seq' is used in 'annote_cleavage_site()' to check whether the
  'base' makes sense with the orientation of the transcript
- find_polyA_cleavage() will return candidates from both ends now. Before it checks if
  candidates are found from both ends, and if so, discard the contig because a contig 
  cannot have 'real' candidates from both ends.  But scenario occurs where candidates are
  captured from both ends yet only 1 end has valid support and the other doesn't - just a 
  1-base bridge read, which is allowed because during the 'discovery' phase the criteria is 
  intended to be loose.  So hopefully all the false positives will be weeded out during the 
  filtering phase when the maximum bridge reads is required to have at least 2 As
- output_bridge_reads() is changed so that it will no longer capitalize the clipped bases while
  the rest of the read bases are in lowercase.  Because of the complexity involving trimming,
  possibility of clipped bases overlapping tail sequence, extended bridge reads, it takes 
  quite an undertaking to make sure the representation is correct.  In light of time 
  constraint, I abandon this, albeit useful, feature for now.
- I separated the track into a plus-strand track and a minus-strand track so that it's 
  better for visualizing nearby transcripts from both strands
"""
__version__ = '0.02'

import re
import sys
import os
import glob
import subprocess
from utilities import track, assembly
from optparse import OptionParser, OptionGroup
from utilities.overlap_coord import OverlapCoord
from analysis.annotations import est, repeat, ensembl, knownGene, refGene, aceview, ensg
#from analysis.transcript import Transcript
from utilities.tools import get_refseq_from_2bit, reverse_complement, proper_chrom, ucsc_chroms, compare_chr
from utilities.intspan import subsume, overlap, cardinality
from utilities.align_parsers import psl, sam
from utilities.bam import BAM
from time import time
#<add>#
from alignment import Alignment
import pysam
import gmap
import transcript
from itertools import groupby
import pickle
#</add>#

PACKAGE_DIR = '/'.join(os.path.abspath(__file__).split('/')[:-2])

class PolyAFinder:
    """Finds polyA cleavage sites by examining clipped read alignments and read pairs"""

    output_fields = ['gene',
                     'transcript',
                     'transcript_strand',
                     'coding',
                     'contig',
                     'chromosome',
                     'cleavage_site',
                     'within_UTR',
                     'distance_from_annotated_site',
                     'ESTs',
                     'length_of_tail_in_contig',
                     'number_of_tail_reads',
                     'number_of_bridge_reads',
                     'max_bridge_read_tail_length',
                     'bridge_read_identities',
                     'tail+bridge_reads',
                     'number_of_link_pairs',
                     'max_link_pair_length',
                     'link_pair_identities',
                     ]

    def __init__(self, genome, gene_models, out_file, bam_file=None, trim_reads=None, overlap_est=False, no_link=False, use_tmp=False, filters=None):  
        if (genome == 'hg19'):
            genome_file = '/projects/btl/trans-abyss/public_releases/v1.4.8/annotations/hg19/hg19.fa'
        # genome reference
        self.genome = genome
        #<add>#
        self.refseq = pysam.Fastafile(genome_file)
        #print self.refseq.fetch('chr1',45600,45700)
        #self.refseq = get_refseq_from_2bit(genome)
        #print self.refseq.GetSequence('chr1',45600,45700)

        annotation_file = '/projects/clingenetics/reference_data_sets/RefSeq/bgzip_tabix/refseq_hg19.refGene.genePred.utr_cds.sorted.2014-08-21.gtf.gz'
        if not (os.path.isfile(os.path.join(os.path.dirname(bam_file),'transcripts.pkl'))):
            print "transcripts not saved! Saving now..."
            self.transcripts = transcript.Transcript.extract_transcripts(annotation_file)
            with open(os.path.join(os.path.dirname(bam_file),'transcripts.pkl'), 'wb') as f:
                pickle.dump(self.transcripts, f, pickle.HIGHEST_PROTOCOL)
        else:
            print "Loading previously saved transcripts..."
            with open(os.path.join(os.path.dirname(bam_file),'transcripts.pkl'), 'rb') as f:
                self.transcripts = pickle.load(f)
        print "Done with transcripts"

        self.annot = pysam.Tabixfile(annotation_file, parser=pysam.asGTF())

        #</add>#
        #<remove>#
        #self.refseq = get_refseq_from_2bit(genome)
        #</remove>#
        self.chrom_proper = ucsc_chroms(genome)
        
        # gene models
        self.gene_models = gene_models
        #self.txt_overlaps = transcript.Transcript.prepare_overlap(genome, gene_models)
        
        # ESTs
        if overlap_est:
            self.est_overlap = est.prepare_overlap(genome)
        else:
            self.est_overlap = None
            
        # reads-to-contigs BAM
        if bam_file is not None:
            #self.bam = BAM(bam_file)
            self.bam = pysam.Samfile(bam_file, "rb")
                    
        # trim read sequences
        self.trim_reads = trim_reads
        self.set_trimming()
        
        # output
        self.out_file = out_file
        
        # find link pair
        self.no_link = no_link
        
        # use /tmp for Blat alignments
        self.use_tmp = use_tmp
        
        # filters
        self.filters = filters

        # aligner
        self.aligner = 'gmap'

    def prepare_output(self, output_reads=False):
        """Creates file handle(s) for output"""
        if os.path.exists(self.out_file):
            os.remove(self.out_file)
        self.out_result = open(self.out_file, 'a')
        self.out_result.write('%s\n' % '\t'.join(self.output_fields))
        
        if output_reads:
            prefix = os.path.splitext(self.out_file)[0]
            out_bridge_reads_file = prefix + '-reads.fa'
            out_link_pairs_file = prefix + '-link.fa'
            if os.path.exists(out_bridge_reads_file):
                os.remove(out_bridge_reads_file)
            if os.path.exists(out_link_pairs_file):
                os.remove(out_link_pairs_file)
                
            self.out_bridge_reads = open(out_bridge_reads_file, 'a')
            self.out_link_pairs = open(out_link_pairs_file, 'a')            
        else:
            self.out_bridge_reads = self.out_link_pairs = None
            
    def output(self, report_lines, bridge_lines=None, link_lines=None):
        """Writes output lines to files"""      
        if self.out_result is not None and report_lines:
            for line in report_lines:
                self.out_result.write(line)
        self.out_result.close()
                
        if bridge_lines and self.out_bridge_reads is not None:
            for line in bridge_lines:
                self.out_bridge_reads.write(line)
        if self.out_bridge_reads is not None:
            self.out_bridge_reads.close()
            
        if link_lines and self.out_link_pairs is not None:
            for line in link_lines:
                self.out_link_pairs.write(line)
        if self.out_link_pairs is not None:
            self.out_link_pairs.close()
        
            
    def find_polyA_cleavage(self, align):
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
        print 'checking', align.query

        min_len = 1
        if self.filters and self.filters.has_key('min_at'):
            min_len = self.filters['min_at']
        mismatch = [1, 1]
        if self.filters and self.filters.has_key('max_diff'):
            mismatch = self.filters['max_diff']
        print "finding bridge reads..."
        tail = self.find_bridge_reads(align, min_len, mismatch, tail=self.find_tail_contig(align, min_len, mismatch))
        
        results = []
        for clipped_pos in tail.keys():
            for event in tail[clipped_pos]:
                # contig coordinate of cleavage site                
                last_matched, cleavage_site, base, tail_seq, num_tail_reads, bridge_reads, bridge_clipped_seq = event
                
                result = self.annotate_cleavage_site(align, cleavage_site, clipped_pos, base)
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
                        
    def find_link_pairs(self, align, clipped_pos, last_matched, cleavage_site, txt, homo_len=20, max_mismatch=0):
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
        print "clipped_pos: {}".format(clipped_pos)
        anchor_read_strand = None
        if clipped_pos == 'start':
            anchor_read_strand = '-'
        elif clipped_pos == 'end':
            anchor_read_strand = '+'
        
        if anchor_read_strand is None:
            return []
        
        # check if genomic region has polyA - if so, no good
        genome_buffer = 200
        if txt.strand == '-':
            span = (int(cleavage_site) - genome_buffer, int(cleavage_site))
        else:
            span = (int(cleavage_site), int(cleavage_site) + genome_buffer)
        #<remove>#
        #genome_seq = self.refseq.GetSequence(align.target, span[0], span[1])
        #</remove>#
        #<add>#
        genome_seq = self.refseq.fetch(align.target, span[0]-1, span[1])
        #<add>#
        if re.search('A{%s,}' % (homo_len), genome_seq, re.IGNORECASE) or re.search('T{%s,}' % (homo_len), genome_seq, re.IGNORECASE):
            sys.stdout.write('genome sequence has polyAT tract - no reliable link pairs can be retrieved %s %s %s:%s-%s\n' % 
                             (align.query, cleavage_site, align.target, span[0], span[1]))
            return []

        mate_loc = {}
        for read in self.bam.bam.fetch(align.query):
            # skip when both mates mapped to same contig
            if not read.mate_is_unmapped and read.tid == read.rnext:
                continue
            
            # anchor read must be mapped entirely within contig
            if len(read.cigar) > 1:
                continue
            
            if anchor_read_strand == '+' and (read.is_reverse or read.pos + 1 > last_matched):
                continue
            if anchor_read_strand == '-' and (not read.is_reverse or read.pos + read.rlen < last_matched):
                continue
            
            if read.rnext >= 0:
                mate_contig = self.bam.bam.getrname(read.rnext)             
                if not mate_loc.has_key(mate_contig):
                    mate_loc[mate_contig] = {}
                mate_loc[mate_contig][read.qname] = read
            else:
                print 'cannot find unmapped mate %s' % read.qname
        #print "{} reads total".format(k)
        #print "{} reads excluded because mate is not unmapped and tid == rnext".format(i)
        #print "{} reads excluded because cigar length > 1".format(j)
        #print "{} reads excluded because anchor_read_strand == + and either read.is_reverse or read.pos+1 > last_matched".format(y)
        #print "{} reads excluded because anchor_read_strand == - and either read is NOT reverse or read.pos+read.rlen < last_matched".format(x)
        #print "\t* Out of {} reads, {} were because NOT read.is_reverse".format(x,a)
        #print "\t* Out of {} reads, {} were because read.pos + read.rlen < last_matched".format(x,b)
        #sys.exit()
                
        link_pairs = []
        for contig in mate_loc.keys():
            for read in self.bam.bam.fetch(contig):
                if not mate_loc[contig].has_key(read.qname):
                    continue
                
                trimmed_seq = read.seq
                if self.trim_reads:
                    trimmed_seq = self.trim_bases(read.seq, read.qual)

                if trimmed_seq:
                    for base in ('A', 'T'):
                        if self.is_bridge_read_good(trimmed_seq, base, len(trimmed_seq) - max_mismatch, mismatch=[max_mismatch, len(trimmed_seq)]):
                            link_pairs.append([read, mate_loc[contig][read.qname], trimmed_seq])
                            break
        
        return link_pairs
                
    def set_trimming(self):
        """Sets poor base quality string for trimming"""
        # poor quality string used for trimming
        self.poor_quals = ''
        
        if self.trim_reads:
            for i in range(self.trim_reads + 1):
                poor_qual = BAM.int_to_base_qual(i, 33)
                if poor_qual is not None:
                    self.poor_quals += poor_qual
                    
    def trim_bases(self, seq, qual, end=None):
        """Trim poor quality bases from read sequence"""
        if self.poor_quals is None or self.poor_quals == '':
            return seq
        
        match_end = match_start = None
        match_end = re.search(r'[%s]+$' % self.poor_quals, qual)            
        match_start = re.search(r'^[%s]+' % self.poor_quals, qual)
        
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
        
        
    def show_trimmed_bases(self, seq, trimmed_bases):
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
        
        
    def annotate_cleavage_site(self, align, cleavage_site, clipped_pos, base, min_txt_match_percent=0.6):
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
        
        chrom = proper_chrom(align.target, chrom_proper=self.chrom_proper)
        
        # determine which transcript strand can the cleavage site come from
        if clipped_pos == 'start':
            if align.query_strand == '+':
                txt_strand = '-'
            else:
                txt_strand = '+'
        else:
            if align.query_strand == '+':
                txt_strand = '+'
            else:
                txt_strand = '-'
                
        if (txt_strand == '+' and base == 'A') or\
           (txt_strand == '-' and base == 'T'):     
            ests = []
            
            txts_screened = self.get_txts_with_min_match(align, cleavage_site, txt_strand)      
            if txts_screened:
                ## discard transcripts if proposed cleavage site is before 3'UTR
                #txts_screened = [t for t in txts_screened if t['within_utr']]
                #if txts_screened:
                # if a transcript whose transcript end is exactly the same as proposed cleavage site,
                # discard other candidates
                exact_txts = [t for t in txts_screened if t['identical']]
                if exact_txts:
                    txts_screened = exact_txts
                
                # pick transcript closest to proposed cleavage site
                txts_screened.sort(key=lambda t: abs(t['distance_from_end']))
                closest_txt = txts_screened[0]['txt']
                
                # find ESTs that end in exactly same site
                if self.est_overlap:
                    if closest_txt.strand == '+':
                        ests = self.extract_est(chrom, int(closest_txt.txStart), cleavage_site)
                    else:
                        ests = self.extract_est(chrom, cleavage_site, int(closest_txt.txEnd))
                else:
                    ests = None
                    
                result = {
                    'ests': ests,
                    'txt': closest_txt,
                    'novel': not txts_screened[0]['identical'],
                    'within_utr': txts_screened[0]['within_utr'],
                    'coord': '%s:%d' % (align.target, cleavage_site),
                    'cleavage_site': cleavage_site,
                    'from_end': abs(txts_screened[0]['distance_from_end']),
                    'which_end': clipped_pos,
                    'base': base
                }
            else:
                print "%s : %s:%s clipped is not 3' to CDS end of all matching transcripts" % (align.query, align.target, cleavage_site)
                
        else:
            print '%s : %s:%s clipped base(%s) does not match transcript strand(%s)' % (align.query, align.target, cleavage_site, base, txt_strand)
            
        return result
    
    def find_tail_contig(self, align, min_len, mismatch):
        """Finds contigs that have polyA tail reconstructed"""
        # size of stretch of contig sequence to be added to polyA tail
        # against trasncript sequences to see it polyA tail is genomic
        junction_buffer=50
        results = {}
        
        clipped = {'start':False, 'end':False}
        if int(align.qstart) > 1:
            clipped['start'] = True
            
        if int(align.qend) < int(align.query_len):
            clipped['end'] = True
            
        for clipped_pos in ('start', 'end'):
            if clipped[clipped_pos]:
                if clipped_pos == 'start':
                    last_matched = int(align.qstart)
                    clipped_seq = align.contig.sequence[:int(align.qstart)-1]
                    junction_seq = align.contig.sequence[:int(align.qstart) - 1 + junction_buffer]
                else:
                    last_matched = int(align.qend)
                    clipped_seq = align.contig.sequence[int(align.qend):]
                    junction_seq = align.contig.sequence[int(align.qend) - junction_buffer:]
                                        
                cleavage_site = align.qpos_to_tpos(last_matched)
                clipped_seq_genome = clipped_seq
                if align.query_strand == '-':
                    clipped_seq_genome = reverse_complement(clipped_seq)
                
                matched_transcript = in_homopolymer = False
                for base in ('A', 'T'):
                    if matched_transcript or in_homopolymer:
                        continue
                                            
                    perfect = self.is_polyA_tail(clipped_seq_genome, base, min_len=1, max_nonAT_allowed=[0, 1])
                    #imperfect = self.is_polyA_tail(clipped_seq_genome, base, min_len=4, max_nonAT_allowed=[1, 4])
                    imperfect = self.is_polyA_tail(clipped_seq_genome, base, min_len=min_len, max_nonAT_allowed=mismatch)
                    
                    if perfect or imperfect:
                        # don't need to do the following 2 checks if it's not a potential tail
                        if len(clipped_seq) == 1 and self.in_homopolymer_neighbor(align.target, cleavage_site, clipped_seq_genome, clipped_seq[0]):
                            print '%s : clipped seq in middle of homopolyer run (%s) %s' % (align.query, clipped_pos, clipped_seq)
                            in_homopolymer = True
                            continue
                        
                        entirely_mapped = self.align_transcript_seq(align, {align.query:junction_seq}, 'junction', self.get_full_blat_aln)
                        if entirely_mapped.has_key(align.query):
                            print '%s : clipped seq is just transcript seq (%s) %s %d' % (align.query, clipped_pos, clipped_seq, len(clipped_seq))
                            matched_transcript = True
                            continue
                        
                        # find reads corresponding to tail
                        num_tail_reads = self.get_num_tail_reads(align, last_matched)
                                                                        
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
    
    
    def get_num_tail_reads(self, align, last_matched):
        """Reports number of reads spanning cleavage site in contig"""
        num = 0
        for read in self.bam.bam.fetch(align.query):
            if not read.cigar or len(read.cigar) != 1:
                continue
            
            if read.pos + 1 <= last_matched and read.pos + read.alen > last_matched:
                num += 1
        
        return num
                
                
    def check_freq(self, seq):
        """Returns frequency of each base in given sequence"""
        freq = {}
        for nt in seq:
            if not freq.has_key(nt.upper()):
                freq[nt.upper()] = 0
            freq[nt.upper()] += 1
        
        return freq
        
    def is_polyA_tail(self, seq, expected_base, min_len, max_nonAT_allowed):
        """Determines if sequence can be possible tail
        
        min_len = minimum number of expected base in seq
        max_nonAT_allowed(N,M) = N base(s) other than expected base allowed
                                 per stretch of M bases
        """
        result = True
    
        if seq is None or seq == '' or expected_base is None or not expected_base in seq:
            return False
            
        # minimum number of As or Ts
        freq = self.check_freq(seq)
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
            
    def is_bridge_read_good(self, clipped_seq, base, min_len, mismatch):
        """Determines if clipped sequence is possible polyA tail
        
        If clipped_seq is composed of base only, then it is automatically 
        considered a potential bridge read regardless of length
        Otherwise will check the frequecy of 'the other bases' using is_polyA_tail()
        to determine whether it's acceptable
        """
        good = False
        if clipped_seq[0].upper() == base and len(re.sub(r'(.)\1+', r'\1', clipped_seq)) == 1:
            good = True
        elif self.is_polyA_tail(clipped_seq, base, min_len=min_len, max_nonAT_allowed=mismatch):
            good = True

        return good
                        
    def find_bridge_reads(self, align, min_len, mismatch, genome_buffer=1000, tail=None):
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
        for read in self.bam.bam.fetch(align.query):
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
                if self.trim_reads:
                    clipped_seq = self.trim_bases(clipped_seq, clipped_qual, clipped_pos)
                    
                if len(clipped_seq) < 1:
                    continue
                if self.filters is not None and self.filters.has_key('min_bridge_size') and len(clipped_seq) < self.filters['min_bridge_size']:
                    continue
                
                # reverse complement to be in agreement with reference instead of contig
                clipped_seq_genome = clipped_seq
                if align.query_strand == '-':
                    clipped_seq_genome = reverse_complement(clipped_seq)
                
                # check for possible tail (stretch of A's or T's)
                pos_genome = align.qpos_to_tpos(last_matched)
                picked = False
                for base in ('A', 'T'):
                    if self.is_bridge_read_good(clipped_seq_genome, base, min_len, mismatch):
                        if not clipped_reads[clipped_pos].has_key(last_matched):
                            clipped_reads[clipped_pos][last_matched] = {}
                        if not clipped_reads[clipped_pos][last_matched].has_key(base):
                            clipped_reads[clipped_pos][last_matched][base] = []
                        clipped_reads[clipped_pos][last_matched][base].append([read, clipped_seq_genome, pos_genome])
                        picked = True
                        
                if not picked:
                    second_round[clipped_pos].append(read)
                        
        extended_clipped_reads = self.find_extended_bridge_reads(align, second_round, min_len, mismatch)
        self.merge_clipped_reads(clipped_reads, extended_clipped_reads)
                    
        # filter events against reference sequence
        self.filter_vs_reference(align, clipped_reads)
        
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
    
    def merge_clipped_reads(self, clipped_reads, extended_clipped_reads):
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
                            
    def find_extended_bridge_reads(self, align, reads_to_screen, min_len, mismatch, genome_buffer=1000):
        """Finds bridge reads where only ending portion represent polyA tail"""
        query_seqs = {}
        for reads in reads_to_screen.values():
            for read in reads:
                query_seqs[read.qname] = read.seq
        
        entirely_mapped = self.align_transcript_seq(align, query_seqs, 'extended', self.get_full_blat_aln)
        
        clipped_reads = {'start':{}, 'end':{}}
        for clipped_pos, reads in reads_to_screen.iteritems():
            if reads:
                if (clipped_pos == 'start' and align.query_strand == '+') or\
                   (clipped_pos == 'end' and align.query_strand == '-'):
                    target_coord = [int(align.tstart) - genome_buffer, int(align.tend)]
                else:
                    target_coord = [int(align.tstart), int(align.tend) + genome_buffer]
                
                query_seqs = dict((read.qname, read.seq) for read in reads if not entirely_mapped.has_key(read.qname))
                partial_aligns = self.align_genome_seq(align, query_seqs, target_coord, 'extended-bridge-genome', self.get_partial_blat_aln)
                
                read_objs = dict((read.qname, read) for read in reads)
                
                for read_name, mapped_coord in partial_aligns.iteritems():
                    if mapped_coord[0] == 0:
                        clipped_seq = read_objs[read_name].seq[mapped_coord[1]:]
                    else:
                        clipped_seq = read_objs[read_name].seq[:mapped_coord[0]]
                        
                    if self.filters is not None and self.filters.has_key('min_bridge_size') and len(clipped_seq) < self.filters['min_bridge_size']:
                        continue
                                            
                    # reverse complement to be in agreement with reference instead of contig
                    clipped_seq_genome = clipped_seq
                    if align.query_strand == '-':
                        clipped_seq_genome = reverse_complement(clipped_seq)
                        
                    if mapped_coord[0] == 0:
                        last_matched = read_objs[read_name].pos + mapped_coord[1]
                        pos_genome = target_coord[0] + mapped_coord[3] - 1
                    else:
                        last_matched = read_objs[read_name].pos - mapped_coord[0]
                        pos_genome = target_coord[0] + mapped_coord[2]
                                            
                    for base in ('A', 'T'):
                        if self.is_bridge_read_good(clipped_seq, base, min_len, mismatch):
                            #print 'possible extended bridge reads', align.query, read_objs[read_name].qname, read_objs[read_name].seq, is_seed, clipped_seq_genome
                            if not clipped_reads[clipped_pos].has_key(last_matched):
                                clipped_reads[clipped_pos][last_matched] = {}
                            if not clipped_reads[clipped_pos][last_matched].has_key(base):
                                clipped_reads[clipped_pos][last_matched][base] = []
                            clipped_reads[clipped_pos][last_matched][base].append([read_objs[read_name], clipped_seq_genome, pos_genome])
                            
        return clipped_reads
                                                                                                    
    def align_genome_seq(self, align, query_seqs, coord, label, parse_fn):
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
        target_seq = self.refseq.fetch(align.target, coord[0]-1, coord[1])
        
        #for cleanup
        tmp_files = []
        
        if self.use_tmp:
            path = '/tmp'
        else:
            path = os.path.dirname(os.path.abspath(self.out_file))
        
        # create genome reference file for Blatting
        target_file = '%s/%s-genome-%s.fa' % (path, align.query, label)
        if os.path.exists(target_file):
            os.remove(target_file)
        out = open(target_file, 'w')
        out.write('>%s:%d-%d\n%s\n' % (align.target, coord[0], coord[1], target_seq))
        out.close()
        tmp_files.append(target_file)

        # create query for Blatting
        query_file = '%s/%s-query-%s.fa' % (path, align.query, label)
        if os.path.exists(query_file):
            os.remove(query_file)
        out = open(query_file, 'w')
        for query, seq in query_seqs.iteritems():
            out.write('>%s\n%s\n' % (query, seq))
        out.close()
        tmp_files.append(query_file)

        # align query against target
        aln_file = '%s/%s-%s.psl' % (path, align.query, label)
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
        
    def align_transcript_seq(self, align, query_seqs, label, parse_fn):
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
        
        if self.use_tmp:
            path = '/tmp'
        else:
            path = os.path.dirname(os.path.abspath(self.out_file))
        
        # create transcript sequence file for Blatting
        target_file = '%s/%s-target-%s.fa' % (path, align.query, label)
        if os.path.exists(target_file):
            os.remove(target_file)
        target_empty = self.extract_transcript_seq(align, target_file)
        #target_empty = self.map_contigs_to_transcripts()
        tmp_files.append(target_file)

        # create query for Blatting
        query_file = '%s/%s-query-%s.fa' % (path, align.query, label)
        if os.path.exists(query_file):
            os.remove(query_file)
        out = open(query_file, 'w')
        for query, seq in query_seqs.iteritems():
            out.write('>%s\n%s\n' % (query, seq))
        out.close()
        tmp_files.append(query_file)

        # align query against target
        aln_file = '%s/%s-%s.psl' % (path, align.query, label)
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
    
    def filter_vs_reference(self, align, clipped_reads):
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
        
        if self.use_tmp:
            path = '/tmp'
        else:
            path = os.path.dirname(os.path.abspath(self.out_file))
            
        # create genome sequence file for Blatting
        target_file = '%s/%s-genome.fa' % (path, align.query)
        if os.path.exists(target_file):
            os.remove(target_file)
        out = open(target_file, 'w')
        genome_start, genome_end = int(align.tstart) - genome_buffer, int(align.tend) + genome_buffer
        #seq = self.refseq.GetSequence(align.target, genome_start, genome_end)
        seq = self.refseq.fetch(align.target, genome_start-1, genome_end)
        out.write('>%s:%s-%s\n%s\n' % (align.target, genome_start, genome_end, seq))
        out.close()
        tmp_files.append(target_file)
        

        # create reads for Blatting
        reads_file = '%s/%s-bridge.fa' % (path, align.query)
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
        aln_file = '%s/%s-bridge-genome.psl' % (path, align.query)
        if os.path.exists(aln_file):
            os.remove(aln_file)
        tmp_files.append(aln_file)
        try:
            subprocess.call(['blat', target_file, reads_file, aln_file])
        except CalledProcessError as err:
            sys.stderr.write('error running blat:%s' % err.cmd)
        else:
            if os.path.exists(aln_file):
                fully_aligned_genome = self.get_full_blat_aln(aln_file)
                
        # create transcript sequence file for Blatting
        target_file = '%s/%s-txt.fa' % (path, align.query)
        if os.path.exists(target_file):
            os.remove(target_file)
        target_empty = self.extract_transcript_seq(align, target_file)
        tmp_files.append(target_file)
        
        # align reads against transcripts
        fully_aligned_txt = None
        aln_file = '%s/%s-bridge-txt.psl' % (path, align.query)
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
                    fully_aligned_txt = self.get_full_blat_aln(aln_file)
                                        
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
                            
                        if read[2] is not None and self.in_homopolymer_neighbor(align.target, read[2], read[1], base):
                            bad_neighbor = True
                            break
                                            
                    if is_genomic or is_transcriptome or bad_neighbor:
                        if is_genomic:
                            reason = 'genome seq'
                        elif is_transcriptome:
                            reason = 'transcriptome seq'
                        else:
                            reason = 'neighborhood'
                            
                        print 'failed %s check %s %s:%s-%s %s %s %s %s %s' % (reason, align.query, align.target, align.tstart, align.tend, pos, read[2], read[0].qname, read[0].seq, read[1])
                        
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
                                                
    def get_full_blat_aln(self, aln_file):
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
    
    def get_partial_blat_aln(self, aln_file):
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
    
    def get_blat_matches(self, aln_file):
        """Reports match length of every BLAT aligment
        
        This is used for calculating the percentage of contig that is mapped to a transcript
        """
        matches = {}
        for line in open(aln_file, 'r'):
            if not re.search('^\d', line):
                continue
            cols = line.rstrip('\n').split('\t')
            query, qsize, qstart, qend, target = cols[9:14]
            
            match_len = int(qend) - int(qstart)
            if not matches.has_key(target) or matches[target] < match_len:
                matches[target] = match_len
            
        return matches

    def extract_aligns(self, alns):
        try:
            chimeric_aligns = {'gmap': gmap.find_chimera,}[self.aligner](alns, self.bam)
            if chimeric_aligns:
                return chimeric_aligns
            else:
                return [{'gmap': gmap.find_single_unique,}[self.aligner](alns, self.bam)]
        except:
            sys.exit("Can't convert \"%s\" alignments - abort" % aligner)

    def extract_transcript_seq(self, align, out_file):
#    def map_contigs_to_transcripts(self):
#        for transcript in self.transcripts:
#        """Extracts transcripts overlapping an alignment and outputs their
#        sequences to the given file
#        
#        Returns empty=True if it fails to find overlapping transcripts
#        """
#
#        for contig, group in groupby(self.bam.fetch(until_eof=True), lambda x: x.qname):
#            print "analyzing {} contig".format(contig)
#            alns = list(group)
#            aligns = self.extract_aligns(alns)
#            if aligns is None:
#                sys.stdout.write('no valid alignment: %s' % contig)
#                continue
#
#            chimera = True if len(aligns) > 1 else False
#            chimera_block_matches = []
#            for align in aligns:
#                if align is None:
#                    sys.stdout.write('bad alignment: %s\n' % contig)
#                    continue
#                if re.search('[._Mm]', align.target):
#                    sys.stdout.write('skip target:%s %s\n' % (contig, align.target))
#                    continue
#
#                within_intron = []
#                within_exon = []
#                transcripts_mapped = Set()
#                events = []
#                for gtf in self.annot.fetch(align.target, align.tstart, align.tend):
#                    if gtf.feature == 'exon':
#                        transcripts_mapped.add(gtf.transcript_id) 
#                    elif (gtf.feature == 'intron') and (not chimera) and ((align.tstart >= gtf.start) and (align.tend <= gtf.end)):
#                        match = self.match_exon((align.tstart, align.tend), (gtf.start, gtf.end))
#                        within_intron.append((gtf, match))
#
#                if transcripts_mapped:
#                    mappings = []
#                    # key = transcript name, value = "matches"
#                    all_block_matches = {}
#                    # Transcript objects that are fully matched
#                    full_matched_transcripts = []
#                    for txt in transcripts_mapped:
#                        block_matches = self.map_exons(align.blocks, transcripts[txt].exons)
#                        all_block_matches[txt] = block_matches
#                        mappings.append((transcripts[txt], block_matches))
#                
#                        if not chimera and self.is_full_match(block_matches):
#                            full_matched_transcripts.append(transcripts[txt])
#                    
#                    # report mapping
#                    best_mapping = Mapping.pick_best(mappings, align, debug=self.debug)
#                    self.mappings.append(best_mapping)
#
#                    if not full_matched_transcripts:
#                        best_transcript = best_mapping.transcripts[0]
#                        events = self.find_events({best_transcript.id:all_block_matches[best_transcript.id]},align,{best_transcript.id:best_transcript})
#                        for event in events:
#                            event.contig_sizes.append(len(contig_seq))
#                        if events:
#                            self.events.extend(events)
#                        elif self.debug:
#                            sys.stdout.write('%s - partial but no events\n' % align.query)
#                    if chimera:
#                        chimera_block_matches.append(all_block_matches)
#                elif not chimera:
#                    if within_exon:
#                        sys.stdout.write("contig mapped within single exon: %s %s:%s-%s %s\n" % (contig,align.target, align.tstart, align.tend, within_exon[0]))
#                    elif within_intron:
#                        sys.stdout.write("contig mapped within single intron: %s %s:%s-%s %s\n" % (contig,align.target, align.tstart, align.tend, within_exon[0]))
#            if chimera and chimera_block_matches:
#                if len(chimera_block_matches) == len(aligns):
#                    fusion = FusionFinder.find_chimera(chimera_block_matches, transcripts, aligns, contig_seq,exon_bound_only=self.fusion_conditions['exon_bound_only'],coding_only=self.fusion_conditions['coding_only'],sense_only=self.fusion_conditions['sense_only'])
#                    if fusion:
#                        homol_seq, homol_coords = None, None
#                        if self.aligner.lower() == 'gmap':
#                            homol_seq, homol_coords = gmap.find_microhomology(alns[0], contig_seq)
#                        if homol_seq is not None:
#                            fusion.homol_seq.append(homol_seq)
#                            fusion.homol_coords.append(homol_coords)
#                        fusion.contig_sizes.append(len(contig_seq))
#                        self.events.append(fusion)
#                
#        # expand contig span
#        for event in self.events:
#            # if contig_support_span is not defined (it can be pre-defined in ITD)
#            # then set it
#            if not event.contig_support_span:
#                event.contig_support_span = event.contig_breaks
#                if event.rearrangement == 'ins' or event.rearrangement == 'dup':
#                    expanded_contig_breaks = expand_contig_breaks(event.chroms[0], 
#                                                                  event.breaks, 
#                                                                  event.contigs[0], 
#                                                                  [event.contig_breaks[0][0] + 1, event.contig_breaks[0][1] - 1], 
#                                                                  event.rearrangement, 
#                                                                  self.ref_fasta,
#                                                                  self.contigs_fasta,
#                                                                  self.debug)
#                    if expanded_contig_breaks is not None:
#                        event.contig_support_span = [(expanded_contig_breaks[0] - 1, expanded_contig_breaks[1] + 1)]
                        

        chrom = proper_chrom(align.target, chrom_proper=self.chrom_proper)
        out = open(out_file, 'w')
        empty = True
        for txt in self.transcripts:
            transcript_seq = self.transcripts[txt].get_sequence(self.refseq, chrom=align.target).upper()
            out.write('>%s\n%s\n' % (self.transcripts[txt].id, transcript_seq))
            empty = False             
        out.close()
        
        return empty
    
    def in_homopolymer_neighbor(self, chrom, pos, tail, base):
        """Checks if tail is juxtaposed against a homopolymer of the given base in genome
        
        The window size for checking is 5 or 2*length of the tail sequence (tail), whichever larger.
        The homopolyer run is constructed by making a string of the given base (base) with length 
        equal to the frequency of that base in the tail sequence.
        If the homopolyer is found immedicately before or after the given position (pos),
        then the result is True.
        """
        min_len = 5
        length = max(min_len, 2 * len(tail))
        #neighbor_seq = self.refseq.GetSequence(chrom, pos-length, pos+length)
        neighbor_seq = self.refseq.fetch(chrom, pos-length-1, pos+length)
        homo = ''
        freq = self.check_freq(tail)
        
        if freq.has_key(base):
            for i in range(length):
                homo += base.upper()
                    
            m = re.search(homo, neighbor_seq, re.IGNORECASE)
            if m:               
                # if homopolymer is touching the middle base of neighbor_seq
                if m.start() <= len(neighbor_seq)/2 + 1 and m.end() - 1 >= len(neighbor_seq)/2 - 1:
                    return True
                        
        return False
    
    def extract_est(self, chrom, start, end):
        """Extracts ESTs given a coordinate
        
        If the EST is on + strand, 'end' is expected to be the cleavage site
        If the EST in on - strand, 'start' is expected to be the cleavage site
        """
        print 'finding ESTs for %s:%s-%s' % (chrom, start, end)
        ests = []
        lines = self.est_overlap.overlap(chrom, min(start, end), max(start, end))
        
        for line in lines:
            e = est.parse_line(line)
            if (e.strand == '-' and int(e.txStart) == int(start)) or (e.strand == '+' and int(e.txEnd) == int(end)):
                ests.append(e)
                
        return ests
      
    def extract_txts(self, chrom, start, end, cleavage_site, txt_strand):
        """Extracts transcripts overlapping given region and strand
        
        This method first finds all transcripts that overlap the given region (chrom, start, end).
        It will skip ignore transcripts that don't match the given strand.
        For each kept transcript, it will keep 4 pieces of info:
        1. the transcript object itself (txt)
        2. whether the given cleavage site is identical to transcript end (identical)
        3. whether the cleavage site is 3' to the beginning of the 3'UTR of the transcript (within_utr)
        4. distance of cleavage site from transcript end (distance_from_end)
        The 4 pieces of info are kept in a dictionary for each matching transcript,
        and the final return variable is a list of dictionary
        """     
        txts = transcript.find_overlaps(self.txt_overlaps, self.gene_models.split(','), chrom, start, end)
                
        results = []
        for txt in txts:
            # only keep transcripts matching the given strand
            if txt.strand == txt_strand:
                if (txt.strand == '+' and int(txt.txEnd) == cleavage_site) or\
                   (txt.strand == '-' and int(txt.txStart) + 1 == cleavage_site):
                    identical = True
                    within_utr = True
                    distance_from_end = 0
                else:
                    identical = False
                    
                    within_utr = False
                    utr3 = txt.utr(3)
                    # not actually within 3UTR, but overlaps
                    if utr3:
                        if (txt.strand == '+' and cleavage_site >= utr3[0]) or\
                           (txt.strand == '-' and cleavage_site <= utr3[1]):
                            within_utr = True                       
                    
                    if txt.strand == '+':
                        annotated_end = int(txt.txEnd)
                    else:
                        annotated_end = int(txt.txStart) + 1
                    distance_from_end = annotated_end - cleavage_site
                    
                results.append({'txt': txt,
                                'identical': identical,
                                'within_utr': within_utr,
                                'distance_from_end': distance_from_end
                                })
        return results
            
                            
    def get_txts_with_min_match(self, align, cleavage_site, txt_strand):
        """Extracts transcripts matching alignment and strand with a minimum match percentage
        
        It extracts transcripts matching the alignment region and given strand using extract_txts().
        It asks align_transcript_seq() to align the contig sequence to any transcript overlapping the 
        alignment region (regardless of strand).
        It then screens all the transcripts returned from extract_txts() to make sure 
        a) at least N% of the contig is aligned to the transcript in question
        b) at least N% of the span of the alignment's genomic target span is accounted for by the 
           overlap between the transcript span and the alignment target span
        N = min_txt_match_percent, and set as 60%
        Condition b) is designed to screen out wrong mapping of a contig to a much smaller transcript
        i.e. you want
                     |-1-|------>-------|-2-| contig)
                |--------|--------------|--|  transcript1
             but not
                     |---|--------------|---| contig
                   |-<-|                      trancript2
        (transcript1 may not be able to map to contig because the clipped base indicates the transcript should b
         on -ve strand, which only transcript2 fulfills, but the clipped base is way over at block 2)
        The return variable has the same structure of that from extract_txts(): a list of dictionaries.
        """
        min_txt_match_percent=0.6
        
        chrom = proper_chrom(align.target, chrom_proper=self.chrom_proper)
        matches = self.align_transcript_seq(align, {align.query:align.contig.sequence}, 'contig', self.get_blat_matches)
        #txts = self.extract_txts(chrom, int(align.tstart), int(align.tend), cleavage_site, txt_strand)
        txts_screened = []
        # check if polyA base is matching transcript strand
        for txt in self.transcripts:
            if matches.has_key(self.transcripts[txt].id):
                match_percent = float(matches[self.transcripts[txt].id])/float(align.query_len)
                span_overlap_percent = float(overlap([int(self.transcripts[txt].cds_start), int(self.transcripts[txt].cds_end)], [align.tstart, align.tend]))/float(cardinality([align.tstart, align.tend]))            
                if match_percent >= min_txt_match_percent and span_overlap_percent >= min_txt_match_percent:
                    txts_screened.append(self.transcripts[txt])
                else:
                    print '%s : %s skipped - match percent too low %d/%d=%2f' % (align.query, self.transcripts[txt].id, matches[self.transcripts[txt].id], align.query_len, match_percent)
                    continue
                
        return txts_screened

    def contains_utr(self, align):
        """Extract transcript that best matches alignment if contig reconstructs some part of 3'UTR
        
        This method is designed for link pairs retrieval.  Link pairs are only attempted if the end of 
        a contig can be mapped to 3'UTR of an annoated trancript.
        It first looks for transcripts that have their 3'UTRs overlapping either the start or end of
        the contig alignment.
        For each end (start or end), it keeps the transcript where the alignment end is closest to 
        the transcript end.  If transcripts can be found at both the start and end, again, the one closest
        to the alignment end is chosen.
        The return dictionary has the same fields as that returned when an actual tail is found, with
        certain fields set to None or empty.
        """
        
        txts_start = [t for t in self.get_txts_with_min_match(align, int(align.tstart), '-') if t['within_utr']]
        txts_end = [t for t in self.get_txts_with_min_match(align, int(align.tend), '+') if t['within_utr']]
                
        txt_start = None
        cleavage_site_guess = None
        if txts_start:
            txts_start.sort(key=lambda t: abs(t['distance_from_end']))
            txt_start = txts_start[0]
            
        txt_end = None
        if txts_end:
            txts_end.sort(key=lambda t: abs(t['distance_from_end']))
            txt_end = txts_end[0]
                
        txt = None
        if txt_start is not None or txt_end is not None:
            if txt_start is not None and (txt_end is None or txt_start['distance_from_end'] < txt_end['distance_from_end']):
                txt = txt_start
                cleavage_site_guess = int(align.tstart)
                if align.query_strand == '+':
                    clipped_pos = 'start'
                    last_matched = int(align.qstart)
                else:
                    clipped_pos = 'end'
                    last_matched = int(align.qend)
                    
            else:
                txt = txt_end
                cleavage_site_guess = int(align.tend)
                if align.query_strand == '+':
                    clipped_pos = 'end'
                    last_matched = int(align.qend)
                else:
                    clipped_pos = 'start'
                    last_matched = int(align.qstart)
            
        if txt is not None:
            return {
                'ests': [],
                'txt': txt['txt'],
                'novel': True,
                'within_utr': True,
                'coord': "%s:%d" % (align.target, cleavage_site_guess),
                'from_end': txt['distance_from_end'],
                'clipped_pos': clipped_pos,
                'last_matched': last_matched,
                'cleavage_site': cleavage_site_guess,
                'tail_seq': None,
            }
                
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

    def output_bridge_reads(self, align, result):
        """Outputs bridges reads in FASTA format"""
        out = ''
        if result.has_key('bridge_reads'):
            for i in range(len(result['bridge_reads'])):
                read = result['bridge_reads'][i]                
                out += ">%s %s %s\n%s\n" % (read.qname, align.query, result['coord'], read.seq.upper())                 
        return out
    
    def output_link_pairs(self, align, reads):
        """Outputs link pairs in FASTA format"""
        out = ''
        for r in reads:
            trimmed_seq = self.show_trimmed_bases(r[0].seq, r[2])
            num_trimmed_bases = r[0].rlen - len(r[2])
            if r[1].is_reverse:
                anchor_direction = 'L'
            else:
                anchor_direction = 'R'
                
            if r[0].is_unmapped:
                mate_contig = 'unmapped'
            else:
                mate_contig = self.bam.bam.getrname(r[0].tid)
                
            out += '>%s %s %s %s %s trimmed:%s\n%s\n' % (r[0].qname, align.query, r[1].pos, anchor_direction, 
                                                         mate_contig, num_trimmed_bases, trimmed_seq)                   
        return out
    
    @classmethod
    def group_and_filter(cls, path, out_file, filters=None, make_track=None, rgb='0,0,0'):
        """Consolidates contig-centric results into coordinate-centric results
        
        path = directory of tab-delimited results files
        """
        groups = {}
        
        results_files = [f for f in glob.glob(os.path.join(path, '*')) if not '-link' in f and not '-reads' in f]
        for results_file in sorted(results_files):
            print results_file
            for line in open(results_file, 'r'):
                cols = line.rstrip('\n').split('\t')
                if cols[0] == cls.output_fields[0]:
                    continue
                
                chrom, cleavage_site = cols[5], cols[6]
                if not groups.has_key(chrom):
                    groups[chrom] = {}
                if not groups[chrom].has_key(cleavage_site):
                    groups[chrom][cleavage_site] = []
                groups[chrom][cleavage_site].append(cols)
                
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
        out.write('%s\n' % '\t'.join(cls.output_fields))
        
        track_plus = []
        track_minus = []
        for chrom in sorted(groups.keys(), cmp=compare_chr):
            for cleavage in sorted(groups[chrom].keys()):
                results = groups[chrom][cleavage]
                # if more than one contig reports same cleavage site, add up the support numbers
                if len(results) > 1:
                    result = cls.merge_results(results)
                else:
                    result = results[0]
                                        
                if filters is not None and filters.has_key('min_bridge_size') and int(result[13]) < filters['min_bridge_size']:
                    continue
                    
                out.write('%s\n' % '\t'.join(result))
                if make_track is not None:
                    if result[2] == '+':
                        track_plus.append(cls.show_expression(result))
                    else:
                        track_minus.append(cls.show_expression(result))
                
                # stats
                cls.update_stats(stats, result)
                    
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
        
    @classmethod
    def output_track(cls, out_file, name, desc, rgb='0,0,0', plus=[], minus=[]):
        """Outputs bedgraph, separate tracks for plus and minus strands"""
        out = open(out_file, 'w')
        # plus strand
        out.write('%s\n' % cls.prepare_track_header(name + ' (+)', desc + ' plus strand', rgb))
        for line in plus:
            out.write('%s\n' % line)
        # minus strand
        out.write('%s\n' % cls.prepare_track_header(name + '(-)', desc + ' minus strand', rgb))
        for line in minus:
            out.write('%s\n' % line)
        out.close()
                            
    @classmethod
    def merge_results(cls, results):        
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
    
    @classmethod
    def update_stats(cls, stats, result):
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

    @classmethod
    def output_stats(cls, stats, out_file):
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

    @classmethod
    def prepare_track_header(cls, name, desc, rgb):
        """Creates header for track"""
        return 'track type=bedGraph name="%s" description="%s" visibility=full color=%s' % (name, desc, rgb)
    
    @classmethod
    def show_expression(cls, result):
        """Creates bed-graph line depicting expression of cleavage site"""
        return '%s %s %s %s' % (result[5], int(result[6]) - 1, result[6], result[15])
            
def main(args, options):
    out_file = options.out_file
    filters = {}
    
    # group and filter
    if options.combine:
        filters['min_bridge_size'] = options.min_bridge_size
        PolyAFinder.group_and_filter(options.combine, out_file, filters=filters, make_track=options.track, rgb=options.rgb)
        
    # find cleavage sites
    else:
        aln_file, contigs_file, genome = options.aln_file, options.ctg_file, options.genome
        filters['min_at'] = options.min_at
        filters['max_diff'] = options.max_diff
        filters['max_diff_link'] = options.max_diff_link
                
        pf = PolyAFinder(genome, options.models, out_file, 
                         bam_file=options.bam_file, 
                         trim_reads=options.trim_reads, 
                         overlap_est=options.overlap_est, 
                         no_link=options.no_link,
                         use_tmp=options.use_tmp,
                         filters=filters)
        
        # parse alignments
        filters = {'unique':True, 'bestn':1, 'match':60.0, 'identity':90.0}
        ext = os.path.splitext(aln_file)[1]
        if ext == '.psl':
            aligns = psl.parse(aln_file, filters=filters, noline=False)
        elif ext == '.sam':
            aligns = sam.parse(aln_file, filters=filters, header=True)
        #<add>#
        elif ext == '.bam':
            aligns = pysam.Samfile(aln_file, "rb")
            alist = []
            for aln in aligns:
                if (not aln.is_unmapped):
                    alist.append(Alignment.from_alignedRead(aln, aligns))
            aligns = alist
        else:
            print "Unknown file format: {}, exiting".format(aln_file)
            sys.exit()
        #</add>#
        
        # parse and store contig sequences
        ass = assembly.Assembly(None, fasta=contigs_file)
        contigs = ass.get_contigs(ids=[a.query for a in aligns], sequence=True)
        contig_dict = dict((c.num, c) for c in contigs)
        for align in aligns:
            print "Looking at alignment: {} {} {}".format(align.target, align.tstart, align.tend)
            if contig_dict.has_key(align.query):
                align.contig = contig_dict[align.query]
        
        # open output streams
        pf.prepare_output(output_reads=options.output_reads)
        lines_result = lines_bridge = lines_link = ''
        
        for align in aligns:
            # find link pair support
            link_pairs = None
            result_link = None
            if not pf.no_link:
                result_link = pf.contains_utr(align)
                if result_link:
                    if pf.bam:
                        max_mismatch = 0
                        if pf.filters is not None and pf.filters.has_key('max_diff_link'):
                            max_mismatch = pf.filters['max_diff_link']
                        link_pairs = pf.find_link_pairs(align, result_link['clipped_pos'], result_link['last_matched'], result_link['coord'].split(':')[1], result_link['txt'], max_mismatch=max_mismatch)
                        if link_pairs and options.output_reads:
                            lines_link += pf.output_link_pairs(align, link_pairs)
                       
            # find tail and bridge support
            results = pf.find_polyA_cleavage(align)
            if results:
                for result in results:
                    lines_result += pf.output_result(align, result, link_pairs=link_pairs)
                    
                    if options.output_reads and result.has_key('bridge_reads') and result['bridge_reads']:
                        lines_bridge += pf.output_bridge_reads(align, result)
                        
            elif result_link is not None:
                lines_result += pf.output_result(align, result_link, link_pairs=link_pairs)
                      
        # close output streams
        pf.output(lines_result, lines_bridge, lines_link)
                                 

if __name__ == '__main__':
    usage = "Usage: %prog <options>"

    parser = OptionParser(usage=usage, version="%prog " + iddddd inputs = OptionGroup(parser, 'inputs')
    inputs.add_option('-a', '--aln', dest='aln_file', help='alignment file')
    inputs.add_option('-c', '--ctg', dest='ctg_file', help='contig file')
    inputs.add_option('-g', '--genome', dest='genome', help='genome')
    inputs.add_option('-m', '--models', dest='models', help='gene models:e,r,a,k', default='e,r,a,k')
    inputs.add_option('-b', '--bam_file', dest='bam_file', help='reads-to-contigs BAM file')
    parser.add_option_group(inputs)
    
    outputs = OptionGroup(parser, 'outputs')
    outputs.add_option('-o', '--out', dest='out_file', help='output file')
    outputs.add_option('-r', '--output_reads', dest='output_reads', help='output read sequences', action='store_true', default=False)
    parser.add_option_group(outputs)
    
    computation = OptionGroup(parser, 'computation')
    computation.add_option('-t', '--trim', dest='trim_reads', help='trim poor bases off quality equal or below. Default: 3 (0 = no trimming)', default=3, type='int')
    computation.add_option('--use_tmp', dest='use_tmp', help='use /tmp space', action='store_true', default=False)
    computation.add_option('--no_link', dest='no_link', help='do not find link pairs', action='store_true', default=False)
    computation.add_option('-e', '--overlap_est', dest='overlap_est', help='overlap ESTs', action='store_true', default=False)
    parser.add_option_group(computation)
    
    read_filters = OptionGroup(parser, "read_filters")
    read_filters.add_option('--min_at', dest='min_at', help='minimum number of A|T base in tail. Default:4', type='int', default=4)
    read_filters.add_option('--max_diff', dest='max_diff', help='maximum rate of non A|T base allowed in tail Default:1 5', 
                       nargs=2, default=[1, 5])
    read_filters.add_option('--max_diff_link', dest='max_diff_link', help='maximum number of non A|T bases in entire link read. Default:2', 
                       default=2, type='int')
    parser.add_option_group(read_filters)
    
    event_filters = OptionGroup(parser, "event_filters")
    event_filters.add_option('-n', '--min_bridge_size', dest='min_bridge_size', help='minimum size of bridge read. Default:1', 
                             default=1, type='int')
    parser.add_option_group(event_filters)
    
    combine = OptionGroup(parser, "combine")
    combine.add_option('-x', '--combine', dest='combine', help='path of discovery phase results for combining')
    combine.add_option('-k', '--track', dest='track', help='bed graph track name, desc', nargs=2)
    combine.add_option('--rgb', dest='rgb', help='RGB value of bed graph. Default:0,0,255', default='0,0,255')
    parser.add_option_group(combine)
    
    (options, args) = parser.parse_args()
    main(args, options)
        
