import pysam
import sys
import os
import re
import logging
import subprocess
import vcf_utils

#
# Class definitions
#
class breakpoint_filter:

    def __init__(self, max_depth, min_clip_size, junc_num_thres, mapq_thres, header_flag, exclude_sam_flags, reference_genome):
        self.reference_genome = reference_genome
        self.max_depth = max_depth
        self.min_clip_size = min_clip_size
        self.junc_num_thres = junc_num_thres
        self.mapq_thres = mapq_thres
        self.header_flag = header_flag
        self.exclude_sam_flags = exclude_sam_flags


    def write_result_file(self, line, file_handle, dist, max_junc_cnt):
        print >> file_handle, (line +"\t"+ str(max_junc_cnt) +"\t"+str(dist)) 
        
        
    def filter_main(self, chr, start, end, samfile):

        if samfile.count(chr, start, (start+1)) >= self.max_depth:
            return ('---', '---')

        bp_dict = {}
        max_junc_pos = int(0)
        max_junc_cnt_p = int(0)
        max_junc_cnt_m = int(0)
            
        ####
        for read in samfile.fetch(chr, max(0, int(start-(self.min_clip_size))), int(end+(self.min_clip_size))):

            # get the flag information
            read_flag = int(read.flag)

            if 0 != int(bin(self.exclude_sam_flags & read_flag),2): continue

            flags = format(read_flag, "#014b")[:1:-1]

            # no clipping
            if len(read.cigar) == 1: continue
                
            # skip low mapping quality
            if read.mapq < self.mapq_thres: continue

            left_clipping = (read.cigar[0][1] if read.cigar[0][0] in [4, 5] else 0)
            right_clipping = (read.cigar[len(read.cigar) - 1][1] if read.cigar[len(read.cigar) - 1][0] in [4, 5] else 0)

            # get strand info
            strand = "-" if flags[4] == "1" else "+"

            # when the right side is clipped...
            if right_clipping >= self.min_clip_size:
                    
                juncPos_current = str(int(read.pos + 1) + read.alen - 1)
                key = juncPos_current +"\tR"
        
                if key not in bp_dict:
                    bp_dict[key] = {"+":0, "-":0 }
                bp_dict[key][strand] += 1

            # when the left side is clipped...
            if left_clipping >= self.min_clip_size:

                juncPos_current = str(int(read.pos + 1))
                juncPos_current = str(read.pos)
                key = juncPos_current +"\tL"

                if key not in bp_dict:
                    bp_dict[key] = {"+":0, "-":0 }
                bp_dict[key][strand] += 1

        dist = 0
        for key in bp_dict:
            juncPos_current = (key.split("\t")[0])

            if ((int(juncPos_current) - self.min_clip_size) <= start <= int(juncPos_current)
            or (start <= int(juncPos_current) <= end)
            or (int(juncPos_current) <= end <= int(juncPos_current) + self.min_clip_size)):
                
                if (bp_dict[key]["+"] + bp_dict[key]["-"]) > (max_junc_cnt_p + max_junc_cnt_m):
                    max_junc_cnt_p = bp_dict[key]["+"]
                    max_junc_cnt_m = bp_dict[key]["-"]
                    max_junc_pos = juncPos_current

            sdist = abs(start - int(max_junc_pos))
            edist = abs(end   - int(max_junc_pos))
            dist = sdist if sdist < edist else edist

        if (max_junc_cnt_p + max_junc_cnt_m) == 0:
            dist = 0
        return (dist, (max_junc_cnt_p + max_junc_cnt_m))

                
    def filter(self, in_mutation_file, in_bam, output):
   
        seq_filename, seq_ext = os.path.splitext(in_bam)

        if seq_ext == ".cram": 
            samfile = pysam.AlignmentFile(in_bam, "rc", reference_filename=self.reference_genome)
        else:
            samfile = pysam.AlignmentFile(in_bam, "rb")

        srcfile = open(in_mutation_file,'r')
        hResult = open(output,'w')
        if self.header_flag:
            header = srcfile.readline().rstrip('\n')  
            newheader = "bp_mismatch_count\tdistance_from_breakpoint"
            print >> hResult, (header +"\t"+ newheader)
        
        for line in srcfile:
            line = line.rstrip()
            itemlist = line.split('\t')
    
            # input file is annovar format (not zero-based number)
            chr = itemlist[0]
            start = (int(itemlist[1]) - 1)
            end = int(itemlist[2])

            dist = "---"
            junction_num = "---"
            if samfile.count(chr, start, (start+1)) < self.max_depth:
                dist, junction_num = self.filter_main(chr, start, end, samfile)

            ####
            if junction_num =="---" or  junction_num >= self.junc_num_thres:
                self.write_result_file(line, hResult, dist, junction_num)

        ####
        hResult.close()
        srcfile.close()


    def filter_vcf(self, in_mutation_file, in_bam, output, tumor_sample, normal_sample):
  
        import collections
        import vcf
        import copy

        vcf_reader = vcf.Reader(filename = in_mutation_file)
        f_keys = vcf_reader.formats.keys() #it's an ordered dict
 
        vcf_reader.formats['NB'] = vcf.parser._Format('NB', 1, 'Integer', "The number of breakpoint-containing reads around ALT by the matched normal sample") 
        vcf_reader.formats['LB'] = vcf.parser._Format('LB', 1, 'Integer', "The genome length from breakpoint-position to ALT position by matched normal sample") 
        new_keys = vcf_reader.formats.keys()
        sample_list = vcf_reader.samples

        # handle output vcf file
        vcf_writer = vcf.Writer(open(output, 'w'), vcf_reader)

        seq_filename, seq_ext = os.path.splitext(in_bam)

        if seq_ext == ".cram": 
            samfile = pysam.AlignmentFile(in_bam, "rc", reference_filename=self.reference_genome)
        else:
            samfile = pysam.AlignmentFile(in_bam, "rb")

        for record in vcf_reader:
            # input file is annovar format (not zero-based number)
            new_record = copy.deepcopy(record)
            chr, start, end, ref, alt, is_conv = vcf_utils.vcf_fields2anno(record.CHROM, record.POS, record.REF, record.ALT[0])

            dist = "."
            junction_num = "."
            if samfile.count(chr, start, (start+1)) < self.max_depth:
                dist, junction_num = self.filter_main(chr, start, end, samfile)   

            if junction_num =="." or junction_num >= self.junc_num_thres:
                # Add FPRMAT
                new_record.FORMAT = new_record.FORMAT+":NB:LB"
                ## tumor sample
                sx = sample_list.index(tumor_sample)
                new_record.samples[sx].data = collections.namedtuple('CallData', new_keys)
                f_vals = [record.samples[sx].data[vx] for vx in range(len(f_keys))]
                handy_dict = dict(zip(f_keys, f_vals))
                handy_dict['NB'] = junction_num
                handy_dict['LB'] = dist
                new_vals = [handy_dict[x] for x in new_keys]
                new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)
                ## normal sample
                sx = sample_list.index(normal_sample)
                new_record.samples[sx].data = collections.namedtuple('CallData', new_keys)
                f_vals = [record.samples[sx].data[vx] for vx in range(len(f_keys))]
                handy_dict = dict(zip(f_keys, f_vals))
                handy_dict['NB'] = "."
                handy_dict['LB'] = "."
                new_vals = [handy_dict[x] for x in new_keys]
                new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)

                vcf_writer.write_record(new_record)

        ####
        vcf_writer.close()
        samfile.close()
        
