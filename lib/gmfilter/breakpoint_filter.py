import pysam
import sys
import os
import re
import logging

#
# Class definitions
#
class breakpoint_filter:

    def __init__(self, max_depth, min_clip_size, junc_num_thres):
        self.max_depth = max_depth
        self.min_clip_size = min_clip_size
        self.junc_num_thres = junc_num_thres


    ############################################################
    def filter(self, in_mutation_file, in_bam, output):
        logging.info( 'filter start')
    
        samfile = pysam.Samfile(in_bam, "rb")

        chrIndex = 0
        srcfile = open(in_mutation_file,'r')
        header = srcfile.readline()  
        headerlist = header.split('\t')
        for colname in headerlist:
            if (colname == 'Chr'): 
                break
            chrIndex += 1
        
        for line in srcfile:
            line = line.rstrip()
            itemlist = line.split('\t')
    
            # input file is annovar format (not zero-based number)
            chr = itemlist[chrIndex]
            start = (int(itemlist[chrIndex + 1]) - 1)
            end = int(itemlist[chrIndex + 2])
            
            max_junc_pos = ""
            max_junc_cnt_p = int(0)
            max_junc_cnt_m = int(0)
            
            if samfile.count(chr, start, (start+1)) >= self.max_depth:
                print line +"\t\t0\t0"
                continue

            bp_dict = {}

            ####
            for read in samfile.fetch(chr, max(0, int(start-(self.min_clip_size))), int(end+(self.min_clip_size))):

                # get the flag information
                flags = format(int(read.flag), "#014b")[:1:-1]

                # skip read unmapped
                if flags[2] == "1": continue

                # skip supplementary alignment
                if flags[8] == "1" or flags[11] == "1": continue

                # skip duplicated reads
                if flags[10] == "1": continue

                # no clipping
                if len(read.cigar) == 1: continue
                
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
                    key = juncPos_current +"\tL"

                    if key not in bp_dict:
                        bp_dict[key] = {"+":0, "-":0 }
                    bp_dict[key][strand] += 1


            for key in bp_dict:
                juncPos_current = (key.split("\t")[0])
                
                if ((int(juncPos_current) - self.min_clip_size) <= start <= int(juncPos_current)
                or (start <= int(juncPos_current) <= end)
                or (int(juncPos_current) <= end <= int(juncPos_current) + self.min_clip_size)):
                
                    if (bp_dict[key]["+"] + bp_dict[key]["-"]) > (max_junc_cnt_p + max_junc_cnt_m):
                        max_junc_cnt_p = bp_dict[key]["+"]
                        max_junc_cnt_m = bp_dict[key]["-"]
                        max_junc_pos = juncPos_current

            ####
            if(max_junc_cnt_p + max_junc_cnt_m) >= self.junc_num_thres:
                sdist = abs(start - int(max_junc_pos))
                edist = abs(end - int(max_junc_pos))
                dist = sdist if sdist < edist else edist
                print line +"\t"+ str(dist) +"\t"+ str(max_junc_cnt_p) +"\t"+ str(max_junc_cnt_m)) 
        
                
