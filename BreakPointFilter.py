import pysam
import sys
import os
import re
import logging

#
# Class definitions
#
class BreakPointFilter:

    def __init__(self, mutationFileInfo):
        logging.info( 'create BreakPointFilter')
        self.mfi = mutationFileInfo
        self.max_depth = int(1000)
        self.min_major_clip_size = int(20)
        self.junc_num_thres = int(3)


    ############################################################
    def set_filter_info(self, mutation_line_info, max_junc_pos, max_junc_cnt_p, max_junc_cnt_m):
        if ((max_junc_cnt_p + max_junc_cnt_m) <= 0):
            mutation_line_info[2].append(0)
            mutation_line_info[3].append("0\t0\t0")

        else:
            if(max_junc_cnt_p + max_junc_cnt_m) >= self.junc_num_thres:
                mutation_line_info[2].append(1)
            else:
                mutation_line_info[2].append(0)
            
            sdist = abs(start - int(max_junc_pos))
            edist = abs(end - int(max_junc_pos))
            dist = sdist if sdist < edist else edist
            mutation_line_info[3].append(str(dist)+"\t"+ str(max_junc_cnt_p) +"\t"+ str(max_junc_cnt_m)) 


    ############################################################
    def filter(self, in_bam):
        logging.info( 'filter start')
    
        samfile = pysam.Samfile(in_bam, "rb")
        for mutation_line_info in self.mfi['line_info']:
            chr = mutation_line_info[1]['chr']
            start = mutation_line_info[1]['start']
            end = mutation_line_info[1]['end']
            
            max_junc_pos = ""
            max_junc_cnt_p = int(0)
            max_junc_cnt_m = int(0)
            
            if samfile.count(chr, start, (start+1)) >= self.max_depth:
                self.set_filter_info(mutation_line_info, max_junc_pos, max_junc_cnt_p, max_junc_cnt_m)
                continue

            bp_dict = {}

            ####
            for read in samfile.fetch(chr, max(0, int(start-(self.min_major_clip_size))), int(end+(self.min_major_clip_size))):

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
                if right_clipping >= self.min_major_clip_size:
                    
                    juncPos_current = str(int(read.pos + 1) + read.alen - 1)
                    key = juncPos_current +"\tR"
        
                    if key not in bp_dict:
                        bp_dict[key] = {"+":0, "-":0 }
                    bp_dict[key][strand] += 1

                # when the left side is clipped...
                if left_clipping >= self.min_major_clip_size:

                    juncPos_current = str(int(read.pos + 1))
                    key = juncPos_current +"\tL"

                    if key not in bp_dict:
                        bp_dict[key] = {"+":0, "-":0 }
                    bp_dict[key][strand] += 1


            for key in bp_dict:
                    
                juncPos_current = (key.split("\t")[0])
                
                if ((int(juncPos_current) - self.min_major_clip_size) <= start <= int(juncPos_current)
                or (start <= int(juncPos_current) <= end)
                or (int(juncPos_current) <= end <= int(juncPos_current) + self.min_major_clip_size)):
                
                    if (bp_dict[key]["+"] + bp_dict[key]["-"]) > (max_junc_cnt_p + max_junc_cnt_m):
                        max_junc_cnt_p = bp_dict[key]["+"]
                        max_junc_cnt_m = bp_dict[key]["-"]
                        max_junc_pos = juncPos_current

            ####
            self.set_filter_info(mutation_line_info, max_junc_pos, max_junc_cnt_p, max_junc_cnt_m)
        
        logging.info( 'filter end')

                
