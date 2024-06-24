#! /usr/bin/env python

import re, sys, math, pysam
import os
from scipy import stats
from . import utils
import subprocess
import numpy as np


#
# Class definitions
#
class Position_filter:

    def __init__(self,ref_genome, samtools_path, mpileup_params):
        self.ref_genome = ref_genome
        self.samtools_path = samtools_path
        self.mpileup_params = mpileup_params
        self.target = re.compile( '([\+\-])([0-9]+)([ACGTNRMacgtnrm]+)' )
        self.remove_chr = re.compile( '\^.' )

 
    def parse_bases(self, bases, positions, qnames):

        var2num = {}
        var2pos = {}
        var2num_plus = {}
        var2qname = {}
    
        l_positions = positions.split(',')
        l_qnames = qnames.split(',')
        base_ind = 0
        depth_p, depth_n = 0, 0
    
        while bases != '':
            if bases[0] in ['>', '<', '*']: 
                base_ind = base_ind + 1
                bases = bases[1:]
    
            elif bases[0] in '^':
                bases = bases[2:]
            elif bases[0] in '$':
                bases = bases[1:]
            elif bases[0] in ['.', ',', 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n']:
                if bases[0] not in ['.', ',']: 
                    var_original = bases[0]
                    var = var_original.upper()
                    if var not in var2num:
                        var2num[var], var2pos[var], var2num_plus[var], var2qname[var] = 0, [], 0, []
                    var2num[var] = var2num[var] + 1
                    var2pos[var].append(l_positions[base_ind])
                    var2qname[var].append(l_qnames[base_ind])
                    if var == var_original: 
                        var2num_plus[var] = var2num_plus[var] + 1
    
                if bases[0] in ['.', 'A', 'C', 'G', 'T', 'N']:
                    depth_p = depth_p + 1
                else:
                    depth_n = depth_n + 1
    
                bases = bases[1:]
    
                if len(bases) > 0 and bases[0] in ['+', '-']:
    
                    match = re.search(r'^[\+\-](\d+)', bases)
                    indel_size = int(match.group(1))
                    var_original = bases[0] + bases[(len(str(indel_size)) + 1):(len(str(indel_size)) + indel_size + 1)]
                    var = var_original.upper()
                    if var not in var2num:
                        var2num[var], var2pos[var], var2num_plus[var], var2qname[var] = 0, [], 0, []
                    var2num[var] = var2num[var] + 1
                    var2pos[var].append(l_positions[base_ind])
                    var2qname[var].append(l_qnames[base_ind])
                    if var == var_original: var2num_plus[var] = var2num_plus[var] + 1
    
                    bases = bases[(len(str(indel_size)) + indel_size + 1):]
                base_ind = base_ind + 1
    
        if len(l_positions) != base_ind:
            print("Error???")
            sys.exit(1)
    
        return depth_p, depth_n, var2num, var2pos, var2qname
    
    
    def call_mpileup(self, reg, bam, FNULL):

        l_ret = None

        # samtools mpileup 
        mpileup_cmd = [self.samtools_path, "mpileup", "-r", reg, "-f", self.ref_genome]
        mpileup_cmd.extend(self.mpileup_params.split(" "))
        mpileup_cmd.extend([bam])

        # print mpileup_cmd
        with subprocess.Popen(mpileup_cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr = FNULL) as pileup:
            for mpileup in pileup.stdout:
                l_ret = mpileup.rstrip('\n').split('\t')

        return l_ret

     
    def prepare_mpileup_params(self, chrom, start, end, alt): 

        ret = ""
        if alt == "-":
            ret = chrom + ":" + str(start-1) +"-"+ str(start-1) 
        else:
            ret = chrom + ":" + str(start) +"-"+ str(end) 

        return ret


    def pysam_fetch(self, chrom, pos1, pos2, bam):

        d_ret = {}
        samfile = pysam.AlignmentFile(bam)

        for read in samfile.fetch(chrom,pos1,pos2):
            d_ret[read.qname] = (read.cigar,read.query_length)

        return d_ret


    def prepare_pysam_params(self, chrom, start, end, alt): 

        pos1, pos2 = None, None

        if alt == "-":
            pos1, pos2 = start-2, start-1
        else:
            pos1, pos2 = start-1, end

        return pos1, pos2


    def get_alt_pileup_key(self, ref, alt):
    
        ret = None

        # for insertion
        if len(ref) < len(alt) and len(ref) == 1 and alt[0:1] == ref:
            ret = "+"+alt[1:]
    
        # for deletion
        elif len(ref) > len(alt) and len(alt) == 1 and ref[0:1] == alt:
            ret = "-"+ref[1:]
    
        # for SNV
        elif len(ref) == 1 and len(alt) == 1:
            ret = alt
    
        # for block substitution
        else:
            ret = None
    
        return ret


    def get_cigar_size(self, cigar):

        cigar_left = 0
        cigar_right = 0

        if cigar[0][0] == 4:
            cigar_left=cigar[0][1]
        elif cigar[-1][0] == 4:
            cigar_right=cigar[-1][1]

        return cigar_left, cigar_right

    
    def filter(self, in_mutation_file, bam, output):
    
        with open(in_mutation_file, "r") as srcfile, open(output,'w') as hout, open(os.devnull, 'w') as FNULL:
            for line in srcfile:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    print(line, file=hout)
                    continue
                elif line.startswith("Chr"):
                    print(line+"\tleft_read_position_mean\tleft_read_positon_sd\tright_read_position_mean\tright_read_position_sd", file=hout)
                    continue

                F = line.split('\t')

                # annovar input file (not zero-based number)
                chrom,start,end,ref,alt, is_conv = utils.vcf_fields2anno(F[0], int(F[1]), F[2], F[3]) 
                pileup_key = self.get_alt_pileup_key(F[2], F[3]) 

                l_left_position = []
                l_right_position = []

                # block substitution not suppport
                if pileup_key != None:

                    reg = self.prepare_mpileup_params(chrom, start, end, alt) 
                    l_mp = self.call_mpileup(reg, bam, FNULL)
                    depth_p, depth_n, var2num, var2pos, var2qname = self.parse_bases(l_mp[4], l_mp[6], l_mp[7])

                    pos1, pos2 = self.prepare_pysam_params(chrom, start, end, alt) 
                    d_qname_pysam = self.pysam_fetch(chrom, pos1, pos2, bam)

                    for idx, qname in enumerate(var2qname[pileup_key]):
                        mut_position = var2pos[pileup_key][idx]
                        cigar, query_length = d_qname_pysam[qname]

                        cigar_left, cigar_right = self.get_cigar_size(cigar)

                        l_left_position.append(int(mut_position) - int(cigar_left))
                        l_right_position.append(int(query_length) - int(cigar_right) - int(mut_position) + 1)

                    left_mean = math.floor(np.average(l_left_position) * 10000) / 10000
                    left_std = math.floor(np.std(l_left_position) * 10000) / 10000
                    right_mean = math.floor(np.average(l_right_position) * 10000) / 10000
                    right_std = math.floor(np.std(l_right_position) * 10000) / 10000
                print(line+"\t"+str(left_mean)+"\t"+str(left_std)+"\t"+str(right_mean)+"\t"+str(right_std), file=hout)

    def filter_vcf(self, bam, target_file):
        return None

