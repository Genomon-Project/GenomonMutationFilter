#! /usr/bin/env python

import re, sys, math, pysam
import os
import subprocess
import numpy as np
import vcf
import copy
import collections
import multiprocessing
from . import utils


#
# Class definitions
#
class Position_filter:

    def __init__(self,ref_genome, samtools_path, mpileup_params, thread_num):
        self.ref_genome = ref_genome
        self.samtools_path = samtools_path
        self.mpileup_params = mpileup_params
        self.thread_num = thread_num

    def parse_bases(self, bases, positions, qnames, flags, ref):

        var2num = {}
        var2pos = {}
        var2qname = {}
        var2flag = {}
    
        l_positions = positions.split(',')
        l_qnames = qnames.split(',')
        l_flags = flags.split(',')
        base_ind = 0
    
        while bases != '':
            if bases[0] in ['>', '<', '*']: 
                base_ind = base_ind + 1
                bases = bases[1:]
    
            elif bases[0] in '^':
                bases = bases[2:]
            elif bases[0] in '$':
                bases = bases[1:]
            elif bases[0] in ['.', ',', 'A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n']:
                var_original = bases[0]
                if var_original in ['.', ',']: 
                    var_original = ref
                var = var_original.upper()
                if var not in var2num:
                    var2num[var], var2pos[var], var2qname[var], var2flag[var]  = 0, [], [], []
                var2num[var] = var2num[var] + 1
                var2pos[var].append(l_positions[base_ind])
                var2qname[var].append(l_qnames[base_ind])
                var2flag[var].append(l_flags[base_ind])
    
                bases = bases[1:]
    
                if len(bases) > 0 and bases[0] in ['+', '-']:
    
                    match = re.search(r'^[\+\-](\d+)', bases)
                    indel_size = int(match.group(1))
                    var_original = bases[0] + bases[(len(str(indel_size)) + 1):(len(str(indel_size)) + indel_size + 1)]
                    var = var_original.upper()
                    if var not in var2num:
                        var2num[var], var2pos[var], var2qname[var], var2flag[var]  = 0, [], [], []
                    var2num[var] = var2num[var] + 1
                    var2pos[var].append(l_positions[base_ind])
                    var2qname[var].append(l_qnames[base_ind])
                    var2flag[var].append(l_flags[base_ind])
    
                    bases = bases[(len(str(indel_size)) + indel_size + 1):]
                base_ind = base_ind + 1
    
        if len(l_positions) != base_ind:
            print("Error???")
            sys.exit(1)
    
        return var2num, var2pos, var2qname, var2flag
    
    
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

     
    def pysam_fetch(self, chrom, pos1, pos2, samfile):

        d_ret = {}

        for read in samfile.fetch(chrom,int(pos1)-1,int(pos2)):
            d_ret[read.qname +"\t"+ str(read.flag)] = (read.cigar,read.query_length,read.tags)

        return d_ret


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
    
        # for MNV (same processing as SNV)
        elif len(ref) > 1 and len(alt) > 1 and len(ref) == len(alt):
            ret = alt[0]

        # for block substitution
        else:
            ret = None
    
        return ret


    def get_cigar_size(self, cigar):

        # softclip:4, hardclip:5
        cigar_left = 0
        cigar_right = 0

        if cigar[0][0] == 4 or cigar[0][0] == 5:
            cigar_left=cigar[0][1]
        if cigar[-1][0] == 4 or cigar[-1][0] == 5:
            cigar_right=cigar[-1][1]

        return cigar_left, cigar_right

    
    def get_nm(self, tags):

        ret = 0
        for tag, val in tags:
            if tag == "NM":
                ret = val
        return ret


    def filter(self, in_mutation_file, bam, output):

        pysam_file = pysam.AlignmentFile(bam)
    
        with open(in_mutation_file, "r") as srcfile, open(output,'w') as hout, open(os.devnull, 'w') as FNULL:
            for line in srcfile:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    print(line, file=hout)
                    continue
                elif line.startswith("Chr"):
                    print(line+"\tleft_read_position_mean\tleft_read_positon_sd\tright_read_position_mean\tright_read_position_sd\tNM_mean_without_ALT_len", file=hout)
                    continue

                F = line.split('\t')

                # annovar input file (not zero-based number)
                pileup_key = self.get_alt_pileup_key(F[2], F[3]) 
                nm = abs(len(F[2]) - len(F[3])) if len(F[2]) != len(F[3]) else len(F[3])

                left_mean = ""
                left_std = ""
                right_mean = ""
                right_std = ""
                nm_without_alt = ""

                # block substitution not suppport
                if pileup_key != None:
                
                    l_mp = self.call_mpileup(f"{F[0]}:{F[1]}-{F[1]}", bam, FNULL)
                    var2num, var2pos, var2qname, var2flag  = self.parse_bases(l_mp[4], l_mp[6], l_mp[7], l_mp[8], l_mp[2])
                    if pileup_key in var2num:

                        l_left_position = []
                        l_right_position = []
                        l_alt_mismatch = []

                        d_qname_pysam = self.pysam_fetch(F[0], F[1], F[1], pysam_file)

                        for qname, mp_flag, mut_position in zip(var2qname[pileup_key], var2flag[pileup_key], var2pos[pileup_key]):

                            if qname +"\t"+ mp_flag not in d_qname_pysam:
                                continue 

                            cigar, query_length, tags = d_qname_pysam[qname +"\t"+ mp_flag]
                            cigar_left, cigar_right = self.get_cigar_size(cigar)

                            l_left_position.append(int(mut_position) - int(cigar_left))
                            l_right_position.append(int(query_length) - int(cigar_right) - int(mut_position) + 1)
                            l_alt_mismatch.append(int(self.get_nm(tags)))

                        fct = 1000
                        left_mean = math.floor(np.average(l_left_position) * fct) / fct
                        left_std = math.floor(np.std(l_left_position) * fct) / fct
                        right_mean = math.floor(np.average(l_right_position) * fct) / fct
                        right_std = math.floor(np.std(l_right_position) * fct) / fct
                        nm_without_alt =  math.floor((np.average(l_alt_mismatch) - float(nm))  * fct) / fct

                print(line+"\t"+str(left_mean)+"\t"+str(left_std)+"\t"+str(right_mean)+"\t"+str(right_std)+"\t"+str(nm_without_alt), file=hout)

        pysam_file.close()


    def add_meta_vcf(self, vcf_reader):
        vcf_reader.formats['LPM'] = vcf.parser._Format('LPM', 1, 'Float', "Mean position of the mismatches in bases starting from left end")
        vcf_reader.formats['LPS'] = vcf.parser._Format('LPS', 1, 'Float', "Standard devition position of the mismatches in bases starting from left end")
        vcf_reader.formats['RPM'] = vcf.parser._Format('RPM', 1, 'Float', "Mean position of the mismatches in bases starting from right end")
        vcf_reader.formats['RPS'] = vcf.parser._Format('RPS', 1, 'Float', "Standard devition position of the mismatches in bases starting from right end")
        vcf_reader.formats['NMA'] = vcf.parser._Format('NMA', 1, 'Float', "Mean number of the mismathes without ALT length")


    def filter_main_vcf(self, in_mutation_file, bam_tumor, output, tumor_sample, normal_sample):

        pysam_file = pysam.AlignmentFile(bam_tumor)
    
        with open(in_mutation_file, "r") as srcfile, open(output,'w') as hout, open(os.devnull, 'w') as FNULL:

            vcf_reader = vcf.Reader(srcfile)
            self.add_meta_vcf(vcf_reader)
            sample_list = vcf_reader.samples

            vcf_writer = vcf.Writer(hout, vcf_reader)

            for record in vcf_reader:
                new_record = copy.deepcopy(record)

                # annovar input file (not zero-based number)
                pileup_key = self.get_alt_pileup_key(record.REF, str(record.ALT[0])) 
                nm = abs(len(record.REF) - len(record.ALT[0])) if len(record.REF) != len(record.ALT[0]) else len(record.ALT[0])

                left_mean = "."
                left_std = "."
                right_mean = "."
                right_std = "."
                nm_without_alt = "."

                # block substitution not suppport
                if pileup_key != None:

                    l_mp = self.call_mpileup(f"{record.CHROM}:{record.POS}-{record.POS}", bam_tumor, FNULL)
                    var2num, var2pos, var2qname, var2flag  = self.parse_bases(l_mp[4], l_mp[6], l_mp[7], l_mp[8], l_mp[2])
                    if pileup_key in var2num:
                    
                        l_left_position = []
                        l_right_position = []
                        l_alt_mismatch = []

                        d_qname_pysam = self.pysam_fetch(record.CHROM, record.POS, record.POS, pysam_file)

                        for qname, mp_flag, mut_position in zip(var2qname[pileup_key], var2flag[pileup_key], var2pos[pileup_key]):

                            if qname +"\t"+ mp_flag not in d_qname_pysam:
                                continue 

                            cigar, query_length, tags = d_qname_pysam[qname +"\t"+ mp_flag]
                            cigar_left, cigar_right = self.get_cigar_size(cigar)

                            l_left_position.append(int(mut_position) - int(cigar_left))
                            l_right_position.append(int(query_length) - int(cigar_right) - int(mut_position) + 1)
                            l_alt_mismatch.append(int(self.get_nm(tags)))
                       
                        fct = 1000
                        left_mean = math.floor(np.average(l_left_position) * fct) / fct
                        left_std = math.floor(np.std(l_left_position) * fct) / fct
                        right_mean = math.floor(np.average(l_right_position) * fct) / fct
                        right_std = math.floor(np.std(l_right_position) * fct) / fct
                        nm_without_alt =  math.floor((np.average(l_alt_mismatch) - float(nm))  * fct) / fct

                # Add FPRMAT
                new_record.FORMAT = new_record.FORMAT+":LPM:LPS:RPM:RPS:NMA"
                ## tumor sample
                sx = sample_list.index(tumor_sample)
                f_keys = record.samples[sx].data._fields
                f_vals = [record.samples[sx].data[vx] for vx in range(len(f_keys))]
                handy_dict = dict(zip(f_keys, f_vals))
                handy_dict['LPM'] = left_mean
                handy_dict['LPS'] = left_std
                handy_dict['RPM'] = right_mean
                handy_dict['RPS'] = right_std
                handy_dict['NMA'] = nm_without_alt
                new_record.samples[sx].data = collections.namedtuple('CallData', f_keys+("LPM","LPS","RPM","RPS","NMA",))
                new_vals = [handy_dict[x] for x in f_keys+("LPM","LPS","RPM","RPS","NMA",)]
                new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)
                ## normal sample
                if normal_sample != None:
                    sx = sample_list.index(normal_sample)
                    f_vals = [record.samples[sx].data[vx] for vx in range(len(f_keys))]
                    handy_dict = dict(zip(f_keys, f_vals))
                    handy_dict['LPM'] = "."
                    handy_dict['LPS'] = "."
                    handy_dict['RPM'] = "."
                    handy_dict['RPS'] = "."
                    handy_dict['NMA'] = "."
                    new_record.samples[sx].data = collections.namedtuple('CallData', f_keys+("LPM","LPS","RPM","RPS","NMA",))
                    new_vals = [handy_dict[x] for x in f_keys+("LPM","LPS","RPM","RPS","NMA",)]
                    new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)

                vcf_writer.write_record(new_record)

        pysam_file.close()


    def filter_vcf(self, in_mutation_file, bam_tumor, output, tumor_sample, normal_sample):

        thread_num_mod = 1
        #
        # multi thread
        #             
        if self.thread_num > 1:
            thread_num_mod = utils.partition_vcf(in_mutation_file, self.thread_num)
            jobs = []
            for idx in range(1, thread_num_mod+1): 
                proc = multiprocessing.Process(target = self.filter_main_vcf, \
                    args = (in_mutation_file +"."+ str(idx), bam_tumor, output +"."+ str(idx), tumor_sample, normal_sample))
                jobs.append(proc)
                proc.start()

            for idx in range(0, thread_num_mod): 
                jobs[idx].join() 
                if jobs[idx].exitcode != 0:
                    raise RuntimeError('There was an error!')

            with open(in_mutation_file, 'r') as hin:
                vcf_reader = vcf.Reader(hin)
                self.add_meta_vcf(vcf_reader)
                with open(output, 'w') as hout:
                    vcf_writer = vcf.Writer(hout, vcf_reader)
                    for idx in range(1, thread_num_mod+1): 
                        with open(output +"."+ str(idx), 'r') as hin_tmp:
                            vcf_reader_tmp = vcf.Reader(hin_tmp)
                            for record in vcf_reader_tmp:
                                vcf_writer.write_record(record)
                vcf_writer.close()

        #
        # single thread
        # 
        else:
            self.filter_main_vcf(in_mutation_file, bam_tumor, output, tumor_sample, normal_sample)

        ####
        for idx in range(1, thread_num_mod+1): 
            if os.path.exists(in_mutation_file +"."+str(idx)): os.unlink(in_mutation_file +"."+str(idx))
            if os.path.exists(output +"."+str(idx)): os.unlink(output +"."+str(idx))


        return None

