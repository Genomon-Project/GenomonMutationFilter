import sys
import os
import re
import logging
import subprocess
import math
import vcf
import copy
import collections
import multiprocessing
from . import utils

#
# Class definitions
#
class Oxog_filter:

    def __init__(self, samtools_path, mpileup_params, thread_num):
        self.samtools_path = samtools_path
        self.mpileup_params = mpileup_params
        self.thread_num = thread_num
    

    def parse_bases(self, bases, qual_list, flags):

        var2num_1st = {}
        var2num_2nd = {}
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
                var = bases[0]
                if var not in var2num_1st:
                    var2num_1st[var] = 0
                if var not in var2num_2nd:
                    var2num_2nd[var] = 0

                if int(l_flags[base_ind]) & 64 == 64:
                   var2num_1st[var] = var2num_1st[var] + 1
                elif int(l_flags[base_ind]) & 128 == 128:
                   var2num_2nd[var] = var2num_2nd[var] + 1

                bases = bases[1:]
    
                if len(bases) > 0 and bases[0] in ['+', '-']:
    
                    match = re.search(r'^[\+\-](\d+)', bases)
                    indel_size = int(match.group(1))
                    bases = bases[(len(str(indel_size)) + indel_size + 1):]
                base_ind = base_ind + 1
    
        if len(qual_list) != base_ind:
            print("Error???")
            sys.exit(1)
    
        return var2num_1st, var2num_2nd
    
    
    def flag_oxog(self, ref, alt, alt_F1R2, alt_F2R1):
        oxog_flag = 0
        if ref == "C" and alt == "A":
            if int(alt_F1R2) < 2 or float(alt_F2R1) / (float(alt_F1R2) + float(alt_F2R1)) >= 0.9: oxog_flag = 1
        elif ref == "G" and alt == "T":
            if int(alt_F2R1) < 2 or float(alt_F1R2) / (float(alt_F1R2) + float(alt_F2R1)) >= 0.9: oxog_flag = 1
        return oxog_flag


    def call_mpileup(self, reg, bam_tumor, FNULL):

        #prepare mpileup params 
        m_params = self.mpileup_params.split(" ")

        d_first_pair_bases = {}
        d_second_pair_bases = {}

        # samtools mpileup 
        mpileup_cmd = [self.samtools_path, "mpileup", "-r", reg]
        mpileup_cmd.extend(m_params)
        mpileup_cmd.append(bam_tumor)

        # print mpileup_cmd
        with subprocess.Popen(mpileup_cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr = FNULL) as pileup:
            for mpileup in pileup.stdout:
                mp_list = mpileup.rstrip('\n').split('\t')
                # Prepare mpileup data
                d_first_pair_bases, d_second_pair_bases = self.parse_bases(mp_list[4], mp_list[5], mp_list[6])

        return d_first_pair_bases, d_second_pair_bases


    def filter(self, in_mutation_file, bam_tumor, output):

        with open(in_mutation_file, "r") as srcfile, open(output,'w') as hout, open(os.devnull, 'w') as FNULL:
            for line in srcfile:
                line = line.rstrip('\n')
                if line.startswith("#"):
                    print(line, file=hout)
                    continue
                elif line.startswith("Chr"):
                    print(line+"\talt_F1R2\talt_F2R1\toxog_flag", file=hout)
                    continue

                F = line.split('\t')

                if len(F[2]) == 1 and len(F[3]) == 1:
                    d_first_pair_bases, d_second_pair_bases = self.call_mpileup(f"{F[0]}:{F[1]}-{F[1]}", bam_tumor, FNULL)

                    #  F1R2 (forward 1st, reverse 2nd)
                    #  F2R1 (forward 2nd, reverse 1st)
                    f1 = d_first_pair_bases[F[3].upper()] if F[3].upper() in d_first_pair_bases else 0
                    r1 = d_first_pair_bases[F[3].lower()] if F[3].lower() in d_first_pair_bases else 0
                    f2 = d_second_pair_bases[F[3].upper()] if F[3].upper() in d_second_pair_bases else 0
                    r2 = d_second_pair_bases[F[3].lower()] if F[3].lower() in d_second_pair_bases else 0
                    f_oxog = self.flag_oxog(F[2], F[3], f1+r2, f2+r1)

                    print(f"{line}\t{f1+r2}\t{r1+f2}\t{f_oxog}", file=hout)

                else:
                    print(line+"\t\t\t",file=hout)


    def add_meta_vcf(self, vcf_reader):
        vcf_reader.formats['OF1R2'] = vcf.parser._Format('OF1R2', 1, 'Integer', "Number of F1R2 reads with ALT mismatch for oxoG")
        vcf_reader.formats['OF2R1'] = vcf.parser._Format('OF2R1', 1, 'Integer', "Number of F2R1 reads with ALT mismatch for oxoG")
        vcf_reader.infos['OXOG'] = vcf.parser._Info('OXOG', 0, 'Flag', "OxoG mutation pattern", "MutationFilter", "")


    def filter_main_vcf(self, in_mutation_file, bam_tumor, output, tumor_sample, normal_sample):

        with open(in_mutation_file, "r") as srcfile, open(output,'w') as hout, open(os.devnull, 'w') as FNULL:

            vcf_reader = vcf.Reader(srcfile)
            self.add_meta_vcf(vcf_reader)
            sample_list = vcf_reader.samples

            vcf_writer = vcf.Writer(hout, vcf_reader)

            for record in vcf_reader:
                new_record = copy.deepcopy(record)

                f_oxog = 0
                f1r2 = "."
                f2r1 = "."

                if len(record.REF) == 1 and len(str(record.ALT[0])) == 1:
                    d_first_pair_bases, d_second_pair_bases = self.call_mpileup(f"{record.CHROM}:{record.POS}-{record.POS}", bam_tumor, FNULL)

                    #  F1R2 (forward 1st, reverse 2nd)
                    #  F2R1 (forward 2nd, reverse 1st)
                    f1 = d_first_pair_bases[str(record.ALT[0]).upper()] if str(record.ALT[0]).upper() in d_first_pair_bases else 0
                    r1 = d_first_pair_bases[str(record.ALT[0]).lower()] if str(record.ALT[0]).lower() in d_first_pair_bases else 0
                    f2 = d_second_pair_bases[str(record.ALT[0]).upper()] if str(record.ALT[0]).upper() in d_second_pair_bases else 0
                    r2 = d_second_pair_bases[str(record.ALT[0]).lower()] if str(record.ALT[0]).lower() in d_second_pair_bases else 0
                    f_oxog = self.flag_oxog(record.REF, str(record.ALT[0]), f1+r2, f2+r1)
                    f1r2 = f1+r2
                    f2r1 = f2+r1

                new_record.INFO['OXOG'] = True if f_oxog == 1 else False

                # Add FPRMAT
                new_record.FORMAT = new_record.FORMAT+":OF1R2:OF2R1"
                ## tumor sample
                sx = sample_list.index(tumor_sample)
                f_keys = record.samples[sx].data._fields
                f_vals = [record.samples[sx].data[vx] for vx in range(len(f_keys))]
                handy_dict = dict(zip(f_keys, f_vals))
                handy_dict['OF1R2'] = f1r2
                handy_dict['OF2R1'] = f2r1
                new_record.samples[sx].data = collections.namedtuple('CallData', f_keys+("OF1R2","OF2R1",))
                new_vals = [handy_dict[x] for x in f_keys+("OF1R2","OF2R1",)]
                new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)
                ## normal sample
                if normal_sample != None:
                    sx = sample_list.index(normal_sample)
                    f_vals = [record.samples[sx].data[vx] for vx in range(len(f_keys))]
                    handy_dict = dict(zip(f_keys, f_vals))
                    handy_dict['OF1R2'] = "."
                    handy_dict['OF2R1'] = "."
                    new_record.samples[sx].data = collections.namedtuple('CallData', f_keys+("OF1R2","OF2R1",))
                    new_vals = [handy_dict[x] for x in f_keys+("OF1R2","OF2R1",)]
                    new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)

                vcf_writer.write_record(new_record)

        vcf_writer.close()


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


