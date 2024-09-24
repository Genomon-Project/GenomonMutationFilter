import sys
import os
import re
import logging
import subprocess
import math


#
# Class definitions
#
class Oxog_filter:


    def __init__(self, samtools_path, mpileup_params):
        self.samtools_path = samtools_path
        self.mpileup_params = mpileup_params
    

    def parse_bases(self, bases, qual_list):

        var2num = {}
        var2pos = {}
    
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
                if var not in var2num:
                    var2num[var] = 0
                var2num[var] = var2num[var] + 1
    
                bases = bases[1:]
    
                if len(bases) > 0 and bases[0] in ['+', '-']:
    
                    match = re.search(r'^[\+\-](\d+)', bases)
                    indel_size = int(match.group(1))
                    bases = bases[(len(str(indel_size)) + indel_size + 1):]
                base_ind = base_ind + 1
    
        if len(qual_list) != base_ind:
            print("Error???")
            sys.exit(1)
    
        return var2num
    
    
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
        mpileup_cmd = [self.samtools_path, "mpileup", "-r", reg, "--rf", "64"]
        mpileup_cmd.extend(m_params)
        mpileup_cmd.append(bam_tumor)

        # print mpileup_cmd
        with subprocess.Popen(mpileup_cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr = FNULL) as pileup:
            for mpileup in pileup.stdout:
                mp_list = mpileup.rstrip('\n').split('\t')
                # Prepare mpileup data
                d_first_pair_bases = self.parse_bases(mp_list[4], mp_list[5])

        # samtools mpileup 
        mpileup_cmd = [self.samtools_path, "mpileup", "-r", reg, "--rf", "128"]
        mpileup_cmd.extend(m_params)
        mpileup_cmd.append(bam_tumor)

        # print mpileup_cmd
        with subprocess.Popen(mpileup_cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr = FNULL) as pileup:
            for mpileup in pileup.stdout:
                mp_list = mpileup.rstrip('\n').split('\t')
                # Prepare mpileup data
                d_second_pair_bases = self.parse_bases(mp_list[4], mp_list[5])
        
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



    def filter_vcf(self, in_mutation_file, bam_tumor, output):

        return None

