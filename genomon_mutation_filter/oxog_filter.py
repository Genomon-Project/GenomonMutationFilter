import sys
import os
import re
import logging
import subprocess
import numpy as np
from . import utils
from scipy.stats import fisher_exact as fisher
import math
import gzip


#
# Class definitions
#
class Oxog_filter:


    def __init__(self,ref_genome, samtools_path, mpileup_params):
        self.ref_genome = ref_genome
        self.samtools_path = samtools_path
        self.mpileup_params = mpileup_params
        self.target = re.compile( '([\+\-])([0-9]+)([ACGTNRMacgtnrm]+)' )
        self.remove_chr = re.compile( '\^.' )
    

    def flag_oxog(self, ref, alt, alt_F1R2, alt_F2R1):
        oxog_flag = 0
        if ref == "C" and alt == "A":
            if int(alt_F1R2) < 2: oxog_flag = 1
            if float(alt_F2R1) / (float(alt_F1R2) + float(alt_F2R1)) >= 0.9: oxog_flag = 1
        elif ref == "G" and alt == "T":
            if int(alt_F2R1) < 2: oxog_flag = 1
            if float(alt_F1R2) / (float(alt_F1R2) + float(alt_F2R1)) >= 0.9: oxog_flag = 1
        return oxog_flag


    def proxG(self, in_x1, in_x2, in_isCtoA):

        x1 = np.array(in_x1)
        x2 = np.array(in_x2)
        isCtoA = np.array(in_isCtoA)

        theta = 0.1
        for var in range(0, 100):
            r1 = isCtoA * theta * 0.95**x1 * 0.05**x2
            r2 = (1 - theta) * ((5/float(6)) - (2/float(3)) * isCtoA) * 0.5**(x1 + x2)
            r = r1 / (r1 + r2)
            # p1 = np.log(isCtoA) + np.log(theta) + x1*np.log(0.95) + x2*np.log(0.05)
            # p2 = np.log(1 - theta) + np.log((5/float(6)) - (2/float(3)) * isCtoA) + (x1 + x2) * np.log(0.5)
            # r = 1 / (1 + np.exp(p2 - p1))
            theta = sum(r) / float(len(x1))
        score = -10 * np.log10(1-r)
        return score


    def get_fisher_pvalue(self, ref1, ref2, alt1, alt2):
        odds_ratio, fisher_pvalue = fisher(
        ((int(ref1), int(ref2)),
         (int(alt1), int(alt2))),
          alternative='two-sided'
        )
        val = float(0.0)
        if fisher_pvalue < 10**(-60):
            val = float(60.0)
        elif fisher_pvalue  > 1.0 - 10**(-10) :
            val = float(0.0)
        else:
            val = -math.log( fisher_pvalue, 10 )
        return val


    def get_ref_gene(self, chrom, start, end, FNULL):

        # samtools faidx
        reg = chrom + ":" + str(start-2) +"-"+ str(end+2) 
        faidx_cmd = [self.samtools_path, "faidx", self.ref_genome, reg]
        # print mpileup_cmd
        with subprocess.Popen(faidx_cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr = FNULL) as faidx:
            for faidx_line in faidx.stdout:
                if faidx_line.startswith(">"): continue
                ret = faidx_line.rstrip("\n")
        return ret


    def bases_format_process(self, read_bases, qual_list):

        deleted = 0
        iters = self.target.finditer( read_bases )
        for m in iters:
            site = m.start()
            type = m.group( 1 )
            num = m.group( 2 )
            bases = m.group( 3 )[ 0:int( num ) ]
            read_bases = read_bases[ 0:site - deleted ] + read_bases[ site + int( num ) + len( num ) + 1 - deleted: ]
            deleted += 1 + len( num ) + int( num )

        # Remove '^.' and '$'
        read_bases = self.remove_chr.sub( '', read_bases )
        read_bases = read_bases.replace( '$', '' )
        qual_list = qual_list.rstrip('\n')

        # Error check
        if len( read_bases ) != len( qual_list ):
            print("mpileup data is not good: [{0}], [{1}]".format( read_bases, qual_list ), file=sys.stderr)
            return None
        # Count mismatch
        return read_bases


    def set_mpileup_data(self, mp_list, d_bases1, d_bases2):

        # Prepare mpileup data
        bases_bam1 = self.bases_format_process(mp_list[4], mp_list[5])
        for base in bases_bam1:
            if base in 'ATGCatgc': d_bases1[base] += 1

        if len(mp_list) > 7:
            bases_bam2 = self.bases_format_process(mp_list[7], mp_list[8])
            for base in bases_bam2:
                if base in 'ATGCatgc': d_bases2[base] += 1


    def call_mpileup(self, reg, bam_tumor, bam_rna, FNULL):

        #prepare mpileup params 
        m_params = self.mpileup_params.split(" ")

        d_f1_bases = {"A": 0, "C": 0, "G": 0, "T": 0, "a": 0, "c": 0, "g": 0, "t": 0}
        d_f1_bases_rna = {"A": 0, "C": 0, "G": 0, "T": 0, "a": 0, "c": 0, "g": 0, "t": 0}
        d_f2_bases = {"A": 0, "C": 0, "G": 0, "T": 0, "a": 0, "c": 0, "g": 0, "t": 0}
        d_f2_bases_rna = {"A": 0, "C": 0, "G": 0, "T": 0, "a": 0, "c": 0, "g": 0, "t": 0}

        # samtools mpileup 
        mpileup_cmd = [self.samtools_path, "mpileup", "-r", reg, "--rf", "64"]
        mpileup_cmd.extend(m_params)
        l_target_bams = [bam_tumor] if bam_rna == None else [bam_tumor,bam_rna]
        mpileup_cmd.extend(l_target_bams)

        # print mpileup_cmd
        with subprocess.Popen(mpileup_cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr = FNULL) as pileup:
            for mpileup in pileup.stdout:
                mp_list = mpileup.rstrip('\n').split('\t')
                self.set_mpileup_data(mp_list, d_f1_bases, d_f1_bases_rna)

        # samtools mpileup 
        mpileup_cmd = [self.samtools_path, "mpileup", "-r", reg, "--rf", "128"]
        mpileup_cmd.extend(m_params)
        l_target_bams = [bam_tumor] if bam_rna == None else [bam_tumor,bam_rna]
        mpileup_cmd.extend(l_target_bams)

        # print mpileup_cmd
        with subprocess.Popen(mpileup_cmd, encoding='utf-8', stdout=subprocess.PIPE, stderr = FNULL) as pileup:
            for mpileup in pileup.stdout:
                mp_list = mpileup.rstrip('\n').split('\t')
                self.set_mpileup_data(mp_list, d_f2_bases, d_f2_bases_rna)

        return d_f1_bases, d_f1_bases_rna, d_f2_bases, d_f2_bases_rna


    def reverse_complement(self, seq):
        d = {'A':'T','T':'A','G':'C','C':'G', 'N':'N'}
        return "".join([d[base] for base in reversed(seq)])


    def make_record(self, d_f1_bases, d_f1_bases_rna, d_f2_bases, d_f2_bases_rna, chrom, start, end, ref, alt, bases):

        var_depth = sum(d_f1_bases.values()) + sum(d_f2_bases.values())
        var_total = d_f1_bases[alt.upper()] + d_f1_bases[alt.lower()] + d_f2_bases[alt.upper()] + d_f2_bases[alt.lower()]
        ref_f1r2 = d_f1_bases[ref.upper()] + d_f2_bases[ref.lower()]
        ref_f2r1 = d_f1_bases[ref.lower()] + d_f2_bases[ref.upper()]
        alt_f1r2 = d_f1_bases[alt.upper()] + d_f2_bases[alt.lower()]
        alt_f2r1 = d_f1_bases[alt.lower()] + d_f2_bases[alt.upper()]

        rna_depth = sum(d_f1_bases_rna.values()) + sum(d_f2_bases_rna.values())
        rna_total = d_f1_bases_rna[alt.upper()] + d_f1_bases_rna[alt.lower()] + d_f2_bases_rna[alt.upper()] + d_f2_bases_rna[alt.lower()]

        s1 = (alt_f1r2 / float(var_total)) if var_total > 0 else 0.0
        s2 = (alt_f2r1 / float(var_total)) if var_total > 0 else 0.0
        var_mis   = '{0:.3f}'.format(var_total / float(var_depth))

        read_ratio = -1
        if (ref == "C" or ref == "A"):
            r_ref = ref
            r_alt = alt
            r_bases = bases
            if (s1 + s2) > 0:
                read_ratio = (s2 / (s1 + s2))
            fisher_pval = self.get_fisher_pvalue(ref_f2r1, ref_f1r2, alt_f2r1, alt_f1r2)

        elif (ref == "G" or ref == "T"):
            r_ref = self.reverse_complement(ref)
            r_alt = self.reverse_complement(alt)
            r_bases = self.reverse_complement(bases)
            if (s1 + s2) > 0:
                read_ratio = (s1 / (s1 + s2))
            fisher_pval = self.get_fisher_pvalue(ref_f1r2, ref_f2r1, alt_f1r2, alt_f2r1)

        f_oxog = self.flag_oxog(ref, alt, alt_f1r2, alt_f2r1)

        outstr = (
        '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\t{16}\t{17}\t{18}\t{19}'.format(
        chrom, start, end, ref, alt,
        bases, bases[1:2], bases[3:4],
        r_ref, r_alt,
        r_bases, r_bases[1:2], r_bases[3:4],
        var_depth, var_total, var_mis,
        ref_f1r2, ref_f2r1, alt_f1r2, alt_f2r1)
        )
        if read_ratio >= 0:
            outstr = outstr + '\t{0:.3f}'.format(read_ratio)
        else:
            outstr = outstr + "\t---"
        outstr = outstr + '\t{0:.3f}'.format(fisher_pval)
        outstr = outstr + '\t{0}\t{1}'.format(rna_depth, rna_total)
        if rna_depth > 0:
            outstr = outstr + '\t{0:.3f}'.format(rna_total/float(rna_depth))
        else:
            outstr = outstr + "\t---"
        outstr = outstr + "\t" + str(f_oxog)
        return outstr


    def print_anno(self, output):

        x1list = []
        x2list = []
        isCtoAlist = []
        with open(output+".tmp", "r") as h_tmp:
            for line in h_tmp:
                F = line.rstrip('\n').split('\t')
                ref = F[3] # C
                ref_RC = F[8] # A or C
                alt_RC = F[9] # A or C
                alt_F1R2 = int(F[18])
                alt_F2R1 = int(F[19])

                if ref == "C" or ref == "A":
                    x1list.append(alt_F2R1)
                    x2list.append(alt_F1R2)
                elif ref == "G" or ref == "T":
                    x1list.append(alt_F1R2)
                    x2list.append(alt_F2R1)
                else:
                    continue

                if ref_RC == "C" and alt_RC == "A":
                    isCtoAlist.append(1)
                else:
                    isCtoAlist.append(0)
                   
            if len(x1list) != 0:
                score = self.proxG(x1list,x2list,isCtoAlist)


        with open(output,'w') as hResult, open(output+".tmp", "r") as h_tmp:

            newheader = ("chr\tstart\tend"
                     + "\tref\talt\ttwo5and3bases\t5base\t3base"
                     + "\tref_RC\talt_RC\ttwo5and3bases_RC\t5base_RC\t3base_RC"
                     + "\tdepth\tvar\tmisMatch"
                     + "\tref_F1R2\tref_F2R1\talt_F1R2\talt_F2R1"
                     + "\tread_ratio\tfisher_pval"
                     + "\trna_depth\trna_var\trna_misRate\toxog_flag"
                     + "\tscore(oxoG)")
            print(newheader, file=hResult)
    
            num = 0
            for line in h_tmp:
                line = line.rstrip('\n')
                F = line.split('\t')
                if F[3] in ['A','C','G','T']:
                    print(line +'\t{0:.3f}'.format(score[num]), file=hResult)
                    num += 1
                else:
                    print(line, file=hResult)
            

    def prepare_mpileup_params(self, chrom, start, end, alt): 

        ret = ""
        if alt == "-":
            ret = chrom + ":" + str(start-1) +"-"+ str(start-1) 
        else:
            ret = chrom + ":" + str(start) +"-"+ str(end) 

        return ret


    def filter(self, in_mutation_file, bam_tumor, output, bam_rna):

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

                # annovar input file (not zero-based number)
                chrom,start,end,ref,alt, is_conv = utils.vcf_fields2anno(F[0], int(F[1]), F[2], F[3]) 

                if (ref == "-"  or  alt == "-" or is_conv == False):
                    print(line+"\t\t\t",file=hout)

                else:
                    genes = self.get_ref_gene(chrom, start, end, FNULL)

                    reg = self.prepare_mpileup_params(chrom, start, end, alt) 

                    d_f1_bases, d_f1_bases_rna, d_f2_bases, d_f2_bases_rna = self.call_mpileup(reg, bam_tumor, bam_rna, FNULL)

                    record = self.make_record(d_f1_bases, d_f1_bases_rna, d_f2_bases, d_f2_bases_rna, chrom, start, end, ref, alt, genes)

                    l_record = record.split("\t")
                    print(line+"\t"+l_record[18]+"\t"+l_record[19]+"\t"+l_record[25], file=hout)

        # self.print_anno(output)


    def filter_vcf(self, in_mutation_file, bam_tumor, output, bam_rna):

        return None

        '''
        import collections
        import vcf
        import copy

        ####
        with open(in_mutation_file, "r") as srcfile, open(output+".tmp",'w') as h_tmp, open(os.devnull, 'w') as FNULL:
            vcf_reader = vcf.Reader(srcfile)
            f_keys = vcf_reader.formats.keys() #its an ordered dict
            len_f_keys = len(f_keys)

            for record in vcf_reader:
                new_record = copy.deepcopy(record)
                chrom,start,end,ref,alt, is_conv = utils.vcf_fields2anno(record.CHROM,record.POS,record.REF,record.ALT[0])

                if (ref == "-"  or  alt == "-" or is_conv == False): continue
    
                genes = self.get_ref_gene(chrom, start, end, FNULL)

                d_f1_bases, d_f1_bases_rna, d_f2_bases, d_f2_bases_rna = self.call_mpileup(chrom, start, end, bam_tumor, bam_rna, FNULL)

                record = self.make_record(d_f1_bases, d_f1_bases_rna, d_f2_bases, d_f2_bases_rna, chrom, start, end, ref, alt, genes)

                print(record, file=h_tmp)

        self.print_anno(output)
        '''



