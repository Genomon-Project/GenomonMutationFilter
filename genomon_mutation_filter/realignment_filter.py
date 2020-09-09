from __future__ import print_function
import sys
import os
import re
import pysam
import logging
import subprocess
import scipy.special
from scipy.stats import fisher_exact as fisher
import math
from . import vcf_utils
import collections
import vcf
import copy
import multiprocessing
import edlib

#
# Class definitions
#
class Realignment_filter:

    def __init__(self,referenceGenome,tumor_min_mismatch,normal_max_mismatch, search_length, score_difference, blat, header_flag, max_depth, exclude_sam_flags, thread_num, uses_blat):
        self.reference_genome = referenceGenome
        self.window = search_length
        self.score_difference = score_difference
        self.tumor_min_mismatch = tumor_min_mismatch
        self.normal_max_mismatch = normal_max_mismatch
        self.blat_cmds = [blat, '-fine']
        self.complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        self.header_flag = header_flag
        self.max_depth = max_depth
        self.exclude_sam_flags = exclude_sam_flags
        self.thread_num = thread_num
        self.uses_blat = uses_blat
     
    
    ############################################################
    def math_log_fisher_pvalue(self,fisher_pvalue):

        val = float(0.0)
        if fisher_pvalue < 10**(-60):
            val = float(60.0)
        elif fisher_pvalue  > 1.0 - 10**(-10) :
            val = float(0.0)
        else:
            val = -math.log( fisher_pvalue, 10 )
                                                                
        return val


    ############################################################
    def extractRead(self, bamfile, chr,start,end, output):

        with open(output, 'w') as hOUT:
            for read in bamfile.fetch(chr,start,end):

                # get the flag information
                read_flag = int(read.flag)

                if 0 != int(bin(self.exclude_sam_flags & read_flag),2): continue

                flags = format(read_flag, "#014b")[:1:-1]

                tempSeq = ""
                if flags[4] == "1":
                    tempSeq = "".join(self.complement.get(base) for base in reversed(str(read.seq)))
                else:
                    tempSeq = read.seq

                # the first read
                if flags[6] == "1":
                    print('>' + read.qname + '/1',file=hOUT)
                    print(tempSeq, file=hOUT)
                else:
                    print('>' + read.qname + '/2',file=hOUT)
                    print(tempSeq, file=hOUT)
        

    ############################################################
    def checkSecondBestAlignmentOriginal(self, align1, align2):

        align1tmp = list(align1)
        align1tmp.extend(align2)

        if len(align1tmp) == 1:
              return 1
        align1tmp.sort(key=lambda x:x[0],reverse=True)
        # if len(align1tmp) == 3 and (align1tmp[1][0] - align1tmp[2][0]) <= self.score_difference:
        if len(align1tmp) >= 3 and (align1tmp[1][0] - align1tmp[2][0]) <= self.score_difference:
            return 1
        return 0

    ############################################################
    def getScore(self, align):
        if len(align) >= 1:
            align.sort(key=lambda x:x[0],reverse=True)
            return (align[0][0], align[0][1])
        return (-1000,0)


    ############################################################
    def makeTwoReference(self, chr,start,end,ref,alt, output):

        seq = ""
        label = ','.join([chr, str(start), str(end), ref, alt])
        range = chr + ":" + str(int(start) - self.window + 1) +"-"+ str(int(end) + self.window)
        for item in pysam.faidx(self.reference_genome, range):
            # if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()
        seq = seq.replace('>', '')
        seq = seq.replace(range.upper(), '')

        if re.search(r'[^ACGTUWSMKRYBDHVN]', seq) is not None:
            print("The return value in get_seq function includes non-nucleotide characters:", file=sys.stderr)
            print(seq, file=sys.stderr)
            sys.exit(1)

        with open(output, 'w') as hOUT:
            print('>' + label + "_ref", file=hOUT)
            print(seq, file=hOUT)

            # for insertion
            if ref == "-":   seq = seq[0:(self.window + 1)] + alt + seq[-self.window:]
            # for deletion
            elif alt == "-": seq = seq[0:self.window] + seq[-self.window:]
             # for SNV
            else:            seq = seq[0:self.window] + alt + seq[-self.window:]

            print('>' + label + "_alt", file=hOUT)
            print(seq, file=hOUT)


    ############################################################
    def summarizeRefAlt(self, inputPsl):
        
        

        numOther = []
        numAlt = []
        numRef = []

        tempID = ""
        tempAlt = []
        tempRef = []
        ####
        with open(inputPsl, 'r') as hIN:
            for line in hIN:

                F = line.rstrip('\n').split('\t')
                if F[0].isdigit() == False: continue
    
                # remove the read pair num info ("/1", "/2") 
                if tempID != F[9]:
                    if tempID != "":
                        ####
                        if (self.checkSecondBestAlignmentOriginal(tempAlt,tempRef) == 0):
                            tempAltScore, tempAltNM = self.getScore(tempAlt)
                            tempRefScore, tempRefNM = self.getScore(tempRef)
                            if tempAltScore == tempRefScore: numOther.append(tempID[0:-2])
                            elif tempAltScore >  tempRefScore: numAlt.append(tempID[0:-2])
                            elif tempAltScore <  tempRefScore: numRef.append(tempID[0:-2])
    
                    tempID = F[9]
                    tempAlt = []
                    tempRef = []
                
                x = F[18].split(',')
                y = F[19].split(',')
                z = F[20].split(',')
                number_of_mismatch = int(F[1]) + int(F[3])
    
                for i in range(1, int(F[17])):
                
                    ly = int(y[i]) - int(y[i - 1]) - int(x[i - 1]) 
                    lz = int(z[i]) - int(z[i - 1]) - int(x[i - 1]) 
                    if (ly > 0): number_of_mismatch += ly
                    if (lz > 0): number_of_mismatch += lz
    
                my_score = int(int(F[0]) - number_of_mismatch) 
                tNM = int(F[10]) - int(F[0]) + int(F[5]) + int(F[7])
    
                if F[13][-3:] == "alt":
                    tempAlt.append((my_score, tNM))
                elif F[13][-3:] == "ref":
                    tempRef.append((my_score, tNM))

        ####
        if (len(tempAlt) > 0 and len(tempRef) > 0):
            if (self.checkSecondBestAlignmentOriginal(tempAlt,tempRef) == 0):
                tempAltScore, tempAltNM = self.getScore(tempAlt)
                tempRefScore, tempRefNM = self.getScore(tempRef)
                if tempAltScore == tempRefScore: numOther.append(tempID[0:-2])
                elif tempAltScore >  tempRefScore: numAlt.append(tempID[0:-2])
                elif tempAltScore <  tempRefScore: numRef.append(tempID[0:-2])

        return([len(set(numRef)), len(set(numAlt)), len(set(numOther))])
    
    ############################################################
    def blat_read_count(self, samfile,chr,start,end,output, thread_idx):
        ref, alt, other = 0, 0, 0
        if os.path.getsize(output + ".tmp.fa") > 0: 
            # alignment tumor short reads to the reference and alternative sequences
            FNULL = open(os.devnull, 'w')
            subprocess.check_call(self.blat_cmds + [output + ".tmp.refalt.fa", output + ".tmp.fa", output + ".tmp.psl"], 
                              stdout = FNULL, stderr = subprocess.STDOUT)
            FNULL.close()
            # summarize alignment results
            ref, alt, other = self.summarizeRefAlt(output + ".tmp.psl")
        return (ref, alt, other)

    ############################################################
    def edlib_read_count(self, samfile, chr, start, end, output):
        refalt = [None, None]
        with open(output + ".tmp.refalt.fa") as r:
            for (tname, seq) in zip(r, r):
                refalt[tname.rstrip().endswith('_alt')] = seq.rstrip()
        ref, alt = refalt

        refs = set(); alts = set(); others = set()
        with open(output + ".tmp.fa") as r:
            for (qname, seq) in zip(r, r):
                qname = re.sub(r'^>(.*)/[12]$', r'\1', qname.rstrip())
                seq = seq.rstrip()
                ref_nm = edlib.align(seq, ref, mode="HW", task="path")['editDistance']
                alt_nm = edlib.align(seq, alt, mode="HW", task="path")['editDistance']
                if ref_nm < alt_nm:
                    refs.add(qname)
                elif ref_nm > alt_nm:
                    alts.add(qname)
                else:
                    others.add(qname)
        return len(refs), len(alts), len(others)

    ############################################################
    def count_reads(self, samfile, chr, start, end, output, thread_idx):
        # extract short reads from tumor sequence data around the candidate
        self.extractRead(samfile,chr,start,end,output + ".tmp.fa")

        if self.uses_blat:
            return self.blat_read_count(samfile, chr, start, end, output, thread_idx)
        else:
            return self.edlib_read_count(samfile, chr, start, end, output)

    ############################################################
    def calc_fisher_pval(self, tumor_ref, normal_ref, tumor_alt, normal_alt):

        odds_ratio, fisher_pvalue = fisher(((int(tumor_ref),int(normal_ref)),(int(tumor_alt),int(normal_alt))), alternative='two-sided')
        log10_fisher_pvalue = '{0:.3f}'.format(float(self.math_log_fisher_pvalue(fisher_pvalue)))
        return log10_fisher_pvalue

    ############################################################
    def calc_btdtri(self, tumor_ref, tumor_alt):

        beta_01  = '{0:.3f}'.format(float(scipy.special.btdtri( int(tumor_alt) + 1, int(tumor_ref) + 1, 0.1 )))
        beta_mid = '{0:.3f}'.format(float( int(tumor_alt) + 1 ) / float( int(tumor_ref) + int(tumor_alt) + 2 ))
        beta_09  = '{0:.3f}'.format(float(scipy.special.btdtri( int(tumor_alt) + 1, int(tumor_ref) + 1, 0.9 )))
        return (beta_01, beta_mid, beta_09)


    ############################################################
    def partition_anno(self, in_mutation_file):

        # count the number of lines
        with open(in_mutation_file, "r") as hin:
            line_num = sum(1 for line in hin)

        thread_num_mod = min(line_num, self.thread_num)
        if thread_num_mod == 0: thread_num_mod = 1

        each_partition_line_num = line_num / thread_num_mod

        current_line_num = 0
        split_index = 1
        with open(in_mutation_file) as hin:
            hout = open(in_mutation_file +"."+str(split_index), "w")
            for line in hin:
                print(line.rstrip('\n'), file=hout)
                current_line_num = current_line_num + 1
                if current_line_num > each_partition_line_num and split_index < thread_num_mod:
                    current_line_num = 0
                    split_index = split_index + 1
                    hout.close()
                    hout = open(in_mutation_file +"." + str(split_index), "w")

        hout.close()
        return thread_num_mod


    ############################################################
    def Print_header(self, in_mutation_file, hResult, is_TN_pair):
        with open(in_mutation_file,'r') as srcfile:
            header = srcfile.readline().rstrip('\n')  
            if is_TN_pair:
                newheader = ("readPairNum_tumor\tvariantPairNum_tumor\totherPairNum_tumor\treadPairNum_normal\tvariantPairNum_normal\totherPairNum_normal\tP-value(fisher_realignment)")
            else:
                newheader = ("readPairNum\tvariantPairNum\totherPairNum\t10%_posterior_quantile(realignment)\tposterior_mean(realignment)\t90%_posterior_quantile(realignment)")
            print(header +"\t"+ newheader, file=hResult)


    ###########################################################
    def filter_main_pair(self, in_tumor_bam, in_normal_bam, in_mutation_file, output, thread_idx, is_multi_thread):

        seq_filename, seq_ext = os.path.splitext(in_tumor_bam)
        if seq_ext == ".cram":
            tumor_align_file = pysam.AlignmentFile(in_tumor_bam, "rc", reference_filename=self.reference_genome)
            normal_align_file = pysam.AlignmentFile(in_normal_bam, "rc", reference_filename=self.reference_genome)
        else:
            tumor_align_file = pysam.AlignmentFile(in_tumor_bam, "rb")
            normal_align_file = pysam.AlignmentFile(in_normal_bam, "rb")

        srcfile = open(in_mutation_file,'r')
        if is_multi_thread:
            hResult = open(output,'w')
        else:
            hResult = open(output,'a')

        if self.header_flag and thread_idx == 1:
            # skip header
            srcfile.readline()

        ####
        for line in srcfile:
            line = line.rstrip()
            itemlist = line.split('\t')
            # annovar input file (not zero-based number)
            chr, start, end, ref, alt  = (itemlist[0], (int(itemlist[1]) - 1), int(itemlist[2]), itemlist[3], itemlist[4])
                
            tumor_ref, tumor_alt, tumor_other, normal_ref, normal_alt, normal_other, log10_fisher_pvalue= ('---','---','---','---','---','---','---')
            if int(start) >= int(self.window):
                self.makeTwoReference(chr,start,end,ref,alt,output + ".tmp.refalt.fa")

                if tumor_align_file.count(chr,start,end) < self.max_depth:
                    tumor_ref, tumor_alt, tumor_other = self.count_reads(tumor_align_file, chr, start, end, output, thread_idx)

                if normal_align_file.count(chr,start,end) < self.max_depth:
                    normal_ref, normal_alt, normal_other = self.count_reads(normal_align_file, chr, start, end, output, thread_idx)

                if tumor_ref != '---' and  tumor_alt != '---' and  normal_ref != '---' and  normal_alt != '---':
                    log10_fisher_pvalue = self.calc_fisher_pval(tumor_ref, normal_ref, tumor_alt, normal_alt)

                if  ((tumor_alt == '---' or tumor_alt >= self.tumor_min_mismatch) and
                    (normal_alt == '---' or normal_alt <= self.normal_max_mismatch)):
                    print(line +"\t"+ str(tumor_ref)  +"\t"+ str(tumor_alt)  +"\t"+ str(tumor_other)
                               +"\t"+ str(normal_ref) +"\t"+ str(normal_alt) +"\t"+ str(normal_other)
                               +"\t"+ str(log10_fisher_pvalue), file=hResult)
        ####
        hResult.close()
        srcfile.close()
        tumor_align_file.close()
        normal_align_file.close()


    ###########################################################
    def filter_main_single(self, in_tumor_bam, in_mutation_file, output, thread_idx, is_multi_thread):

        seq_filename, seq_ext = os.path.splitext(in_tumor_bam)
        if seq_ext == ".cram":
            tumor_align_file = pysam.AlignmentFile(in_tumor_bam, "rc", reference_filename=self.reference_genome)
        else:
            tumor_align_file = pysam.AlignmentFile(in_tumor_bam, "rb")

        srcfile = open(in_mutation_file,'r')
        if is_multi_thread:
            hResult = open(output,'w')
        else:
            hResult = open(output,'a')

        if self.header_flag and thread_idx == 1:
            # skip header
            srcfile.readline()

        ####
        for line in srcfile:
            line = line.rstrip()
            itemlist = line.split('\t')
            # annovar input file (not zero-based number)
            chr, start, end, ref, alt  = (itemlist[0], (int(itemlist[1]) - 1), int(itemlist[2]), itemlist[3], itemlist[4])

            tumor_ref, tumor_alt, tumor_other, beta_01, beta_mid, beta_09 = ('---','---','---','---','---','---')
               
            if tumor_align_file.count(chr,start,end) < self.max_depth and int(start) >= int(self.window) :

                self.makeTwoReference(chr,start,end,ref,alt,output + ".tmp.refalt.fa")
                tumor_ref, tumor_alt, tumor_other = self.count_reads(tumor_align_file, chr, start, end, output, thread_idx)
                beta_01, beta_mid, beta_09 = self.calc_btdtri(tumor_ref, tumor_alt)

            if (tumor_alt == '---' or tumor_alt >= self.tumor_min_mismatch):
                print(line +"\t"+ str(tumor_ref)  +"\t"+ str(tumor_alt)
                           +"\t"+ str(tumor_other) +"\t"+ str(beta_01) 
                           +"\t"+ str(beta_mid) +"\t"+ str(beta_09), file=hResult)
            
        ####
        hResult.close()
        srcfile.close()
        tumor_align_file.close()


    ############################################################
    def filter(self, in_tumor_bam, in_normal_bam, output, in_mutation_file):

        thread_num_mod = 1
        if in_tumor_bam and in_normal_bam:

            #
            # multi thread
            #             
            if self.thread_num > 1:
                thread_num_mod = self.partition_anno(in_mutation_file)
                jobs = []
                for idx in range(1, thread_num_mod+1): 
                    proc = multiprocessing.Process(target = self.filter_main_pair, \
                        args = (in_tumor_bam, in_normal_bam, in_mutation_file +"."+ str(idx), output +"."+ str(idx), idx, True))
                    jobs.append(proc)
                    proc.start()

                for idx in range(0, thread_num_mod): 
                    jobs[idx].join() 
                    if jobs[idx].exitcode != 0:
                        raise RuntimeError('There was an error!')

                with open(output, 'w') as w:
                    if self.header_flag:
                        self.Print_header(in_mutation_file, w, True)
                    for idx in range(1, thread_num_mod+1): 
                        with open(output +"."+ str(idx), 'r') as hin:
                            for line in hin:
                                print(line.rstrip('\n'), file=w) 

            #
            # single thread
            # 
            else:
                with open(output, 'w') as w:
                    if self.header_flag:
                        self.Print_header(in_mutation_file, w, True)
                self.filter_main_pair(in_tumor_bam, in_normal_bam, in_mutation_file, output, 1, False)

        elif in_tumor_bam:

            #
            # multi thread
            #             
            if self.thread_num > 1:
                thread_num_mod = self.partition_anno(in_mutation_file)
                jobs = []
                for idx in range(1, thread_num_mod+1): 
                    proc = multiprocessing.Process(target = self.filter_main_single, \
                        args = (in_tumor_bam, in_mutation_file +"."+ str(idx), output + "." + str(idx), idx, True))
                    jobs.append(proc)
                    proc.start()

                for idx in range(0,thread_num_mod): 
                    jobs[idx].join() 
                    if jobs[idx].exitcode != 0:
                        raise RuntimeError('There was an error!')

                with open(output, 'w') as w:
                    if self.header_flag:
                        self.Print_header(in_mutation_file, w, False)
                    for idx in range(1,thread_num_mod+1): 
                        with open(output +"."+ str(idx), 'r') as hin:
                            for line in hin:
                                print(line.rstrip('\n'), file=w)

            #
            # single thread
            # 
            else:
                with open(output, 'w') as w:
                    if self.header_flag:
                        self.Print_header(in_mutation_file, w, False)
                self.filter_main_single(in_tumor_bam, in_mutation_file, output, 1, False)

        ####
        for idx in range(1, thread_num_mod+1): 
            if os.path.exists(in_mutation_file +"."+str(idx)): os.unlink(in_mutation_file +"."+str(idx))
            if os.path.exists(output +"."+str(idx)): os.unlink(output +"."+str(idx))
            if os.path.exists(output +"."+str(idx)+ ".tmp.refalt.fa"): os.unlink(output +"."+str(idx)+ ".tmp.refalt.fa")
            if os.path.exists(output +"."+str(idx)+ ".tmp.fa"): os.unlink(output +"."+str(idx)+ ".tmp.fa")
            if os.path.exists(output +"."+str(idx)+ ".tmp.psl"): os.unlink(output +"."+str(idx)+ ".tmp.psl")
        if os.path.exists(output + ".tmp.refalt.fa"): os.unlink(output + ".tmp.refalt.fa")
        if os.path.exists(output + ".tmp.fa"): os.unlink(output + ".tmp.fa")
        if os.path.exists(output + ".tmp.psl"): os.unlink(output + ".tmp.psl")


    ############################################################
    def add_meta_vcf(self, vcf_reader, is_TN_pair):
        if is_TN_pair:
            vcf_reader.infos['FPR'] = vcf.parser._Info('FPR', 1, 'Float', "Minus logarithm of the p-value by Fishers exact test processed with Realignment Filter","MutationFilter","v0.2.0")
        else:
            vcf_reader.infos['B1R'] = vcf.parser._Info('B1R', 1, 'Float', "10% posterior quantile of the beta distribution with Realignment Filter","MutationFilter","v0.2.0")
            vcf_reader.infos['BMR'] = vcf.parser._Info('BMR', 1, 'Float', "Posterior mean processed with Realignmnt Filter","MutationFilter","v0.2.0")
            vcf_reader.infos['B9R'] = vcf.parser._Info('B9R', 1, 'Float', "90% posterior quantile of the beta distribution processed with Realignment Filter","MutationFilter","v0.2.0")
        vcf_reader.formats['NNR'] = vcf.parser._Format('NNR', 1, 'Integer', "Number of non-allelic reads")
        vcf_reader.formats['NAR'] = vcf.parser._Format('NAR', 1, 'Integer', "Number of allelic reads")
        vcf_reader.formats['NOR'] = vcf.parser._Format('NOR', 1, 'Integer', "Number of other reads")
        

    ############################################################
    def partition_vcf(self, in_mutation_file):

        # count the number of lines
        line_num = 0
        hin = open(in_mutation_file, 'r')
        vcf_reader = vcf.Reader(hin)
        for record in vcf_reader:
            line_num = line_num + 1
        hin.close()
        
        thread_num_mod = min(line_num, self.thread_num)
        if thread_num_mod == 0: thread_num_mod = 1

        each_partition_line_num = line_num / thread_num_mod

        current_line_num = 0
        split_index = 1
        hout = open(in_mutation_file +"."+str(split_index), 'w')
        hin = open(in_mutation_file, 'r')
        vcf_reader = vcf.Reader(hin)
        vcf_writer = vcf.Writer(hout, vcf_reader)
        for record in vcf_reader:
            vcf_writer.write_record(record)
            current_line_num = current_line_num + 1
            if current_line_num > each_partition_line_num and split_index < thread_num_mod:
                current_line_num = 0
                split_index = split_index + 1
                vcf_writer.close()
                hout.close()
                hout = open(in_mutation_file +"."+str(split_index), 'w')
                vcf_writer = vcf.Writer(hout, vcf_reader)

        vcf_writer.close()
        hout.close()
        hin.close()
        return thread_num_mod


    ###########################################################
    def filter_main_pair_vcf(self, in_tumor_bam, in_normal_bam, in_mutation_file, output, tumor_sample, normal_sample, thread_idx):

        seq_filename, seq_ext = os.path.splitext(in_tumor_bam)
        if seq_ext == ".cram":
            tumor_align_file = pysam.AlignmentFile(in_tumor_bam, "rc",  reference_filename=self.reference_genome)
            normal_align_file = pysam.AlignmentFile(in_normal_bam, "rc",  reference_filename=self.reference_genome)
        else:
            tumor_align_file = pysam.AlignmentFile(in_tumor_bam, "rb")
            normal_align_file = pysam.AlignmentFile(in_normal_bam, "rb")

        with open(in_mutation_file, 'r') as hin:
            vcf_reader = vcf.Reader(hin)
            f_keys = vcf_reader.formats.keys() #it's an ordered dict
            len_f_keys = len(f_keys)
            self.add_meta_vcf(vcf_reader, True)
            new_keys = vcf_reader.formats.keys()
            sample_list = vcf_reader.samples
            
            hOUT = open(output, 'w')
            vcf_writer = vcf.Writer(hOUT, vcf_reader)
    
            ####
            for record in vcf_reader:
                new_record = copy.deepcopy(record)
                chr, start, end, ref, alt, is_conv = vcf_utils.vcf_fields2anno(record.CHROM, record.POS, record.REF, record.ALT[0])
                   
                tumor_ref, tumor_alt, tumor_other, normal_ref, normal_alt, normal_other, log10_fisher_pvalue= ('.','.','.','.','.','.','.')
                if int(start) >= int(self.window):
                    self.makeTwoReference(chr,start,end,ref,alt,output + ".tmp.refalt.fa")
    
                    if tumor_align_file.count(chr,start,end) < self.max_depth and is_conv:
                        tumor_ref, tumor_alt, tumor_other = self.count_reads(tumor_align_file, chr, start, end, output, thread_idx)
                    
                    if normal_align_file.count(chr,start,end) < self.max_depth and is_conv:
                        normal_ref, normal_alt, normal_other = self.count_reads(normal_align_file, chr, start, end, output, thread_idx)
    
                if tumor_ref != '.' and  tumor_alt != '.' and  normal_ref != '.' and  normal_alt != '.':
                    log10_fisher_pvalue = self.calc_fisher_pval(tumor_ref, normal_ref, tumor_alt, normal_alt)
    
                if  ((tumor_alt == '.' or tumor_alt >= self.tumor_min_mismatch) and
                    (normal_alt == '.' or normal_alt <= self.normal_max_mismatch)):
    
                    # Add INFO
                    new_record.INFO['FPR'] = float(log10_fisher_pvalue)
    
                    # Add FPRMAT
                    new_record.FORMAT = new_record.FORMAT+":NNR:NAR:NOR"
                    ## tumor sample
                    sx = sample_list.index(tumor_sample)
                    new_record.samples[sx].data = collections.namedtuple('CallData', new_keys)
                    f_vals = [record.samples[sx].data[vx] for vx in range(len_f_keys)]
                    handy_dict = dict(zip(f_keys, f_vals))
                    handy_dict['NNR'] = tumor_ref
                    handy_dict['NAR'] = tumor_alt
                    handy_dict['NOR'] = tumor_other
                    new_vals = [handy_dict[x] for x in new_keys]
                    new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)
                    ## normal sample
                    sx = sample_list.index(normal_sample)
                    new_record.samples[sx].data = collections.namedtuple('CallData', new_keys)
                    f_vals = [record.samples[sx].data[vx] for vx in range(len_f_keys)]
                    handy_dict = dict(zip(f_keys, f_vals))
                    handy_dict['NNR'] = normal_ref
                    handy_dict['NAR'] = normal_alt
                    handy_dict['NOR'] = normal_other
                    new_vals = [handy_dict[x] for x in new_keys]
                    new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)
    
                    vcf_writer.write_record(new_record)
             
        ####
        vcf_writer.close()
        hOUT.close()
        tumor_align_file.close()
        normal_align_file.close()


    ###########################################################
    def filter_main_single_vcf(self, in_tumor_bam, in_mutation_file, output, tumor_sample, thread_idx):

        seq_filename, seq_ext = os.path.splitext(in_tumor_bam)
        if seq_ext == ".cram":
            tumor_align_file = pysam.AlignmentFile(in_tumor_bam, "rc", reference_filename=self.reference_genome)
        else:
            tumor_align_file = pysam.AlignmentFile(in_tumor_bam, "rb")

        with open(in_mutation_file, 'r') as hin:
            vcf_reader = vcf.Reader(hin)
            f_keys = vcf_reader.formats.keys() #it's an ordered dict
            len_f_keys = len(f_keys)
            self.add_meta_vcf(vcf_reader, False)
            new_keys = vcf_reader.formats.keys()
            sample_list = vcf_reader.samples
    

            hOUT = open(output, 'w')
            vcf_writer = vcf.Writer(hOUT, vcf_reader)
    
            ####
            for record in vcf_reader:
                new_record = copy.deepcopy(record)
                # chr, start, end, ref, alt  = (rec.chrom (rec.pos - 1), rec.pos, rec.ref, rec.alts[0])
                chr, start, end, ref, alt, is_conv = vcf_utils.vcf_fields2anno(record.CHROM, record.POS, record.REF, record.ALT[0])
                    
                tumor_ref, tumor_alt, tumor_other, beta_01, beta_mid, beta_09 = ('','','','','','')
                   
                if tumor_align_file.count(chr,start,end) < self.max_depth and int(start) >= int(self.window) :
                    self.makeTwoReference(chr,start,end,ref,alt,output + ".tmp.refalt.fa")
                    tumor_ref, tumor_alt, tumor_other = self.count_reads(tumor_align_file, chr, start, end, output, thread_idx)
                    beta_01, beta_mid, beta_09 = self.calc_btdtri(tumor_ref, tumor_alt)
    
                if (tumor_alt == '' or tumor_alt >= self.tumor_min_mismatch):
    
                    # Add INFO
                    new_record.INFO['B1R'] = float(beta_01)
                    new_record.INFO['BMR'] = float(beta_mid)
                    new_record.INFO['B9R'] = float(beta_09)
    
                    # Add FPRMAT
                    new_record.FORMAT = new_record.FORMAT+":NNR:NAR:NOR"
                    ## tumor sample
                    sx = sample_list.index(tumor_sample)
                    new_record.samples[sx].data = collections.namedtuple('CallData', new_keys)
                    f_vals = [record.samples[sx].data[vx] for vx in range(len_f_keys)]
                    handy_dict = dict(zip(f_keys, f_vals))
                    handy_dict['NNR'] = tumor_ref
                    handy_dict['NAR'] = tumor_alt
                    handy_dict['NOR'] = tumor_other
                    new_vals = [handy_dict[x] for x in new_keys]
                    new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)
    
                    vcf_writer.write_record(new_record)
                 
        ####
        vcf_writer.close()
        hOUT.close()
        tumor_align_file.close()


    ############################################################
    def filter_vcf(self, in_tumor_bam, in_normal_bam, output, in_mutation_file, tumor_sample, normal_sample):

        thread_num_mod = 1
        if in_tumor_bam and in_normal_bam:

            #
            # multi thread
            #             
            if self.thread_num > 1:
                thread_num_mod = self.partition_vcf(in_mutation_file)
                jobs = []
                for idx in range(1, thread_num_mod+1): 
                    proc = multiprocessing.Process(target = self.filter_main_pair_vcf, \
                        args = (in_tumor_bam, in_normal_bam, in_mutation_file +"."+ str(idx), output +"."+ str(idx), tumor_sample, normal_sample, idx))
                    jobs.append(proc)
                    proc.start()

                for idx in range(0, thread_num_mod): 
                    jobs[idx].join() 
                    if jobs[idx].exitcode != 0:
                        raise RuntimeError('There was an error!')

                with open(in_mutation_file, 'r') as hin:
                    vcf_reader = vcf.Reader(hin)
                    self.add_meta_vcf(vcf_reader, True)
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
                self.filter_main_pair_vcf(in_tumor_bam, in_normal_bam, in_mutation_file, output, tumor_sample, normal_sample, 1)

        elif in_tumor_bam:

            #
            # multi thread
            #             
            if self.thread_num > 1:
                thread_num_mod = self.partition_vcf(in_mutation_file)
                jobs = []
                for idx in range(1, thread_num_mod+1): 
                    proc = multiprocessing.Process(target = self.filter_main_single_vcf, \
                        args = (in_tumor_bam, in_mutation_file +"."+ str(idx), output +"."+ str(idx), tumor_sample, idx))
                    jobs.append(proc)
                    proc.start()

                for idx in range(0, thread_num_mod): 
                    jobs[idx].join() 
                    if jobs[idx].exitcode != 0:
                        raise RuntimeError('There was an error!')

                with open(in_mutation_file, 'r') as hin:
                    vcf_reader = vcf.Reader(hin)
                    self.add_meta_vcf(vcf_reader, False)
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
                self.filter_main_single_vcf(in_tumor_bam, in_mutation_file, output, tumor_sample, 1)

        ####
        for idx in range(1, thread_num_mod+1): 
            if os.path.exists(in_mutation_file +"."+str(idx)): os.unlink(in_mutation_file +"."+str(idx))
            if os.path.exists(output +"."+str(idx)): os.unlink(output +"."+str(idx))
            if os.path.exists(output +"."+str(idx)+ ".tmp.refalt.fa"): os.unlink(output +"."+str(idx)+ ".tmp.refalt.fa")
            if os.path.exists(output +"."+str(idx)+ ".tmp.fa"): os.unlink(output +"."+str(idx)+ ".tmp.fa")
            if os.path.exists(output +"."+str(idx)+ ".tmp.psl"): os.unlink(output +"."+str(idx)+ ".tmp.psl")
        if os.path.exists(output + ".tmp.refalt.fa"): os.unlink(output + ".tmp.refalt.fa")
        if os.path.exists(output + ".tmp.fa"): os.unlink(output + ".tmp.fa")
        if os.path.exists(output + ".tmp.psl"): os.unlink(output + ".tmp.psl")
