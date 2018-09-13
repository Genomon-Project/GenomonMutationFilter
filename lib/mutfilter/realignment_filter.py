import sys
import os
import re
import pysam
import logging
import subprocess
import scipy.special
from scipy.stats import fisher_exact as fisher
import math
import vcf_utils

#
# Class definitions
#
class realignment_filter:

    def __init__(self,referenceGenome,tumor_min_mismatch,normal_max_mismatch, search_length, score_difference, blat, header_flag, max_depth, exclude_sam_flags):
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
    def makeTwoReference(self, chr,start,end,ref,alt, output):

        hOUT = open(output, 'w')
        
        seq = ""
        label = ','.join([chr, str(start), str(end), ref, alt])
        range = chr + ":" + str(int(start) - self.window + 1) +"-"+ str(int(end) + self.window)
        for item in pysam.faidx(self.reference_genome, range):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()

        print >> hOUT, '>' + label + "_ref"
        print >> hOUT, seq

        # for insertion
        if ref == "-":   seq = seq[0:(self.window + 1)] + alt + seq[-self.window:]
        # for deletion
        elif alt == "-": seq = seq[0:self.window] + seq[-self.window:]
         # for SNV
        else:            seq = seq[0:self.window] + alt + seq[-self.window:]

        print >> hOUT, '>' + label + "_alt"
        print >> hOUT, seq

        hOUT.close()


    ############################################################
    def extractRead(self, bamfile, chr,start,end, output):

        hOUT = open(output, 'w')

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
                print >> hOUT, '>' + read.qname + '/1'
                print >> hOUT, tempSeq
                # print '>' + read.qname + '/1'
                # print tempSeq
            else:
                print >> hOUT, '>' + read.qname + '/2'
                print >> hOUT, tempSeq
                # print '>' + read.qname + '/2'
                # print tempSeq

        hOUT.close()


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
    def checkSecondBestAlignment(self, align):

        # return 1 if there is another alignment whose number of matches if close to the second best alignemt
        if len(align) >= 2:
            align.sort(key=lambda x:x[0],reverse=True)
            if (abs(align[0][0] - align[1][0]) <= range_of_close_best_alignment):
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

        hOUT = open(output, 'w')
        
        seq = ""
        label = ','.join([chr, str(start), str(end), ref, alt])
        range = chr + ":" + str(int(start) - self.window + 1) +"-"+ str(int(end) + self.window)
        for item in pysam.faidx(self.reference_genome, range):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()
            seq = seq.replace('>', '')
            seq = seq.replace(range, '')

        print >> hOUT, '>' + label + "_ref"
        print >> hOUT, seq
        # print '>' + label + "_ref"
        # print seq

        # for insertion
        if ref == "-":   seq = seq[0:(self.window + 1)] + alt + seq[-self.window:]
        # for deletion
        elif alt == "-": seq = seq[0:self.window] + seq[-self.window:]
         # for SNV
        else:            seq = seq[0:self.window] + alt + seq[-self.window:]

        print >> hOUT, '>' + label + "_alt"
        print >> hOUT, seq
        # print '>' + label  + "_alt"
        # print seq

        hOUT.close()

    ############################################################
    def summarizeRefAlt(self, inputPsl):
        
        hIN = open(inputPsl, 'r')

        numOther = []
        numAlt = []
        numRef = []

        tempID = ""
        tempAlt = []
        tempRef = []
        ####
        for line in hIN:
            # print line.rstrip('\n')

            F = line.rstrip('\n').split('\t')
            if F[0].isdigit() == False: continue

            # remove the read pair num info ("/1", "/2") 
            # F[9] = F[9][0:-2]
            if tempID != F[9]:
                if tempID != "":
                
                    ####
                    if (self.checkSecondBestAlignmentOriginal(tempAlt,tempRef) == 0):
                    # if (self.checkSecondBestAlignment(tempAlt) == 0 and self.checkSecondBestAlignment(tempRef) == 0):
                        tempAltScore, tempAltNM = self.getScore(tempAlt)
                        tempRefScore, tempRefNM = self.getScore(tempRef)
                        # print str(tempRefScore) +" " + str(tempAltScore)
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

        hIN.close()

        ####
        if (len(tempAlt) > 0 and len(tempRef) > 0):
            if (self.checkSecondBestAlignmentOriginal(tempAlt,tempRef) == 0):
                # if (self.checkSecondBestAlignment(tempAlt) == 0 and self.checkSecondBestAlignment(tempRef) == 0):
                tempAltScore, tempAltNM = self.getScore(tempAlt)
                tempRefScore, tempRefNM = self.getScore(tempRef)
                # print str(tempRefScore) +" " + str(tempAltScore)
                if tempAltScore == tempRefScore: numOther.append(tempID[0:-2])
                elif tempAltScore >  tempRefScore: numAlt.append(tempID[0:-2])
                elif tempAltScore <  tempRefScore: numRef.append(tempID[0:-2])

        return([len(set(numRef)), len(set(numAlt)), len(set(numOther))])
    

    ############################################################
    def blat_read_count(self, samfile,chr,start,end,output):

        # extract short reads from tumor sequence data around the candidate
        self.extractRead(samfile,chr,start,end,output + ".tmp.fa")
        self.extractRead(samfile,chr,start,end,output + ".tmp.fa")
        # alignment tumor short reads to the reference and alternative sequences
        FNULL = open(os.devnull, 'w')
        subprocess.check_call(self.blat_cmds + [output + ".tmp.refalt.fa", output + ".tmp.fa", output + ".tmp.psl"], 
                              stdout = FNULL, stderr = subprocess.STDOUT)
        FNULL.close()
        # summarize alignment results
        ref, alt, other = self.summarizeRefAlt(output + ".tmp.psl")
        return (ref, alt, other)


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
    def filter(self, in_tumor_bam, in_normal_bam, output, in_mutation_file):

        srcfile = open(in_mutation_file,'r')
        hResult = open(output,'w')

        if in_tumor_bam and in_normal_bam:
            tumor_samfile = pysam.Samfile(in_tumor_bam, "rb")
            normal_samfile = pysam.Samfile(in_normal_bam, "rb")

            if self.header_flag:
                header = srcfile.readline().rstrip('\n')  
                newheader = ("RefNum_tumor\tAltNum_tumor\tOtherNum_tumor"
                         + "\tRefNum_normal\tAltNum_normal\tOtherNum_normal")
                print >> hResult, (header +"\t"+ newheader)

            ####
            for line in srcfile:
                line = line.rstrip()
                itemlist = line.split('\t')
                # annovar input file (not zero-based number)
                chr, start, end, ref, alt  = (itemlist[0], (int(itemlist[1]) - 1), int(itemlist[2]), itemlist[3], itemlist[4])
                
                tumor_ref, tumor_alt, tumor_other, normal_ref, normal_alt, normal_other, log10_fisher_pvalue= ('---','---','---','---','---','---','---')
                self.makeTwoReference(chr,start,end,ref,alt,output + ".tmp.refalt.fa")

                if tumor_samfile.count(chr,start,end) < self.max_depth and int(start) >= int(self.window):
                    tumor_ref, tumor_alt, tumor_other = self.blat_read_count(tumor_samfile, chr, start, end, output)

                if normal_samfile.count(chr,start,end) < self.max_depth and int(start) >= int(self.window):
                    normal_ref, normal_alt, normal_other = self.blat_read_count(normal_samfile, chr, start, end, output)

                if tumor_ref != '---' and  tumor_alt != '---' and  normal_ref != '---' and  normal_alt != '---':
                    log10_fisher_pvalue = self.calc_fisher_pval(tumor_ref, normal_ref, tumor_alt, normal_alt)

                if  ((tumor_alt == '---' or tumor_alt >= self.tumor_min_mismatch) and
                    (normal_alt == '---' or normal_alt <= self.normal_max_mismatch)):
                    print >> hResult, (line +"\t"+ str(tumor_ref)  +"\t"+ str(tumor_alt)  +"\t"+ str(tumor_other)
                                            +"\t"+ str(normal_ref) +"\t"+ str(normal_alt) +"\t"+ str(normal_other)
                                            +"\t"+ str(log10_fisher_pvalue))

            ####
            tumor_samfile.close()
            normal_samfile.close()

        elif in_tumor_bam:
            tumor_samfile = pysam.Samfile(in_tumor_bam, "rb")

            if self.header_flag:
                header = srcfile.readline().rstrip('\n')  
                newheader = ("RefNum_tumor\tAltNum_tumor\tOtherNum_tumor\t0.1\tratio\t0.9")
                print >> hResult, (header +"\t"+ newheader)

            for line in srcfile:
                line = line.rstrip()
                itemlist = line.split('\t')
                # annovar input file (not zero-based number)
                chr, start, end, ref, alt  = (itemlist[0], (int(itemlist[1]) - 1), int(itemlist[2]), itemlist[3], itemlist[4])

                tumor_ref, tumor_alt, tumor_other, beta_01, beta_mid, beta_09 = ('---','---','---','---','---','---')
               
                if tumor_samfile.count(chr,start,end) < self.max_depth and int(start) >= int(self.window) :

                    self.makeTwoReference(chr,start,end,ref,alt,output + ".tmp.refalt.fa")
                    tumor_ref, tumor_alt, tumor_other = self.blat_read_count(tumor_samfile, chr, start, end, output)
                    beta_01, beta_mid, beta_09 = self.calc_btdtri(tumor_ref, tumor_alt)

                if (tumor_alt == '---' or tumor_alt >= self.tumor_min_mismatch):
                    print >> hResult, (line +"\t"+ str(tumor_ref)  +"\t"+ str(tumor_alt)  +"\t"+ str(tumor_other) +"\t"+ str(beta_01) +"\t"+ str(beta_mid) +"\t"+ str(beta_09))
            
            ####
            tumor_samfile.close()

        ####
        hResult.close()
        srcfile.close()

        ####
        if os.path.exists(output + ".tmp.refalt.fa"): os.unlink(output + ".tmp.refalt.fa")
        if os.path.exists(output + ".tmp.fa"): os.unlink(output + ".tmp.fa")
        if os.path.exists(output + ".tmp.psl"): os.unlink(output + ".tmp.psl")


    ############################################################
    def filter_vcf(self, in_tumor_bam, in_normal_bam, output, in_mutation_file, tumor_sample, normal_sample):

        import collections
        import vcf
        import copy

        vcf_reader = vcf.Reader(filename = in_mutation_file)
        f_keys = vcf_reader.formats.keys() #it's an ordered dict
        # add vcf header info
        if in_tumor_bam and in_normal_bam:
            vcf_reader.infos['FPR'] = vcf.parser._Info('FPR', 1, 'Float', "Minus logarithm of the p-value by Fishers exact test processed with Realignment Filter","MutationFilter","v0.2.0")
        elif in_tumor_bam:
            vcf_reader.infos['B1R'] = vcf.parser._Info('B1R', 1, 'Float', "10% posterior quantile of the beta distribution with Realignment Filter","MutationFilter","v0.2.0")
            vcf_reader.infos['BMR'] = vcf.parser._Info('BMR', 1, 'Float', "Posterior mean processed with Realignmnt Filter")
            vcf_reader.infos['B9R'] = vcf.parser._Info('B9R', 1, 'Float', "90% posterior quantile of the beta distribution processed with Realignment Filter","MutationFilter","v0.2.0")

        vcf_reader.formats['NNR'] = vcf.parser._Format('NNR', 1, 'Integer', "Number of non-allelic reads")
        vcf_reader.formats['NAR'] = vcf.parser._Format('NAR', 1, 'Integer', "Number of allelic reads")
        vcf_reader.formats['NOR'] = vcf.parser._Format('NOR', 1, 'Integer', "Number of other reads")
        new_keys = vcf_reader.formats.keys()
        sample_list = vcf_reader.samples

        vcf_writer = vcf.Writer(open(output, 'w'), vcf_reader)

        if in_tumor_bam and in_normal_bam:
            tumor_samfile = pysam.Samfile(in_tumor_bam, "rb")
            normal_samfile = pysam.Samfile(in_normal_bam, "rb")

            ####
            for record in vcf_reader:
                new_record = copy.deepcopy(record)
                chr, start, end, ref, alt, is_conv = vcf_utils.vcf_fields2anno(record.CHROM, record.POS, record.REF, record.ALT[0])
                
                tumor_ref, tumor_alt, tumor_other, normal_ref, normal_alt, normal_other, log10_fisher_pvalue= ('.','.','.','.','.','.','.')
                self.makeTwoReference(chr,start,end,ref,alt,output + ".tmp.refalt.fa")

                if tumor_samfile.count(chr,start,end) < self.max_depth and int(start) >= int(self.window) and is_conv:
                    tumor_ref, tumor_alt, tumor_other = self.blat_read_count(tumor_samfile, chr, start, end, output)
                
                if normal_samfile.count(chr,start,end) < self.max_depth and int(start) >= int(self.window) and is_conv:
                    normal_ref, normal_alt, normal_other = self.blat_read_count(normal_samfile, chr, start, end, output)

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
                    f_vals = [record.samples[sx].data[vx] for vx in range(len(f_keys))]
                    handy_dict = dict(zip(f_keys, f_vals))
                    handy_dict['NNR'] = tumor_ref
                    handy_dict['NAR'] = tumor_alt
                    handy_dict['NOR'] = tumor_other
                    new_vals = [handy_dict[x] for x in new_keys]
                    new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)
                    ## normal sample
                    sx = sample_list.index(normal_sample)
                    new_record.samples[sx].data = collections.namedtuple('CallData', new_keys)
                    f_vals = [record.samples[sx].data[vx] for vx in range(len(f_keys))]
                    handy_dict = dict(zip(f_keys, f_vals))
                    handy_dict['NNR'] = normal_ref
                    handy_dict['NAR'] = normal_alt
                    handy_dict['NOR'] = normal_other
                    new_vals = [handy_dict[x] for x in new_keys]
                    new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)

                    vcf_writer.write_record(new_record)
             
            ####
            tumor_samfile.close()
            normal_samfile.close()

        elif in_tumor_bam:
            tumor_samfile = pysam.Samfile(in_tumor_bam, "rb")

            ####
            for record in vcf_reader:
                # chr, start, end, ref, alt  = (rec.chrom (rec.pos - 1), rec.pos, rec.ref, rec.alts[0])
                chr, start, end, ref, alt, is_conv = vcf_utils.vcf_fields2anno(record.CHROM, record.POS, record.REF, record.ALT[0])
                
                tumor_ref, tumor_alt, tumor_other, beta_01, beta_mid, beta_09 = ('','','','','','')
               
                if tumor_samfile.count(chr,start,end) < self.max_depth and int(start) >= int(self.window) :
                    self.makeTwoReference(chr,start,end,ref,alt,output + ".tmp.refalt.fa")
                    tumor_ref, tumor_alt, tumor_other = self.blat_read_count(tumor_samfile, chr, start, end, output)
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
                    new_record.samples[sx_sample].data = collections.namedtuple('CallData', new_keys)
                    f_vals = [record.samples[sx_sample].data[vx] for vx in range(len(f_keys))]
                    handy_dict = dict(zip(f_keys, f_vals))
                    handy_dict['NNR'] = tumor_ref
                    handy_dict['NAR'] = tumor_alt
                    handy_dict['NOR'] = tumor_other
                    new_vals = [handy_dict[x] for x in new_keys]
                    new_record.samples[sx].data = new_record.samples[sx].data._make(new_vals)

                    vcf_writer.write_record(new_record)
             
            ####
            tumor_samfile.close()

        ####
        vcf_writer.close()

        ####
        if os.path.exists(output + ".tmp.refalt.fa"): os.unlink(output + ".tmp.refalt.fa")
        if os.path.exists(output + ".tmp.fa"): os.unlink(output + ".tmp.fa")
        if os.path.exists(output + ".tmp.psl"): os.unlink(output + ".tmp.psl")

