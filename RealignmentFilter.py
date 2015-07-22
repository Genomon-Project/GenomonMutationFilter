import sys
import os
import re
import pysam
import scipy.special
import argparse
import logging
import subprocess
import math
from AutoVivification import AutoVivification

#
# Class definitions
#
class RealignmentFilter:

    def __init__(self, mutationFileInfo):
        self.mfi = mutationFileInfo
        self.max_depth = int(10000)
        self.max_mismatch = int(1000)

        self.reference_genome = '/home/ogawaprj/ngs/ref/GRCh37-lite_PCAWG_bwa-0.7.10/GRCh37-lite_PCAWG.fa'
        self.window = int(200)

        self.tumor_min_mismatch = int(3)
        self.normal_max_mismatch = int(1)

        self.complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        # self.blat_cmds = ['/usr/local/bin/blat', '-stepSize=5', '-repMatch=2253']
        self.blat_cmds = ['/home/ogawaprj/ngs/bin/blat_x86_64/blat', '-fine']
        
        
    ############################################################
    def set_filter_info(self, mutation_line_info, tumor_ref,tumor_alt,tumor_other, normal_ref, normal_alt, normal_other):
            
        if(tumor_alt >= self.tumor_min_mismatch and normal_alt <= self.normal_max_mismatch):
            mutation_line_info[2].append(0)
        else:
            mutation_line_info[2].append(1)
           
        mutation_line_info[3].append(str(tumor_ref)  +"\t"+ str(tumor_alt)  +"\t"+ str(tumor_other) +"\t" +
                                     str(normal_ref) +"\t"+ str(normal_alt) +"\t"+ str(normal_other))


    ############################################################
    def makeTwoReference(self, chr,start,end,ref,alt, outputFilePath):

        hOUT = open(outputFilePath, 'w')
        
        seq = ""
        label = ','.join([chr, str(start), str(end), ref, alt])
        range = chr + ":" + str(int(start) - self.window + 1) +"-"+ str(int(end) + self.window)
        for item in pysam.faidx(self.reference_genome, range):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()

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
    def extractRead(self, bamfile, chr,start,end, outputFilePath):

        hOUT = open(outputFilePath, 'w')

        for read in bamfile.fetch(chr,start,end):
          
            # get the flag information
            flags = format(int(read.flag), "#014b")[:1:-1]

            # skip improper pair
#            if flags[1] == "0": continue 

            # skip unmapped read 
#            if flags[2] == "1" or flags[3] == "1": continue 

            # skip supplementary alignment
            if flags[8] == "1" or flags[11] == "1": continue

            # skip duplicated reads
            if flags[10] == "1": continue

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
    def checkSecondBestAlignment(self, align1, align2):

        align1tmp = list(align1)
        align1tmp.extend(align2)

        if len(align1tmp) == 1:
              return 1
        # return 1 if there is another alignment whose number of matches if close to the second best alignemt
        align1tmp.sort(key=lambda x:x[0],reverse=True)
        # if len(align1tmp) >= 3 and (abs(align1tmp[1][1] - align1tmp[2][1]) <= 5):
        if len(align1tmp) == 3 and (align1tmp[1][0] - align1tmp[2][0]) <= 5:
            return 1
        return 0


    ############################################################
    def checkSecondBestAlignmentxxxxx(self, align):

        # return 1 if there is another alignment whose number of matches if close to the second best alignemt
        if len(align) >= 2:
            align.sort(key=lambda x:x[0],reverse=True)
            if (abs(align[0][0] - align[1][0]) <= 5):
                return 1
        return 0


    ############################################################
    def getScore(self, align):
        if len(align) >= 1:
            align.sort(key=lambda x:x[0],reverse=True)
            return (align[0][0], align[0][1])
        return (0,0)

    ############################################################
    def makeTwoReference(self, chr,start,end,ref,alt, outputFilePath):

        hOUT = open(outputFilePath, 'w')
        
        seq = ""
        label = ','.join([chr, str(start), str(end), ref, alt])
        range = chr + ":" + str(int(start) - self.window + 1) +"-"+ str(int(end) + self.window)
        for item in pysam.faidx(self.reference_genome, range):
            if item[0] == ">": continue
            seq = seq + item.rstrip('\n').upper()

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
            F[9] = F[9][0:-2]
            if tempID != F[9]:
                if tempID != "":
                
                    ####
                    if (self.checkSecondBestAlignment(tempAlt,tempRef) == 0):
                    # if ((self.checkSecondBestAlignment(tempAlt) == 0) and (self.checkSecondBestAlignment(tempRef) == 0)):
                        tempAltScore, tempAltNM = self.getScore(tempAlt)
                        tempRefScore, tempRefNM = self.getScore(tempRef)
                        # if tempAltNM >= self.max_mismatch and tempRefNM >= self.max_mismatch: numOther.append(tempID)
                        if tempAltScore == tempRefScore: numOther.append(tempID)
                        elif tempAltScore >  tempRefScore: numAlt.append(tempID)
                        elif tempAltScore <  tempRefScore: numRef.append(tempID)

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
        if (self.checkSecondBestAlignment(tempAlt,tempRef) == 0):
        # if ((self.checkSecondBestAlignment(tempAlt) == 0) and (self.checkSecondBestAlignment(tempRef) == 0)):
            tempAltScore, tempAltNM = self.getScore(tempAlt)
            tempRefScore, tempRefNM = self.getScore(tempRef)
            # if tempAltNM >= self.max_mismatch and tempRefNM >= self.max_mismatch: numOther.append(tempID)
            if tempAltScore == tempRefScore: numOther.append(tempID)
            elif tempAltScore >  tempRefScore: numAlt.append(tempID)
            elif tempAltScore <  tempRefScore: numRef.append(tempID)


        return([len(set(numRef)), len(set(numAlt)), len(set(numOther))])
               

    ############################################################
    def filter(self, in_tumor_bam, in_normal_bam, outputFilePath):
        logging.info( 'filter start')

        tumor_samfile = pysam.Samfile(in_tumor_bam, "rb")
        normal_samfile = pysam.Samfile(in_normal_bam, "rb")
        
        for mutation_line_info in self.mfi['line_info']:
            chr = mutation_line_info[1]['chr']
            start = mutation_line_info[1]['start']
            end = mutation_line_info[1]['end']
            ref = mutation_line_info[1]['ref']
            alt = mutation_line_info[1]['alt']

            if tumor_samfile.count(chr, start, (start+1)) >= self.max_depth:
                self.set_filter_info(mutation_line_info,0,0,0,0,0,0)
                continue 

            if normal_samfile.count(chr, start, (start+1)) >= self.max_depth:
                self.set_filter_info(mutation_line_info,0,0,0,0,0,0)
                continue 


            self.makeTwoReference(chr,start,end,ref,alt,outputFilePath + "/tmp.refalt.fa")
            
            # extract short reads from tumor sequence data around the candidate
            self.extractRead(tumor_samfile,chr,start,end,outputFilePath + "/tmp.fa")

            # alignment tumor short reads to the reference and alternative sequences
            FNULL = open(os.devnull, 'w')
            subprocess.call(self.blat_cmds + [outputFilePath + "/tmp.refalt.fa", outputFilePath + "/tmp.fa", outputFilePath + "/tmp.psl"], 
                            stdout = FNULL, stderr = subprocess.STDOUT)
            FNULL.close()

            # summarize alignment results
            tumor_ref, tumor_alt, tumor_other = self.summarizeRefAlt(outputFilePath + "/tmp.psl")

            # extract short reads from normal sequence data around the candidate
            self.extractRead(normal_samfile,chr,start,end,outputFilePath + "/tmp.fa")

            # alignment normal short reads to the reference and alternative sequences
            FNULL = open(os.devnull, 'w')
            subprocess.call(self.blat_cmds + [outputFilePath + "/tmp.refalt.fa", outputFilePath + "/tmp.fa", outputFilePath + "/tmp.psl"], 
                            stdout = FNULL, stderr = subprocess.STDOUT)
            FNULL.close()

            # summarize alignment results
            normal_ref, normal_alt, normal_other = self.summarizeRefAlt(outputFilePath + "/tmp.psl")

            self.set_filter_info(mutation_line_info, tumor_ref,tumor_alt,tumor_other, normal_ref, normal_alt, normal_other)

    logging.info( 'filter end')

