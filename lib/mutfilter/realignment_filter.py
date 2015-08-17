import sys
import os
import re
import pysam
import logging
import subprocess

#
# Class definitions
#
class realignment_filter:

    def __init__(self,referenceGenome,tumor_min_mismatch,normal_max_mismatch, search_length, score_difference, blat, header_flag):
        self.reference_genome = referenceGenome
        self.window = search_length
        self.score_difference = score_difference
        self.tumor_min_mismatch = tumor_min_mismatch
        self.normal_max_mismatch = normal_max_mismatch
        self.blat_cmds = [blat, '-fine']
        self.complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
        self.header_flag = header_flag
     
        
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
    def filter(self, in_tumor_bam, in_normal_bam, output, in_mutation_file):

        tumor_samfile = pysam.Samfile(in_tumor_bam, "rb")
        normal_samfile = pysam.Samfile(in_normal_bam, "rb")

        srcfile = open(in_mutation_file,'r')
        header = ""
        if self.header_flag:
            header = srcfile.readline().rstrip('\n')  
        
        else: # no header line
            line = srcfile.readline().rstrip()
            column_len = len(line.split('\t'))
            srcfile.close()
            for num in range (2, column_len):
                header = header + "\t"
            srcfile = open(in_mutation_file,'r')
        
        hResult = open(output,'w')

        # print header
        newheader = ("readPairNum_tumor\tvariantPairNum_tumor\totherPairNum_tumor"
                + "\treadPairNum_normal\tvariantPairNum_normal\totherPairNum_normal")
        print >> hResult, (header +"\t"+ newheader)

        ####
        for line in srcfile:
            line = line.rstrip()
            itemlist = line.split('\t')
    
            # annovar input file (not zero-based number)
            chr = itemlist[0]
            start = (int(itemlist[1]) - 1)
            end = int(itemlist[2])
            ref = itemlist[3]
            alt = itemlist[4]
            
            self.makeTwoReference(chr,start,end,ref,alt,output + ".tmp.refalt.fa")
            
            # extract short reads from tumor sequence data around the candidate
            self.extractRead(tumor_samfile,chr,start,end,output + ".tmp.fa")

            # alignment tumor short reads to the reference and alternative sequences
            FNULL = open(os.devnull, 'w')
            retcode = subprocess.check_call(self.blat_cmds + [output + ".tmp.refalt.fa", output + ".tmp.fa", output + ".tmp.psl"], 
                                            stdout = FNULL, stderr = subprocess.STDOUT)
            FNULL.close()

            # summarize alignment results
            tumor_ref, tumor_alt, tumor_other = self.summarizeRefAlt(output + ".tmp.psl")
            
            # extract short reads from normal sequence data around the candidate
            self.extractRead(normal_samfile,chr,start,end,output + ".tmp.fa")

            # alignment normal short reads to the reference and alternative sequences
            FNULL = open(os.devnull, 'w')
            subprocess.check_call(self.blat_cmds + [output + ".tmp.refalt.fa", output + ".tmp.fa", output + ".tmp.psl"], 
                                  stdout = FNULL, stderr = subprocess.STDOUT)
            FNULL.close()

            # summarize alignment results
            normal_ref, normal_alt, normal_other = self.summarizeRefAlt(output + ".tmp.psl")

            if(tumor_alt >= self.tumor_min_mismatch and normal_alt <= self.normal_max_mismatch):
                print >> hResult, (line +"\t"+ str(tumor_ref)  +"\t"+ str(tumor_alt)  +"\t"+ str(tumor_other)
                                        +"\t"+ str(normal_ref) +"\t"+ str(normal_alt) +"\t"+ str(normal_other))

        ####
        hResult.close()
        srcfile.close()


