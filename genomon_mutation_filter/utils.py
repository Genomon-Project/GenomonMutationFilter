from __future__ import print_function
import os, sys
import vcf
import math

def vcf_fields2anno(chrom, pos_str, ref_sub, alt_sub):
    pos = int(pos_str)
    ref = str(ref_sub)
    alt = str(alt_sub)

    # for insertion
    if len(ref) < len(alt) and len(ref) == 1 and alt[0:1] == ref:
        start = pos -1
        end = pos
        ret = (chrom, start, end, "-", alt[1:], True)

    # for deletion
    elif len(ref) > len(alt) and len(alt) == 1 and ref[0:1] == alt:
        start = pos
        end = pos + len(ref[1:])
        ret = (chrom, start, end, ref[1:], "-", True)

    # for SNV
    elif len(ref) == 1 and len(alt) == 1:
        start = pos -1
        end = pos
        ret = (chrom, start, end, ref, alt, True)

    # for block substitution
    else:
        start = pos - 1
        end = pos
        ret = (chrom, start, end, ref, alt, False)

    return ret


def sort_header(in_header, out_header):

    fileformat_list = []
    filter_list = []
    info_list = []
    format_list = []
    contig_list = []
    reference_list = []
    header_list = []
    others_list = []

    with open(in_header,'r') as hin:
        for line in hin:
            line = line.rstrip()
            if line.startswith("##fileformat="):
                fileformat_list.append(line)
            elif line.startswith("##FILTER="):
                filter_list.append(line)
            elif line.startswith("##INFO="):
                info_list.append(line)
            elif line.startswith("##FORMAT="):
                format_list.append(line)
            elif line.startswith("##contig="):
                contig_list.append(line)
            elif line.startswith("##reference="):
                reference_list.append(line)
            elif line.startswith("#CHROM"):
                header_list.append(line)
            else:
                others_list.append(line)

    out_list = fileformat_list \
             + filter_list \
             + info_list \
             + format_list \
             + contig_list \
             + reference_list \
             + others_list \
             + header_list 

    with open(out_header, "w") as hout:
        for val in out_list:
            print(val, file=hout)


############################################################
def partition_anno(inputFilePath, thread_num):

    partitionNum = thread_num
    outputFilePrefix = inputFilePath 

    # count the number of lines
    with open(inputFilePath, "r") as hin:
        recordNum = sum(1 for line in hin)

    partitionNum_mod = min(recordNum, partitionNum)
    if partitionNum_mod == 0: partitionNum_mod = 1
    eachPartitionNum = math.ceil(recordNum / partitionNum_mod)
    if eachPartitionNum == 0: eachPartitionNum = 1
    partitionNum_mod = math.ceil(recordNum / eachPartitionNum)

    currentPartition = 0
    currentRecordNum = 0
    with open(inputFilePath) as hin:
        hout = open(outputFilePrefix + ".1", 'w')
        for line in hin:
            print(line.rstrip('\n'), file=hout)
            currentRecordNum += 1

            if currentRecordNum >= eachPartitionNum and currentPartition < partitionNum_mod - 1:
                currentRecordNum = 0
                currentPartition += 1
                hout.close()
                hout = open(outputFilePrefix +"."+ str(currentPartition+1), 'w')

        hout.close()
    return partitionNum_mod


############################################################
def partition_vcf(inputFilePath, thread_num):
    partitionNum = thread_num
    outputFilePrefix = inputFilePath 

    hin = open(inputFilePath, 'r')
    vcf_reader1 = vcf.Reader(hin)
    recordNum = 0
    for record in vcf_reader1:
        recordNum += 1
    hin.close()

    partitionNum_mod = min(recordNum, partitionNum)
    if partitionNum_mod == 0: partitionNum_mod = 1
    eachPartitionNum = math.ceil(recordNum / partitionNum_mod)
    if eachPartitionNum == 0: eachPartitionNum = 1
    partitionNum_mod = math.ceil(recordNum / eachPartitionNum)
    
    currentPartition = 0
    currentRecordNum = 0
    
    hin = open(inputFilePath, 'r')
    vcf_reader2 = vcf.Reader(hin)

    hout = open(outputFilePrefix + ".1", 'w')
    vcf_writer = vcf.Writer(hout, vcf_reader2)

    for record in vcf_reader2:
        vcf_writer.write_record(record)
        currentRecordNum += 1

        if currentRecordNum >= eachPartitionNum and currentPartition < partitionNum_mod - 1:
            currentPartition += 1
            currentRecordNum = 0
            vcf_writer.close()
            hout.close()
            hout = open(outputFilePrefix +"."+ str(currentPartition+1), 'w')
            vcf_writer = vcf.Writer(hout, vcf_reader2) 
    vcf_writer.close()
    hout.close()
    hin.close()
    return partitionNum_mod

