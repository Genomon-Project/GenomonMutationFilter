from __future__ import print_function
import os, sys
import vcf

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
def partition_anno(in_mutation_file, thread_num):

    # count the number of lines
    with open(in_mutation_file, "r") as hin:
        line_num = sum(1 for line in hin)

    thread_num_mod = min(line_num, thread_num)
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
def partition_vcf(in_mutation_file, thread_num):

    # count the number of lines
    line_num = 0
    hin = open(in_mutation_file, 'r')
    vcf_reader = vcf.Reader(hin)
    for record in vcf_reader:
        line_num = line_num + 1
    hin.close()
    
    thread_num_mod = min(line_num, thread_num)
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
    
    