from __future__ import print_function
import os, sys

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

