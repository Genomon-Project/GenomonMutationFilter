# GenomonMutationFilter

## Install

```
git clone https://github.com/Genomon-Project/GenomonMutationFilter.git
cd GenomonMutationFilter
python setup.py build
python setup.py install
```

```
wget  http://downloads.sourceforge.net/project/samtools/tabix/tabix-0.2.6.tar.bz2
tar xjvf tabix-0.2.6.tar.bz2
cd tabix-0.2.6
make
```

## Run

```
mutfilter realignment [-h] [-t tumor_min_mismatch]
                           [-n normal_max_mismatch] [-s score_difference]
                           [-w window_size]
                           target.vcf target_tumor.bam target_normal.bam
                           output.vcf ref_genome blat_path
```

```
mutfilter indel [-h] [-s search_length] [-n neighbor]
                     [-b base_qual_thres] [-d min_depth] [-m min_mismatch]
                     [-a allele_frequency_thres]
                     target.vcf target.bam output.vcf
```

```
mutfilter breakpoint [-h] [-d max_depth] [-c min_clip_size]
                          [-j junc_num_thres] [-m mapping_quality_thres]
                          target.vcf target.bam output.vcf
```

```
usage: mutfilter simplerepeat [-h]
                              target.vcf output.vcf simple_repeat_database
```
