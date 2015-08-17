# GenomonMutationFilter

GenomonMutationFilter is a software package for filtering false poistive somatic mutations from cancer genome sequencing data.

## Dependency

```
Python (>= 2.7), pysam, pytabix packages
```

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
                           [-w window_size] [--header]
                           target.anno target_tumor.bam target_normal.bam
                           output.anno ref_genome blat_path
```

```
mutfilter indel [-h] [-s search_length] [-n neighbor]
                     [-b base_qual_thres] [-d min_depth] [-m min_mismatch]
                     [-a allele_frequency_thres] [--header]
                     target.anno target_normal.bam output.anno
```

```
mutfilter breakpoint [-h] [-d max_depth] [-c min_clip_size]
                          [-j junc_num_thres] [-m mapping_quality_thres]
                          [--header]
                          target.anno target_normal.bam output.anno
```

```
mutfilter simplerepeat [-h] [--header]
                            target.anno output.anno simple_repeat_database
```
