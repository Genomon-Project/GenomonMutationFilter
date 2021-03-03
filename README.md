# GenomonMutationFilter

GenomonMutationFilter is a software package for filtering false poistive somatic mutations from cancer genome sequencing data.

## Dependency

Python (>= 2.7, >= 3.7), pysam packages

You can choose between realignment using blat or realignment using edlib.

* [blat](http://genome.ucsc.edu/)
* [edlib](https://pypi.org/project/edlib/)

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

## Preparation

  **target somatic mutation candidats**: the somatic mutation candidates (should be .tsv or .vcf format).  
  **target tumor bam**: the indexed bam file of the target tumor sample.  
  **target normal bam**: the indexed bam file of the target normal sample.  


## Run

```
usage: mutfilter realignment [-h] -t TARGET_MUTATION_FILE -1 BAM1 [-2 BAM2]
                             [-A SAMPLE1] [-B SAMPLE2] -o OUTPUT -r REF_GENOME
                             [--blat] [-b BLAT_PATH] [-m tumor_min_mismatch]
                             [-M normal_max_mismatch] [-s score_difference]
                             [-w window_size] [-d max_depth]
                             [-F exclude_sam_flags] [-O {vcf,anno}] [--header]
                             [-T number_of_threads]
```

```
usage: mutfilter indel [-h] -t TARGET_MUTATION_FILE -2 BAM2 [-A SAMPLE1]
                       [-B SAMPLE2] -o OUTPUT [-l search_length] [-n neighbor]
                       [-d min_depth] [-m min_mismatch]
                       [-a allele_frequency_thres] [--header] -s SAMTOOLS_PATH
                       [-S SAMTOOLS_PARAMS] [-O {vcf,anno}] [-r REF_GENOME]
```

```
usage: mutfilter breakpoint [-h] -t TARGET_MUTATION_FILE -2 BAM2 [-A SAMPLE1]
                            [-B SAMPLE2] -o OUTPUT [-d max_depth]
                            [-c min_clip_size] [-j junc_num_thres]
                            [-m mapping_quality_thres] [-F exclude_sam_flags]
                            [--header] [-O {vcf,anno}] [-r REF_GENOME]
```

```
usage: mutfilter simplerepeat [-h] -t TARGET_MUTATION_FILE -o OUTPUT -S
                              SIMPLE_REPEAT_DB [--header] [-O {vcf,anno}]
```
