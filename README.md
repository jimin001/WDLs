# WDLs

## Illumina reads alignment
bwa_mem2.wdl

## ONT reads alignment
minimap2_ubam.wdl 
  - starting with unaligned bams and saving methylation tags

minimap2_fastq.wdl
  - starting with fastqs

> additional args:
> ```
> -k 17 -y -K 5G --eqx
> ```

## HiFi reads alignment
pbmm2.wdl
