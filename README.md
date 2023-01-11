# Kindred 

Kindred is designed to infer kinship $\phi$ between a pair of samples accounting for inbreeding of the samples. This is done through modelling nine latent Jacquard IBD states. If the pair is one and itself, the kinship is $\phi = (1+F)/2$ where $F$ denotes inbreeding coefficient. 

## Input and options
Kindred takes vcf file or stream as input. It is designed to work with bcftools.  The vcf file must contain "INFO/AF" field and biallelic genotypes. If the AF field is not readily available, it can be populated by bcftools on the fly. 
User can also specify other name of the field, such as EUR_AF or EAS_AF in the 1000 genomes vcf files, via -a option. 
Kindred uses multi-threading to speed up calculation and user can specify the number of threads with -t option. 

## Output
Kindred output two files. One is "pref.grm" the other is "pref.kin", where pref can be specified with -o option. 
In pref.grm, the first line contains individual ID. The rest is an $n\times n$ square matrix with 6 digits accuracy. The $i$-th row and $j$-th column is $2\phi_{ij}$. 

    p1_g1-b1-s1 p1_g1-b1-i1 p1_g2-b1-i1 p1_g2-b2-i1 p2_g1-b1-s1 p2_g1-b1-i1 
    1.016231 0.000000 0.506842 0.435748 0.000000 0.000000 ...
    0.000000 0.988070 0.409146 0.461766 0.000000 0.000000
    0.506842 0.409146 1.015925 0.260444 0.000000 0.000000
    0.435748 0.461766 0.260444 0.989434 0.000000 0.000000
    0.000000 0.000000 0.000000 0.000000 1.002591 0.000000
    0.000000 0.000000 0.000000 0.000000 0.000000 1.011074 ...
    ...

In pref.kin, the first line is the header with space delimit, with $(n+1)n/2$ additional lines. 

    ID1 ID2 phi sumd d1 d2 d3 d4 d5 d6 d7 d8 d9
    p1_g1-b1-s1 p1_g1-b1-s1 0.501266 1.000727 0.00213 0.00033 0.00004 0.00000 0.00004 0.00000 0.99816 0.00000 0.00000 
    p1_g1-b1-s1 p1_g1-b1-i1 0.000000 0.999537 0.00000 0.00000 0.00000 0.00013 0.00000 0.00000 0.00000 0.00000 0.99940 
    p1_g1-b1-s1 p1_g2-b1-i1 0.250565 1.000210 0.00000 0.00000 0.00000 0.00329 0.00000 0.00000 0.01257 0.97710 0.00722
    p1_g1-b1-s1 p1_g2-b2-i1 0.245546 1.000082 0.00000 0.00000 0.00442 0.00000 0.00368 0.00000 0.00000 0.96595 0.02601
    p1_g1-b1-s1 p2_g1-b1-s1 0.002957 1.000164 0.00000 0.00000 0.00000 0.00059 0.00000 0.00000 0.00591 0.00000 0.99365
    ...

ID1 and ID2 are two sample IDs, phi is the kinship, d1, ..., d9 are probabilities of Jacquard latent states, and sumd is the sum of the nine probabilities, should be close to $1$.  
  

## Usage examples

1) If in.vcf.gz contains AF field, one can simply do: 
    
       $ kindred in.vcf.gz -o pref 

2) If in.vcf.gz has no AF field, and one wants to use allele frequencies estimated from genotypes in the vcf file. 

       $ bcftools plugin fill-tags in.vcf.gz | kindred - -o test 

3) If in.vcf.gz contain AF field, but you want to recomputed it. 

       $ bcftools annotate --remove INFO in.vcf.gz | \
          bcftools plugin fill-tags in.vcf.gz | \
          kindred - -o test 

4) You may store precomputed allele frequencies in a vcf file "annotate.vcf.gz". Kindred can use the allele frequencies by  

       $ bcftools annotate --remove INFO -c 'INFO/AF' -a annotate.vcf.gz in.vcf.gz  | \
          kindred - -o test 

5) You may store multiple allele frequencies in an annoation file, and you want use EUR_AF instead of EAS_AF or AFR_AF: 
  
       $ bcftools annotate --remove INFO -c 'INFO/EUR_AF' -a annotate.vcf.gz in.vcf.gz  | \
          kindred - -o test -a EUR_AF

6) If you want to only use markers on chromosome 8 to compute kinship using sample estimated allele frquencies:  

       $ bcftools filter -r 8 in.vcf.gz | bcftools plugin fill-tags  | \
          kindred - -o test.chr8

7) If you want to do it with allele frequencies stored in an annotation file:   

       $ bcftools fitler -r 8 in.vcf.gz | \
          bcftools annotate --remove INFO -c 'INFO/AF' -a annotate.vcf.gz  | \
          kindred - -o test.chr8 


## Protocol to prepare annotation vcf files
To prepare annotation files that contains allele frequencies for a specific population,  
we downloaded 1000 genomes vcf and tbi files, they were arranged in different chromosomes. 
We first obtain all allelic SNPs with minimum of 50 counts of minor alleles (among 2504 samples).  
for each SNP count AN-AC for different populations (the example below is for chr22 of CHB). 

     $ bcftools view -m2 -M2 -v snps -c 50:minor ALL.chr22.phase3.genotypes.vcf.gz | \
        bcftools view -S samples.CHB.txt | bcftools annotate --remove INFO |\
        bcftools plugin fill-AN-AC | \
        bcftools query -f "%CHROM %POS %REF %ALT %INFO/AC %INFO/AN\n" > chb.chr22.an-ac

 
Suppose we obtain for each chrosome the an-ac file for populations CEU and YRI, in addition to CHB.  
For each SNP we  do a chisq test to compute p-values (pv), and write a vcf file with INFO/FST = -10log10(pv). 
We can then use FST to filter SNPs that to be annotated. 

To prepare a vcf file with CEU allele frequencies.  

     $ bcftools view -S ceu.samples.list annotate --remove INFO plugin fill-tags | \
        bcftools view -G > annotation.ceu 

## Other useful tips. 
1) plink can convert plink format to vcf format

       $ plink --bfile prefix --recode vcf --out in

2) Bcftools requires vcf files to be zipped and indexed

       $ bgzip in.vcf 
       $ tabix in.vcf.gz 

You will have in.vcf.gz (instead of in.vcf) and in.vcf.gz.tbi. 
