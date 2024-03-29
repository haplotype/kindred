# Kindred 

Kindred is designed to infer kinship $\phi$ between a pair of samples accounting for inbreeding of the samples. This is done through modelling nine latent Jacquard IBD states. If the pair is one and itself, the kinship is $\phi = (1+F)/2$ where $F$ denotes inbreeding coefficient. 

## Current version

Verion 0.81 was released on 13 March 2023. Exectuables for Linux and Mac can be found here: http://www.haplotype.org/software.html. 
You may also choose to compile from the source code in src/.  The manuscript is here: https://doi.org/10.1093/bioinformatics/btae082. 

## Input and options
Kindred takes vcf file or stream as input with mandatory -i option. The vcf file must contain "INFO/AF" tag and biallelic genotypes. If the AF tag is not readily available, it can be populated by bcftools on the fly. 
User can also specify other tag name, such as EUR_AF or EAS_AF in the 1000 genomes vcf files, via -a option. 
Kindred uses multi-threading to speed up calculation and user can specify the number of threads with -t option. 

## Output
Kindred by default outputs two files: pref.log and pref.rkm.gz (rkm stands for realized kinship matrix),  where pref can be specified with -o option. 
In pref.rkm.gz is an $n\times n$ square matrix of $(2\phi_{ij})$ with 6 digits accuracy. The samples have the same order as in the vcf file.   

    1.000276 0.000000 0.480721 0.473312 0.000000 0.000000 ...
    0.000000 0.998443 0.468201 0.478238 0.000000 0.000000 ...
    0.480721 0.468201 0.999667 0.502964 0.000000 0.000000 ...
    0.473312 0.478238 0.502964 0.999840 0.000000 0.000000 ...
    0.000000 0.000000 0.000000 0.000000 0.997251 0.000000 ...
    0.000000 0.000000 0.000000 0.000000 0.000000 0.999693 ...
    ...

       
With -k option, Kindred can output an additional file pref.kin.gz.  The first line is a space delimited header, and rest are the pairwise kinship (including one and itself) and the detailed parameter estimates, total $(n+1)n/2$ lines. 

    ID1 ID2 phi sumd d1 d2 d3 d4 d5 d6 d7 d8 d9
    pg1 pg1 0.50014 1.00051 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 1.00028 0.00000 0.00023 
    pg1 pg2 0.00000 0.99758 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.99758 
    pg1 pg3 0.24036 1.00040 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.02403 0.91337 0.06299
    pg1 pg4 0.23666 1.00009 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.94662 0.05347
    pg1 pg5 0.00000 0.99922 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000 0.99922
    ...

ID1 and ID2 are two sample IDs, phi is the kinship, d1, ..., d9 are probabilities of Jacquard latent states, and sumd is the sum of the nine probabilities, should be close to $1$. Significant departure from $1$ suggest ill fit, usually caused by misspecification of allele frequencies.  
  
## UKin implementation

To compare Kindred with UKin (https://github.com/zxy22320/UKin), we reimplemented it with multithreading in c and this can be invoked by -u. It will take the same input vcf file as Kindred and output a kinship estimates (a square matrix as in pref.rkm.gz) in a file named pref.ukinGRM.gz.   

## Usage examples

1) If in.vcf.gz contains AF tag, one can simply do: 
    
       kindred -i in.vcf.gz -o pref 

2) If the pairwise kinship computation appears slow and you have many idle cores: 

       kindred -i in.vcf.gz -o pref -t 50

3) If in.vcf.gz has no AF tag, and one wants to use allele frequencies estimated from genotypes in the vcf file: 

       bcftools plugin fill-tags in.vcf.gz | kindred -i - -o pref 

    or equivalently:
    
       bcftools +fill-tags in.vcf.gz | kindred -i - -o pref 

4) If in.vcf.gz contains a AF tag, command in 2) will replace AF values, but if you insist you can remove the original tag first:  

       bcftools annotate --remove INFO in.vcf.gz | bcftools +fill-tags | \
          kindred -i - -o pref 

5) You may store precomputed allele frequencies in a vcf file "annotate.vcf.gz". Kindred can use the allele frequencies by  

       bcftools annotate --remove INFO -c 'INFO/AF' -a annotate.vcf.gz in.vcf.gz  | \
          kindred -i - -o pref 
          
   It's safe to use "--remove" before annotation as you might only have a subset of SNPs annotated and you want ignore the rest. 

6) You may store multiple allele frequencies in an annoation file, and you want use EUR_AF instead of EAS_AF or AFR_AF: 
  
       bcftools annotate -c 'INFO/EUR_AF' -a annotate.vcf.gz in.vcf.gz  | \
          kindred -i - -o pref -a EUR_AF

7) You may use markers on chromosome 8 (or a region) to compute kinship with sample estimated allele frquencies:  

       bcftools filter -r 8 in.vcf.gz | bcftools +fill-tags  | \
          kindred -i - -o test.chr8

8) You may do it with allele frequencies stored in an annotation file:   

       bcftools annotate -c 'INFO/AF' -a annotate.vcf.gz in.vcf.gz filter -r 8 | \
          kindred -i - -o test.chr8 


## Protocol to prepare annotation vcf files
By selecting a subset of SNPs that have similar allele frequenices across population, kindred can infer kinship for admixed samples. Below is a protocol that generates a FST tag that can be used to select such SNPs. 

To prepare annotation files that contains allele frequencies for a specific population,  
we downloaded 1000 genomes vcf and tbi files, they were arranged in different chromosomes. 
We first obtain all allelic SNPs with minimum of 50 counts of minor alleles (among 2504 samples).  
for each SNP count AN-AC for different populations (the example below is for chr22 of CHB). 

     bcftools view -m2 -M2 -v snps -c 50:minor ALL.chr22.phase3.genotypes.vcf.gz | \
        bcftools view -S samples.CHB.txt | bcftools annotate --remove INFO |\
        bcftools plugin fill-AN-AC | \
        bcftools query -f "%CHROM %POS %REF %ALT %INFO/AC %INFO/AN\n" > chb.chr22.an-ac

 
Suppose we obtain for each chrosome the an-ac file for populations CEU and YRI, in addition to CHB.  
For each SNP we  do a chisq test to compute p-values (pv), and write a vcf file with INFO/FST = -10log10(pv). 
We can then use FST to filter SNPs that to be annotated. 

To prepare a vcf file with CEU allele frequencies.  

     bcftools view -S ceu.samples.list annotate --remove INFO +fill-tags | \
        bcftools view -G > annotation.ceu 

## Some tips 
1) plink can convert plink format to vcf format

       plink --bfile prefix --recode vcf --out in

2) Bcftools requires vcf files to be zipped and indexed

       bgzip in.vcf 
       tabix in.vcf.gz 

You will have in.vcf.gz (instead of in.vcf) and in.vcf.gz.tbi. 

3) print all INFO tags in in.vcf.gz 
       
       bcftools query -f '%INFO\n' in.vcf.gz | head -n 2 
       
4) The order of the sample can be obtained from vcf files by  

       bcftools query -l in.vfc.gz
       
5) Rename chromosome name in vcf files.  

       touch chr_name_conv.txt; for i in $(seq 1 22); do echo chr$i $i >> chr_name_conv.txt; done
       bcftools annotate --rename-chrs chr_name_conv.txt input.chr22.vcf.gz | bgzip > output.22.vcf.gz
       
6) Set SNP ID. 
 
       bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT'  in.vcf.gz
       
