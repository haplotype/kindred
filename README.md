# kindred
(Infer realized kinship via latent identity-by-descent states)

Kindred is designed to infer kinship $\phi$ between a pair of samples. If the pair is one and itself, the kinship is $\phi = (1+F)/2$ where $F$ denotes inbreeding coefficient. 

## Input and options
Kindred takes vcf file or stream as input. It is designed to work with bcftools.  The vcf file must contain "INFO/AF" field and biallelic genotypes. If the AF field is not readily available, it can be populated by bcftools on the fly. 
User can also specify other name of the field, such as EUR_AF or EAS_AF in the 1000 genomes vcf files, via -a option. 
Kindred uses multi-threading to do calculation and user can specify the number of threads with -t option. 

## Output
Kindred output two files. One is "pref.grm" the other is "pref.kin", where pref can be specified with -o option. 
In pref.grm, the first line contains individual ID. 
The rest is an $n\times n$ square matrix with 6 digits accuracy. The $i$-th row and $j$-th column is $2\phi_{ij}$. 
In pref.kin, the first line is the header: "ID1 ID2 phi sumd d1 d2 d3 d4 d5 d6 d7 d8 d9", with space delimit. 
ID1 and ID2 are two sample IDs, phi is the kinship, d1, ..., d9 are probabilities of Jacquard latent states, and sumd is the sum of the nine probabilities, and a departure from $1$ suggests ill fit, which happens rarely, but when it happens, it happens on $\phi_{ii}$.   
In addition to header, the file contains $n(n+1)/2$ lines, containing all pairs.   

## Usage examples

1) If in.vcf.gz contains AF field, one can simply do: 

<code> $> kindred in.vcf.gz -o pref </code>

2) If in.vcf.gz has no AF field, and one wants to use allele frequencies estimated from genotypes in the vcf file. 

<code> $> bcftools plugin fill-tags in.vcf.gz | kindred - -o test </code>

3) If in.vcf.gz contain AF field, but you want to recomputed it. 

<code> $> bcftools annotate --remove INFO in.vcf.gz | \
          bcftools plugin fill-tags in.vcf.gz | \
          kindred - -o test </code> 

4) You may store precomputed allele frequencies in a vcf file "annotate.vcf.gz". Kindred can use the allele frequencies by  

<code> $> bcftools annotate --remove INFO -c 'INFO/AF' -a annotate.vcf.gz in.vcf.gz  | \
          kindred -v - -o test </code> 

5) You may store multiple allele frequencies in an annoation file, and you want use EUR_AF instead of EAS_AF or AFR_AF: 
  
<code> $> bcftools annotate --remove INFO -c 'INFO/EUR_AF' -a annotate.vcf.gz in.vcf.gz  | \
          kindred -v - -o test -a EUR_AF </code>

6) If you want to only use markers on chromosome 8 to compute kinship using sample estimated allele frquencies:  

<code> $> bcftools filter -r 8 in.vcf.gz | bcftools plugin fill-tags  | \
          kindred - -o test.chr8 </code>

7) If you want to do it with allele frequencies stored in an annotation file:   

<code> $> bcftools fitler -r 8 in.vcf.gz | \
          bcftools annotate --remove INFO -c 'INFO/AF' -a annotate.vcf.gz  | \
          kindred - -o test.chr8 </code> 


## Protocol to prepare annotation vcf files
To prepare annotation files that contains allele frequencies for a specific population,  
we downloaded 1000 genomes vcf and tbi files, they were arranged in different chromosomes. 
We first obtain all allelic SNPs with minimum of 50 counts of minor alleles (among 2504 samples).  
for each SNP count AN-AC for different populations (the example below is for chr22 of CHB). 

<code>  bcftools view -m2 -M2 -v snps -c 50:minor ALL.chr22.phase3.genotypes.vcf.gz | \
        bcftools view -S samples.CHB.txt | bcftools annotate --remove INFO |\
        bcftools plugin fill-AN-AC | \
        bcftools query -f "%CHROM %POS %REF %ALT %INFO/AC %INFO/AN\n" > chb.chr22.an-ac </code> 

 
Suppose we obtain for each chrosome the an-ac file for populations CEU and YRI, in addition to CHB.  
For each SNP we  do a chisq test to compute p-values (pv), and write a vcf file with INFO/FST = -10log10(pv). 
We can then use FST to filter SNPs that to be annotated. 

To prepare a vcf file with CEU allele frequencies.  

<code> bcftools view -S ceu.samples.list annotate --remove INFO plugin fill-tags | \
 bcftools view -G > annotation.ceu </code> 

## Other useful tips. 
1) plink can convert plink format to vcf format

<code> $> plink --bfile prefix --recode vcf --out in </code>

2) Bcftools requires vcf files to be zipped and indexed

<code> $> bgzip in.vcf </code>
<code> $> tabix in.vcf.gz </code> 

You will have in.vcf.gz (instead of in.vcf) and in.vcf.gz.tbi. 

