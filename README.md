# kindred
(Infer realized kinship via latent identity-by-descent states)

Kindred is designed to infer kinship $\phi$ between a pair of samples. If the pair is one and itself, the kinship is $\phi = (1+F)/2$ where $F$ denotes inbreeding coefficient. 

## Input and options
Kindred takes vcf file as input. The vcf file must contain "INFO/AF" field and genotypes. If the AF field is not readily available, it can be populated by bvftools  on the fly. 
User can also specify other name of the field, such as EUR\_AF or EAS\_AF in 1000 genomes vcf files, via -a option. 
Kindred uses multi-threading to do calculation and one can specify the number of thread with -t option. 

## Output
Kindred output two files. One is "pref.grm" the other is "pref.kin", where pref can be specified with -o option. 
In pref.grm, the first line contains individual ID. 
The rest is an $n\times n$ square matrix with 6 digits accuracy. The $i$-th row and $j$-th column is $2\phi_{ij}$. 
In pref.kin, the first line is the header: "ID1 ID2 phi sumd d1 d2 d3 d4 d5 d6 d7 d8 d9", with space delimit. 
ID1 and ID2 are two sample IDs, phi is the kinship, d1, ..., d9 are probabilities of Jacquard latent states, and sumd is the sum of the nine probabilities, and a departure from $1$ suggests ill fit, which happens rarely, but when it happens, it happens on $\phi_{ii}$.   
In addition to header, the file contains $n(n+1)/2$ lines, containing all pairs.   

## Usage examples

### If in.vcf.gz contains AF field, one can simply do: 
<code> 
$> kindred in.vcf.gz -o pref
</code>




\item If in.vcf.gz has no AF field, and one wants to use allele frequencies estimated from genotypes in the vcf file. 
\begin{verbatim}
$> bcftools plugin fill-tags in.vcf.gz | kindred - -o test
\end{verbatim}

\item If in.vcf.gz contain AF field, but you want to recomputed it. 
\begin{verbatim}
$> bcftools annotate --remove INFO in.vcf.gz | \
 bcftools plugin fill-tags in.vcf.gz | \
 kindred - -o test
\end{verbatim}

\item Suppose you have a precomputed allele frequencies in a vcf file named "annotate.vcf.gz". You can use it to annotate in.vcf.gz so that kindred can use the allele frequencies for computation. 
\begin{verbatim}
$> bcftools annotate --remove INFO -c 'INFO/AF' -a annotate.vcf.gz in.vcf.gz  | \
 kindred -v - -o test
\end{verbatim}

\item Suppose in annotate.vcf.gz you have multiple tags under INFO such as EUR\_AF and EAS\_AF and you want use EUR\_AF in kindred calcultion. 
\begin{verbatim}
$> bcftools annotate --remove INFO -c 'INFO/EUR_AF' -a annotate.vcf.gz in.vcf.gz  | \
 kindred -v - -o test -a EUR_AF
\end{verbatim}
This command first clears INFO fields from in.vcf.gz, then creates EUR\_AF field using information in annotate.vcf.gz.  Kindred takes the output stream to do calculations on the fly.  Those SNPs who are not annotated with an AF field will not go into the kindred calculation. 
%

\item Suppose you want to only use markers on chromosome 8 to do calculation,  you can do 
\begin{verbatim}
$> bcftools filter -r 8 in.vcf.gz | bcftools plugin fill-tags  | \
  kindred - -o test.chr8
\end{verbatim}

\item With annotation file you can do 
\begin{verbatim}
$> bcftools fitler -r 8 in.vcf.gz | \
   bcftools annotate --remove INFO -c 'INFO/AF' -a annotate.vcf.gz  | \
   kindred - -o test.chr8
\end{verbatim}
This will compute kinship using markers on chr10. 
\end{itemize}

\subsection{Protocol to prepare annotation vcf files}
To prepare annotation files that contains allele frequencies for a specific population,  
we downloaded 1000 genomes vcf and tbi files, they were arranged in different chromosomes. 
We first obtain all allelic SNPs with minimum of 50 counts of minor alleles (among 2504 samples).  
for each SNP count AN-AC for different populations (the example below is for chr22 of CHB). 
\begin{verbatim}
 bcftools view -m2 -M2 -v snps -c 50:minor ALL.chr22.phase3.genotypes.vcf.gz | \
 bcftools view -S samples.CHB.txt | bcftools annotate --remove INFO |\
 bcftools plugin fill-AN-AC | \
 bcftools query -f "%CHROM %POS %REF %ALT %INFO/AC %INFO/AN\n" > chb.chr22.an-ac 
 \end{verbatim}
 % bcftools view -m2 -M2 -v snps -c 50:minor raw/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | bcftools view -S samples.CHB.txt
 
Suppose we obtain for each chrosome the an-ac file for populations CEU and YRI, in addition to CHB.  
For each SNP we  do a chisq test to compute p-values (pv), and write a vcf file with INFO/FST = -10log10(pv). 
We can then use FST to filter SNPs that to be annotated. 

 To prepare a vcf file with CEU allele frequencies.  
\begin{verbatim}
bcftools view -S ceu.samples.list annotate --remove INFO plugin fill-tags | \
 bcftools view -G > annotation.ceu 

Below we use bcftools to prepare vcf files, which requires vcf files to be zipped and indexed. 
\begin{verbatim}
$> bgzip in.vcf 
$> tabix in.vcf.gz
\end{verbatim}
You will have in.vcf.gz (instead of in.vcf) and in.vcf.gz.tbi. 

