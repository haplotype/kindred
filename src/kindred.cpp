//////////////////////////////////////////////////////////////////////////////////
//The MIT License (MIT)
//
//Copyright (c) 2023 Yongtao Guan
//  Bug report: ytguan@gmail.com 
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

#include <random> 
#include <iostream> 
#include <fstream> 
#include <algorithm> 
#include <vector> 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 
#include <unistd.h>
//#include <sys/stat.h>
//#include <sys/types.h>
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"              
#include <map> 
#include <math.h>
#include "lsfit.h"
#include <pthread.h>
#include <stdint.h>
#include <zlib.h>
#include "thpool.h"

#include "gsl/gsl_linalg.h"

using namespace std; 

#define VERSION "0.81"
//18 Jan 2023   version 0.7
//01 Feb 2023   version 0.8
//03 Mar 2023   version 0.81
typedef struct KData {
    double * grm; 
    double ** kin;   //kinship 
    uint32_t * gt;    //2 bit hold a genotype
//    double * raf;    //snp allele freq.
    uint16_t * raf;    //snp allele freq.
    int nb;      //number of allele frequency bins. 
    int ni;      //number of individuals
    int ns;      //number of snps. 
    int nth;     //number of threads
    int yesk;    //toggle off write kin. 
    int gsl;     //use gsl fit if gsl == 1; 
    int * beg; 
    int * end; 
    int * ileft; 
    int * iright; 
//    double se; 
    int bychr; 
} KData; 

KData gdadu; 
map<long, long> tr2sq; 
map<long, long> sq2tr; 
//the above is for kindred threads. 

typedef struct TData {
    double ** triM; //nth x upper-triangle matrix  (ni+1) x ni / 2
    vector<uint8_t> v8; 
    vector<float> vaf;
    double * vsumd2; //nth array
} TData; 
TData guru; 
//the additional data structure is for ukin threads. 

//global data for multithreading. 
void jacquard_bychr(void *par);  
void jacquard(void *par);  
void jacquard3(void *par);  
int make_grm(string fn, string pref, int); 
void perSNP_scGRM(void *par); 
void gsl_solve(double ** AtA, double * AtB, double * x, int np);

double ** Allocate2DMatrix(int dim1, int dim2)
{
	int i;
	double ** m;
	
	m = (double **) malloc((size_t)((dim1)*sizeof(double*)));
	m[0] = (double *) malloc((size_t)((dim1*dim2)*sizeof(double)));
	memset(m[0], 0, (dim1*dim2)*sizeof(double)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D double matrix. \n");
		exit(0);
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}

void Free2DMatrix(double ** m)
{
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}

uint8_t ** Allocate2DByteMatrix(int dim1, int dim2)
{
	int i;
	uint8_t ** m;
	
	m = (uint8_t **) malloc((size_t)((dim1)*sizeof(uint8_t *)));
	m[0] = (uint8_t *) malloc((size_t)((dim1*dim2)*sizeof(uint8_t)));
	memset(m[0], 0, (dim1*dim2)*sizeof(uint8_t)); 
	if (!(m && m[0]))
	{
		printf("Error: Problem allocating a 2D double matrix. \n");
		exit(0);
	}
	for(i = 1; i < dim1; i++)
	{
		m[i] = m[i-1] + dim2;
	}
	return (m);
}

void Free2DByteMatrix(uint8_t ** m)
{
	if(m == NULL) return; 
	free(m[0]);
	free(m);
	m = NULL; 
}

void gdm(double p, double * mm) 
{
//    double q = 1 - p; 
//    double p2 = p * p; 
//    double q2 = q * q; 
//    double pq = p * q; 
//    double p3 = p2 * p; 
//    double q3 = q2 * q; 
//    double p2q = p2 * q; 
//    double pq2 = p * q2; 
//
//
//    double tt[9][9] = { {p, p2, p2, p3, p2, p3, p2, p3, p3*p}, 
//	  {0, 0,  pq, 2*p2q, 0, 0, 0, p2q, 2*p2*pq},
//	  {0, pq, 0, pq2, 0, p2q, 0, 0, pq*pq}, 
//	  {0, 0, 0, 0, pq, 2*p2q, 0, p2q, 2*p2*pq}, 
//	  {0, 0, 0, 0, 0, 0, 2*pq, pq, 4*pq*pq}, 
//	  {0, 0, 0, 0, pq, 2*pq2, 0, pq2, 2*q2*pq},
//	  {0, pq, 0, p2q, 0, pq2, 0, 0, pq*pq}, 
//	  {0, 0, pq, 2*pq2, 0, 0, 0, pq2, 2*q2*pq}, 
//	  {q, q2, q2, q3, q2, q3, q2, q3, q3*q} }; 

    double q = 1 - p; 
    double p2 = p * p; 
    double q2 = q * q; 
    double pq = p * q; 
    double p3 = p2 * p; 
    double q3 = q2 * q; 
    double p2q = p2 * q; 
    double pq2 = p * q2; 
    double p4 = p3*p;
    double p3q = p2*pq;
    double p2q2 = pq*pq;
    double pq3 = pq*q2;
    double q4 = q3*q;


    double tt[9][9] = { {p, p2, p2, p3, p2, p3, p2, p3, p4}, 
	  {0, 0,  pq, 2*p2q, 0, 0, 0, p2q, 2*p3q},
	  {0, pq, 0, pq2, 0, p2q, 0, 0, p2q2}, 
	  {0, 0, 0, 0, pq, 2*p2q, 0, p2q, 2*p3q}, 
	  {0, 0, 0, 0, 0, 0, 2*pq, pq, 4*p2q2}, 
	  {0, 0, 0, 0, pq, 2*pq2, 0, pq2, 2*pq3},
	  {0, pq, 0, p2q, 0, pq2, 0, 0, p2q2}, 
	  {0, 0, pq, 2*pq2, 0, 0, 0, pq2, 2*pq3}, 
	  {q, q2, q2, q3, q2, q3, q2, q3, q4} }; 

    for (int i = 0, k=0; i < 9; i++)
	for (int j = 0; j < 9; j++,k++)
	    mm[k] = tt[i][j]; 
}

void igdm(double a, double b, double * mm) //assume p is uniform on (a,b)
{
//    double q = 1 - p; 
//    double p2 = p * p; 
//    double q2 = q * q; 
//    double pq = p * q; 
//    double p3 = p2 * p; 
//    double q3 = q2 * q; 
//    double p2q = p2 * q; 
//    double pq2 = p * q2; 
//    double p4 = p3*p;
//    double p3q = p2*pq;
//    double p2q2 = pq*pq;
//    double pq3 = pq*q2;
//    double q4 = q3*q;

    double a2 = a*a; 
    double a3 = a2*a;
    double b2 = b*b; 
    double b3 = b2*b; 
    double i1 = (b-a); 
    double ip = (b2-a2)/2;
    double iq = i1 - ip; 
    double ip2 = (b3-a3)/3; 	    
    double iq2 = i1 - 2* ip + ip2;  
    double ipq = ip2 - ip;          
    double ip3 = (b3*b - a3*a)/4;   
    double iq3 = i1 - 3*ip + 3*ip2 - ip3;     
    double ip2q = ip2 - ip3;   
    double ipq2 = ip -2*ip2 + ip3;
    double ip4 = (b3*b2 - a3*a2)/5;    
    double ip3q = ip3 - ip4; 
    double ip2q2 = ip2 - 2*ip3 + ip4; 
    double ipq3 = ip-3*ip2+3*ip3-ip4;    
    double iq4 = i1-4*ip + 6*ip2 - 4*ip4 + ip4;   
    ip /=i1;
    iq /=i1;
    ip2 /=i1;     
    iq2 /=i1;
    ipq /=i1;
    ip3 /= i1; 
    iq3 /= i1; 
    ip2q /=i1;
    ipq2 /=i1;
    ip4 /= i1;
    ip3q /= i1; 
    ip2q2 /= i1;
    ipq3 /= i1;
    iq4 /= i1; 
    double tt[9][9] = { {ip, ip2, ip2, ip3, ip2, ip3, ip2, ip3, ip4}, 
	  {0, 0,  ipq, 2*ip2q, 0, 0, 0, ip2q, 2*ip3q},
	  {0, ipq, 0, ipq2, 0, ip2q, 0, 0, ip2q2}, 
	  {0, 0, 0, 0, ipq, 2*ip2q, 0, ip2q, 2*ip3q}, 
	  {0, 0, 0, 0, 0, 0, 2*ipq, ipq, 4*ip2q2}, 
	  {0, 0, 0, 0, ipq, 2*ipq2, 0, ipq2, 2*ipq3},
	  {0, ipq, 0, ip2q, 0, ipq2, 0, 0, ip2q2}, 
	  {0, 0, ipq, 2*ipq2, 0, 0, 0, ipq2, 2*ipq3}, 
	  {iq, iq2, iq2, iq3, iq2, iq3, iq2, iq3, iq4} }; 

    for (int i = 0, k=0; i < 9; i++)
	for (int j = 0; j < 9; j++,k++)
	    mm[k] = tt[i][j]; 
}

void print_progress_num(int last, int p)
{
    if(!last) {
	printf("##processed variants = %d \r", p); 
	fflush(stdout); 	
    } else {
	printf("##processed variants = %d \n", p); 
    }
}

void print_progress_bar(int last, long p, long total)
{
    char str[] = "##processed pairs:"; 
    int progress = (int) (100.0 * p / total); 
//    int barsize = (int) (progress / 2.0); 
//    char bar[100];
//    memset(bar, '\0', 100); 
    if(!last) {
//	    for (int i = 0; i < barsize; i++)
//		    bar[i] = '>'; 
	    printf("%s %ld or %d%%\r", str, p, progress); 
	    fflush(stdout); 	
    } else {
//	    for (int i = 0; i < barsize; i++)
//		    bar[i] = '>'; 
//	    printf("%s [%-50s] 100%%\n", str, bar); 
	    printf("%s %ld or %d%%\n", str, p, progress); 
    }
}

int usage()
{
	fprintf(stdout, "\n");
	fprintf(stdout, "Kindred (infer realized kinship and inbreeding coefficients)\n");
	fprintf(stdout, "Version: %s\n", VERSION);
	fprintf(stdout, "Usage:   kindred -i in.vcf.gz [-adfgjkotu]\n");
	fprintf(stdout, "Options: \n");
	fprintf(stdout, "         -a str        INFO tag for allele frequency [AF]\n");
//	fprintf(stdout, "         -b int        number of allele freq bins (max 128) [100]\n");
	fprintf(stdout, "         -d int        thin marker by a min neighor distance [1] \n");
//	fprintf(stdout, "         -e int        effective n for AF estimates [100]\n");    
	fprintf(stdout, "         -f flt        minor allele freq threshold [0.01]\n");    
	fprintf(stdout, "         -g str        convert square rkm to binary gcta format\n");
	fprintf(stdout, "         -i str        (indexed) bgzipped vcf file\n");
	fprintf(stdout, "         -j num        number of Jacquard states (3 or 9) [9]\n");
	fprintf(stdout, "         -k            toggle on writing .kin file [off]\n");
//	fprintf(stdout, "         -n str        numeric genotype input file \n");
	fprintf(stdout, "         -o str        output prefix [out]\n");
	fprintf(stdout, "         -t int        number of thread [8]\n");
	fprintf(stdout, "         -u            UKin \n\n");
	fprintf(stdout, "User Manual can be found at github.com/haplotype/kindred\n");
	fprintf(stdout, "Bug report: Yongtao Guan <ytguan@gmail.com>\n\n");
	return 1;
}


int main(int argc, char *argv[])
{
    time_t time0, time1, time2; 
    string fn_vcf; 
    string pref("out");
    string af_tag("AF"); 
    char c;
    int flag_vcf = 0; 
    int flag_gcta = 0; 
    int flag_scGRM = 0;    //flag_ukin; 
    int flag_help = 0;
    int flag_gt012 = 0; //input is 012 ns x ni matrix; 

    string fn_grm; 
    string fn_gt; 
    int snp_dist = 0; 

    //optimized values, do not change; 
    int m_nb = 100;   
    double m_se = 0.04;   
    double m_maf = 0.01; 

    int m_nth = 8; 
    int m_jstates = 9; //3 or 9
    int m_bychr = 0; // 0 or 1 or 2; 
    int m_yesk = 0;  // 0 or 1; 
    int m_gsl = 0;   // 0 or 1;



    while ((c = getopt(argc, argv, "a:b:c:d:e:f:g:hi:j:kn:o:t:uy")) >= 0) 
    {
	switch (c) {
	    case 'a': 
		af_tag.assign(optarg); 
		break; 
	    case 'b': 
		m_nb = atoi(optarg); 
		break;
	    case 'c': 
		m_bychr = atoi(optarg); 
		break;
	    case 'd': 
		snp_dist = atoi(optarg); 
		break;
	    case 'e': 
		m_se = atof(optarg); 
		break;
	    case 'f': 
		m_maf = atof(optarg); 
		break;
	    case 'g': 
		flag_gcta = 1; 
		fn_grm.assign(optarg); 
		break;
	    case 'h': 
		flag_help = 1; 
		break; 
	    case 'i': 
		flag_vcf = 1; 
		fn_vcf.assign(optarg); 
		break;
	    case 'j': 
		m_jstates = atoi(optarg); 
		break;
	    case 'k': 
		m_yesk = 1 - m_yesk; 
		break;
	    case 'n': 
		flag_gt012 = 1; 
		fn_gt.assign(optarg); 
		break;
	    case 'o': 
		pref.assign(optarg); 
		break;
	    case 't': 
		m_nth = atoi(optarg); 
		break;
	    case 'u': 
		flag_scGRM = 1; 
		break;
	    case 'y': 
		m_gsl = 1; 
		break;
	    default: 
		break; 
	}
    }

//	std::cout << argc << " " << optind << std::endl; 
//    if(argc > optind) 
//	fn_vcf.assign(argv[argc-1]); 
    //this seems not working with mac version of getopt. 
    if(flag_gcta == 1) return make_grm(fn_grm, pref, 0); 
    if ((flag_vcf == 0 && flag_gt012 == 0) || flag_help == 1) return usage();
    fprintf(stdout, "\n"); 
    fprintf(stdout, "Kindred %s by Yongtao Guan at Framingham Heart Study \n", VERSION); 
    fprintf(stdout, "  National Heart, Lung, and Blood Institute (C) 2023 \n"); 
    if(!(m_jstates == 3 || m_jstates ==9)) 
    {
	fprintf(stdout, "number of latent states has to be eitehr 9 or 3\n"); 
	return usage(); 
    }

    gdadu.nb = m_nb;  
    gdadu.nth = m_nth; 
    gdadu.yesk = m_yesk; 
    gdadu.gsl = m_gsl; 
    gdadu.bychr = m_bychr; 


    gzFile fp1 = NULL;
    gzFile fp2 = NULL;
    gzFile fp3 = NULL; 
    FILE * fplog = NULL;
    string buf; 
    buf.assign(pref); 
    buf.append(".log");  
    fplog = fopen(buf.c_str(), "w");
    if (fplog == NULL) {
	    fprintf(stderr, "can't open file %s\n", buf.c_str()); 
	    exit(EXIT_FAILURE);
    }   // for SNPs. 
    for (int i = 0; i < argc; i++)
	fprintf(fplog, "%s ", argv[i]); 
    fprintf(fplog, "\n"); 
    fprintf(fplog, "##m_nb = %d \n", m_nb); 
    fprintf(fplog, "##m_se = %0.3f \n", m_se); 
    fprintf(fplog, "##m_maf = %0.3f \n", m_maf); 
    fprintf(fplog, "##m_nth = %d \n", m_nth); 
    fprintf(fplog, "##m_jstates = %d \n", m_jstates); 
    fprintf(fplog, "##m_bychr = %d \n", m_bychr); 
    fprintf(fplog, "##m_yesk = %d \n", m_yesk); 
    fprintf(fplog, "##m_gsl = %d \n", m_gsl); 

    time0 = time(NULL); 
//    vector<float> vaf; 
//    vector<uint8_t> v8; 
    //open vcf file; 

    gdadu.ileft = new int[m_nb]; 
    gdadu.iright = new int[m_nb]; 
    for (int i = 0; i < m_nb; i++)
    {
//	double p = double(i)/nb;      
//	double d = 1.96 * sqrt(p * (1-p) / ne); 
//	int wind =  floor(d * nb); 
	int wind = floor(m_se * m_nb); 
	int left = i - wind;
	if(left < 0) left = 0; 
	int right = i + wind;
	if(right > (m_nb-1)) right = m_nb-1; 
	gdadu.ileft[i] = left;  
        gdadu.iright[i] = right; 
    }   //this defines intervals to count freq. and the center is just p = 0.02 ... 0.98. 
//    for (int i = 0; i < m_nb; i+=m_nb/10)
//	fprintf(stdout, "%d < %d > %d \n", gdadu.ileft[i], i, gdadu.iright[i]);

    int snps_by_chr[22]; 
    gdadu.beg = new int[22]; 
    gdadu.end = new int[22]; 
    for (int i = 0; i < 22; i++)
	gdadu.beg[i] = gdadu.end[i] = snps_by_chr[i] = 0; 
    int ns1 = 0; 
    if(flag_vcf == 1) 
    {
	htsFile *fpv = hts_open(fn_vcf.c_str(), "r");   
	bcf_hdr_t *hdr = bcf_hdr_read(fpv);  
	bcf1_t* line = bcf_init();   
	gdadu.ni = bcf_hdr_nsamples(hdr); 
	fprintf(fplog, "##number of samples: %d \n", gdadu.ni); 
	fprintf(stdout, "##number of samples: %d \n", gdadu.ni); 
    //	if(hdr->samples != NULL) {
    //	    for (int i = 0; i < gdadu.ni; i++)       
    //		cout << hdr->samples[i] << endl; 
    //	}
	gdadu.ns = 0; 
	int ns0 = 0;  
	int last_chr = -1; 
	int last_pos = 0; 
	while(bcf_read1(fpv, hdr, line) == 0) {   

	    if(ns1 % 1000 == 0) 
	        print_progress_num(0, ns1);  
	    ns1++; 

	    bcf_get_variant_types(line); 
	    if(line->d.var_type != VCF_SNP) continue; 
	    if(line->n_allele > 2) continue;       
	    if(snp_dist > 0) {
               if((last_chr != line->rid) || (line->pos - last_pos > snp_dist))
	       {
		   last_chr = line->rid; 
		   last_pos = line->pos; 
	       }
	       else 
		   continue; 
	    }  //this is to thin SNPs by bp distance.  simple, but not optimal. 

	    int dp_status = 0; 
	    float * dst = NULL;
	    int ndst = 0;
	    double af = 0; 
	    dp_status = bcf_get_info_values(hdr, line, af_tag.c_str(),  (void**) &dst, &ndst, BCF_HT_REAL);
	    if(dp_status > 0) {
		af = dst[0]; 
		free(dst); 
		if(af < m_maf || af > 1-m_maf)
		    continue; 
		//filter by maf. 
	    } else {
		ns0++; 
		free(dst); 
		continue; 
		//filter by af_tag. 
	    }
			 
	    int32_t *gt_arr = NULL, ngt_arr = 0;
            int ngt;
	    ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
	    if ( ngt<=0 ) 
	    {
		cout << "no genotypes at " << bcf_seqname(hdr, line) << " " << line->pos << endl; 
		continue;
		//fiter by gt presence. 
	    }

	    guru.vaf.push_back(af); 
	    for (int i=0; i<gdadu.ni; i++)
	    {
		float gt = 3; 
		if (gt_arr[2*i+0] != bcf_gt_missing && gt_arr[2*i+1] != bcf_gt_missing)
		    gt = bcf_gt_allele(gt_arr[2*i+0]) + bcf_gt_allele(gt_arr[2*i+1]); 
		guru.v8.push_back(gt); 
	   }
	   free(gt_arr);
	   gdadu.ns++; 
           snps_by_chr[line->rid]++;  
	   //all satisfied, taking in the data; 
       }

       fprintf(fplog, "##number of biallelic SNPs: %d \n", gdadu.ns); 
       fprintf(fplog, "##number of biallelic SNPs without tag: %s = %d \n", af_tag.c_str(), ns0); 
       fprintf(stdout, "##number of biallelic SNPs: %d \n", gdadu.ns); 
       fprintf(stdout, "##number of biallelic SNPs without tag: %s = %d \n", af_tag.c_str(), ns0); 
       for (int i = 0; i < 22; i++)
       	   fprintf(fplog, "%d %d \n", i+1, snps_by_chr[i]); 
       print_progress_num(1, ns1);  
       fflush(fplog); 
       hts_close(fpv); 
       if(gdadu.ns == 0 || ns0 == gdadu.ns) 
       {
	   fprintf(fplog, "##no SNPs with AF tag. abort. \n"); 
	   fprintf(stdout, "##no SNPs with AF tag. abort. \n"); 
	   fflush(fplog); 
           fclose(fplog); 
	   exit(0); 
       }
    }

    gdadu.beg[0] = 0; 
    gdadu.end[0] = snps_by_chr[0]; 
    for (int i = 1; i < 22; i++) 
    {
	gdadu.end[i] = gdadu.end[i-1] + snps_by_chr[i]; 
	gdadu.beg[i] = gdadu.end[i-1]; 
    }
    //end is a cumulative sum of snp_by_chr; 

    if(flag_gt012 == 1)  //this is to read in genotype 012 matrix for debugging. 
    {
	ifstream fin; // indata is like cin
	float num; 
	fin.open(fn_gt); // opens the file
	if(!fin) { // file couldn't be opened
	  cerr << "Error: file could not be opened" << endl;
	  exit(1);
	}                     
	string line; 

	gdadu.ns = 0; 
        while (std::getline(fin, line))
             ++gdadu.ns;
        //cout lines; 
	fin.clear();
	fin.seekg(0, ios::beg);
	//rewind; 

	fin >> num;
	while ( !fin.eof() ) { // keep reading until end-of-file
	//      cout << "The next number is " << num << endl;
	      guru.v8.push_back(num); 
	  fin >> num; // sets EOF flag if no value found
	}
	fin.close();

	gdadu.ni = guru.v8.size() / gdadu.ns; 
	for(int m = 0; m < gdadu.ns; m++)
	{
	    float sum = 0; 
	    for (int i = 0; i < gdadu.ni; i++)
	    {
		float gt = guru.v8.at(m*gdadu.ni+i); 
//		if(m == 0) 
//		    cout << gt << " " << endl;  
		sum += gt; 
	    }
	    float af = sum/(2.0*gdadu.ni); 
//	    cout << af << endl; 
	    guru.vaf.push_back(af); 
	}
//	cout << guru.v8.size() << endl; 
	cout << gdadu.ni << endl; 
	cout << gdadu.ns << endl; 
//	cout << guru.vaf.size() << endl; 
    }

    if(flag_scGRM == 1)  //ukin method
    {
	int ni = gdadu.ni; 
	int ns = gdadu.ns; 
	double ** scM = Allocate2DMatrix(ni,ni); 
	for (int i = 0; i < ni; i++)
	    for (int j = 0; j < ni; j++)
		scM[i][j] = 0; 
//   	compute_scGRM(scM); 
	long len = (ni+1)*ni/2; 
	guru.triM = Allocate2DMatrix(gdadu.nth, len); 
	for (int i = 0; i < gdadu.nth; i++)
	    for (int j = 0; j < len; j++)
		guru.triM[i][j] = 0; 

	guru.vsumd2 = new double[gdadu.nth]; 
	for (int i = 0; i < gdadu.nth; i++)
	    guru.vsumd2[i] = 0; 

////////////////////////////////////////////////////////////
	threadpool thpool = thpool_init(gdadu.nth);
	for(long m = 0; m < (long) ns; m++)
	{
	    thpool_add_work(thpool, perSNP_scGRM, (void*)m);
	    if(m % 1000 == 0) {
	       thpool_wait(thpool); 
	       print_progress_num(0, m); 
	    }
	}
	print_progress_num(1, (long) ns); 

	thpool_wait(thpool);
	thpool_destroy(thpool);
///////////////////////////////////////////////////////////

	for (int i = 0; i < ni; i++)
	    for (int j = i; j < ni; j++)
		scM[i][j] = 0; 
	double localsumd2 = 0; 
	for (int p = 0; p < gdadu.nth; p++)
	{
	    for (int i = 0, k = 0; i < ni; i++)
		for (int j = i; j < ni; j++, k++)
		    scM[i][j] += guru.triM[p][k]; 
	    localsumd2 += guru.vsumd2[p]; 
	}
//	cout << localsumd2 << endl; 
//	for (int i = 0; i < ni; i++)
//	    cout << scM[0][i] << " "; 
//	cout << endl; 
	for (int i = 0; i < ni; i++)
	    for (int j = i; j < ni; j++)
	    {
		scM[i][j] /= localsumd2; 
		scM[j][i] = scM[i][j]; 
	    }

	double * row = new double[ni]; 
	for (int i = 0; i < ni; i++)
	{
	    row[i] = 0; 
	    for (int j = 0; j < ni; j++)
		if(j != i) row[i] += scM[i][j]; 
	}
	for (int i = 0; i < ni; i++)
	    row[i] /= 2.0;
	
	for (int i = 0; i < ni; i++)
	    for (int j = 0; j < ni; j++)
		scM[i][j] += (row[i] + row[j] + 1);
	delete[] row; 

	delete[] guru.vsumd2; 
	Free2DMatrix(guru.triM); 

	buf.assign(pref); 
	buf.append(".ukinGRM.gz");  
	fp3 = gzopen(buf.c_str(), "w");
	if (fp3 == NULL) {
		fprintf(stderr, "can't open file %s\n", buf.c_str()); 
		exit(EXIT_FAILURE);
	}   // for SNPs. 

	for (int ii = 0; ii < ni; ii++)
	{
	   for (int jj = 0; jj < ni; jj++)
	   {
	       gzprintf(fp3, "%8.6f ", scM[ii][jj]); 
	   }
	   gzprintf(fp3, "\n"); 
	}
	gzclose(fp3); 
        return 1;
    }


   gdadu.raf = new uint16_t[gdadu.ns]; //rounded allele frequency. 
   for (int m = 0; m < gdadu.ns; m++)
   {
       gdadu.raf[m] = floor(gdadu.nb * guru.vaf.at(m));  
       if(gdadu.raf[m] == gdadu.nb) 
	   gdadu.raf[m] = gdadu.nb-1; 
   }
   
//   cout << gdadu.ns << endl; 
//   cout << gdadu.ni << endl; 
   long blen = ((gdadu.ns >> 2) +1) * ((gdadu.ni >> 2) +1) ; 
//   cout << blen << endl;
   gdadu.gt = new uint32_t[blen]; 
   for (long i = 0; i < blen; i++)
       gdadu.gt[i] &= ~0xFFFFFFFF; 
   unsigned long qk = 0; 
   uint8_t rk = 0; 
   for (int i = 0; i < gdadu.ni; i++)
   {
       unsigned long pos = i; 
//       cout << pos << endl; 
       for (int j = 0; j < gdadu.ns; j++)
       {
//	   gdadu.gt[k/16] |= ((v8.at(j*gdadu.ni+i) & 3) << ((k%16)*2)); 
	   gdadu.gt[qk] |= ((guru.v8.at(pos) & 3) << ((rk)*2)); 
	   rk++;  if((rk & 0xF) == 0) { rk &= ~0xFF; qk++;}       
	   pos += gdadu.ni; 
       }
   }
   //using two bits to store one genotype. 11 is for missing. 
   //qk is quotient k/16; rk is remainder k%16. this is faster. 
   fprintf(fplog, "##stored genotypes in %ld uint32_t\n", qk);         

   vector<uint8_t> ().swap(guru.v8); 

    gdadu.grm = new double[gdadu.ni*(gdadu.ni+1)/2];
    if(gdadu.yesk == 1) 
    	gdadu.kin = Allocate2DMatrix(gdadu.nth*1000,m_jstates+2); 
    else 
	gdadu.kin = NULL; 

   for (int i=0, k=0; i<gdadu.ni; i++) 
       for (int j=i; j<gdadu.ni; j++, k++) 
       {
		   long tr = (long) k; 
		   long sq = (long) (i*gdadu.ni+j); 
		   tr2sq[tr] = sq; 
		   sq2tr[sq] = tr;
       }
   //global map for exchange indices; upper triangular to full matrix. 

   time1 = time(NULL); 

//   double wt[9] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1}; 
   if(gdadu.yesk == 1) 
   {
	string buf; 
	buf.assign(pref); 
	buf.append(".kin.gz");  
	fp1 = gzopen(buf.c_str(), "w");
	if (fp1 == NULL) {
		fprintf(stderr, "can't open file %s\n", buf.c_str()); 
		exit(EXIT_FAILURE);
	}   // for SNPs. 
   	gzprintf(fp1, "ID1 ID2 phi sumd d1 d2 d3 d4 d5 d6 d7 d8 d9\n"); 
   }

    buf.assign(pref); 
    buf.append(".rkm.gz");  
    fp2 = gzopen(buf.c_str(), "w");
    if (fp2 == NULL) {
	    fprintf(stderr, "can't open file %s\n", buf.c_str()); 
	    exit(EXIT_FAILURE);
    }   // for SNPs. 

   //if i save upper triangular matrix linearly. 
   //how to efficienty recover coordinates from linear poistion?  
    //for (int i =0, k=0; i<ni; i++)
    //    for (int j = i; j< ni; j++, k++)
    //          square[i][j] = linear[k]; 


    fprintf(fplog, "##init thread pool n = %d\n", gdadu.nth); 
    fprintf(stdout, "##init thread pool n = %d\n", gdadu.nth); 
    threadpool thpool = thpool_init(gdadu.nth);

    int load = 0; 
    long tiktok = 0; 
    long total = gdadu.ni * (gdadu.ni+1) / 2; 
    for (int ii = 0; ii < gdadu.ni; ii++)
       for (int jj = ii; jj < gdadu.ni; jj++)
       {
	   if(gdadu.bychr > 0) 
	   	thpool_add_work(thpool, jacquard_bychr, (void*)tiktok);
	   else if(m_jstates == 3) 
	   	thpool_add_work(thpool, jacquard3, (void*)tiktok);
	   else //if(m_jstates == 9) 
	   	thpool_add_work(thpool, jacquard, (void*)tiktok);
	   load++; 
	   tiktok++; 
	   if(load == gdadu.nth * 1000) {
	       thpool_wait(thpool); 
//	       fprintf(stdout, "processed %ld many pair of samples \n", tiktok); 
	       print_progress_bar(0, tiktok, total); 
	       if(gdadu.yesk == 1) 
	       {
		   for (long i = tiktok - load;  i < tiktok; i++)
		   {
		       long idx = tr2sq[i]; 
		       int i1 = idx / gdadu.ni; 
		       int i2 = idx % gdadu.ni; 
//		       if(hdr->samples != NULL) 
//			    fprintf(fp1, "%s %s ", hdr->samples[i1], hdr->samples[i2]); 
//		       else 
			    gzprintf(fp1, "%d %d ", i1+1, i2+1);  //using integer is R friendly; sample ID can be obtain by bcftools. 
		       for (int c = 0; c < m_jstates+2; c++)
			   gzprintf(fp1, "%.6f ", gdadu.kin[i%(gdadu.nth*1000)][c]); 
		       gzprintf(fp1, "\n"); 
		   }
	       }
	       load = 0; 
	   }
       }
       print_progress_bar(1, tiktok, total); 

///////////////////////////////////////////////////////
    thpool_wait(thpool);
    thpool_destroy(thpool);
    fprintf(fplog, "##drain thread pool. \n"); 
    fprintf(stdout, "##drain thread pool.\n"); 
///////////////////////////////////////////////////////
    if(gdadu.yesk == 1) 
    {
       for (long i = tiktok - load;  i < tiktok; i++)
       {
	   long idx = tr2sq[i]; 
	   int i1 = idx / gdadu.ni; 
	   int i2 = idx % gdadu.ni; 
//	   if(hdr->samples != NULL) 
//		fprintf(fp1, "%s %s ", hdr->samples[i1], hdr->samples[i2]); 
//	   else 
		gzprintf(fp1, "%d %d ", i1+1, i2+1); 
	   for (int c = 0; c < m_jstates+2; c++)
	       gzprintf(fp1, "%.6f ", gdadu.kin[i%(gdadu.nth*1000)][c]); 
	   gzprintf(fp1, "\n"); 
       }
       gzclose(fp1); 
       string buf1(pref); 
       buf1.append(".kin.gz"); 
       fprintf(stdout, "##wrote %s.\n", buf1.c_str()); 
       fprintf(fplog, "##wrote %s.\n", buf1.c_str()); 
    }

    for (int ii = 0; ii < gdadu.ni; ii++)
    {
       for (int jj = 0; jj < gdadu.ni; jj++)
       {
	   long idx = 0; 
	   if(ii > jj)  idx = sq2tr[(long) jj*gdadu.ni+ii]; 
	   else idx = sq2tr[(long) ii*gdadu.ni+jj]; 
//	   fprintf(fp2, "%8.6f ", gdadu.grm[idx]); 
//	   if(gdadu.grm[idx] < 1e-6) 
//	        gzprintf(fp2, "0 "); 
//	   else 
	   	gzprintf(fp2, "%.6f ", gdadu.grm[idx]); 
       }
       gzprintf(fp2, "\n"); 
    }
    gzclose(fp2); 
    string buf2(pref); 
    buf2.append(".rkm.gz"); 
    fprintf(stdout, "##wrote %s.\n", buf2.c_str()); 
    fprintf(fplog, "##wrote %s.\n", buf2.c_str()); 

    time2 = time(NULL); 
    time_t dt1 = difftime (time1, time0); 
    time_t dt2 = difftime (time2, time1); 
    fprintf(fplog, "##read vcf < %ld seconds\n", dt1+1); 
    fprintf(fplog, "##compute kinship < %ld seconds\n", dt2+1); 
    fprintf(stdout, "##read vcf < %ld seconds\n", dt1+1); 
    fprintf(stdout, "##compute kinship < %ld seconds\n\n", dt2+1); 


    Free2DMatrix(gdadu.kin); 
    delete[] gdadu.grm; 
    delete[] gdadu.gt; 
    delete[] gdadu.raf; 
    delete[] gdadu.ileft; 
    delete[] gdadu.iright; 
    return 0;
}

//work-horse for UKin.
void perSNP_scGRM(void * arg) 
{
    long m = (long) arg; 
    int ni = gdadu.ni; 

    long beg = m * ni; 
    {
	double mu = (double) guru.vaf.at(m)*2; 
	double sum1 = 0; 
	double sum2 = 0; 
	for (int i = 0; i < ni; i++)
	{
	    double ti = (double) guru.v8.at(beg+i);
	    if(ti > 2) {
		sum1 += mu; 
		sum2 += mu * mu; 
	    } else {
		sum1 += ti; 
		sum2 += ti * ti; 
	    }
	}
	sum2 /= ni;
	sum1 /= ni;
	double s2 = (sum2 - sum1 * sum1) * ni / (ni - 1.0); 
	double d2 = mu*(2-mu) - s2 * (1.0 - 1.0/ni); 

	double cotable[6]; 
	cotable[0] = mu * mu; 
	cotable[1] = -mu * (1-mu); 
	cotable[2] = -mu * (2-mu); 
	cotable[3] = (1-mu) * (2-mu); 
	cotable[4] = (2-mu) * (2-mu);
	cotable[5] = (1-mu) * (1-mu); 

	    
	for (int i = 0, k=0; i < ni; i++)
	    for (int j = i; j < ni; j++, k++)
	    {
		uint8_t ti = guru.v8.at(beg+i);
		uint8_t tj = guru.v8.at(beg+j);
		if(ti < 3 && tj < 3) { 
		    uint8_t key = ti + tj; 
//		    cout << key << "\t " << cotable[key] << endl; 
		    if((key == 2) && (ti == tj))  
			key = 5; 
		    guru.triM[m%gdadu.nth][k] += cotable[key]; 
		}
	    }
	guru.vsumd2[m%gdadu.nth] += d2; 
    }
}

void jacquard(void *par) 
{
    int nb = gdadu.nb; 
//    assert(nb==100); 
    int ns = gdadu.ns; 
    int ni = gdadu.ni; 

    long tiktok = (long) par; 
    long kk = tr2sq[tiktok]; 
    int ii = kk / ni; 
    int jj = kk % ni; 
   
//   double * m3 = new double[81]; 
   double * m1 = new double[81]; 
   double ** sigma = Allocate2DMatrix(9, 9); 
   double * theta = new double[10]; 
   double ** count9 = Allocate2DMatrix(nb,10); 

   int np = 9; 
   double ** AtA = Allocate2DMatrix(9,9); 
   double * AtB = new double[9]; 
   double w[9];
   double zz[9];
   int index[9]; 
   double r2; 
   double x[9]; 

//   clock_t clk0 = clock(); 
//	   count; 
   unsigned long ki = (long) ns * ii; 
   unsigned long kj = (long) ns * jj; 
   unsigned long qi = ki / 16;   
   uint8_t ri = ki % 16; 
   unsigned long qj = kj / 16; 
   uint8_t rj = kj % 16; 

   for (int n = 0; n < nb; n++) 
       for (int j = 0; j < 10; j++)
       	   count9[n][j] = 0; 

   for (int m = 0; m < ns; m++)
   {
//       uint8_t gi = (gdadu.gt[ki/4] >> ((ki%4)*2)) & 3; 
//       uint8_t gj = (gdadu.gt[kj/4] >> ((kj%4)*2)) & 3; 
       uint8_t gi = (gdadu.gt[qi] >> ((ri)*2)) & 3; 
       uint8_t gj = (gdadu.gt[qj] >> ((rj)*2)) & 3; 
       ri++; if((ri& 0xF) == 0) {ri &= ~0xFF; qi++;} 
       rj++; if((rj& 0xF) == 0) {rj &= ~0xFF; qj++;} 

       if((gi^3) == 0 || (gj^3) == 0) continue;     
       //11 is missing

       uint16_t tn = gdadu.raf[m]; 
//       if(tn>nb) { 
//	   tn = nb+nb-tn; 
//	   gi = 2 - gi; 
//	   gj = 2 - gj; 
//       }
       uint8_t cc = ((gi << 2) | gj) - gi; 
       // 00 00    0  -0
       // 00 01    1  -0
       // 00 10    2  -0
       // 01 00    4  -1
       // 01 01    5  -1
       // 01 10    6  -1
       // 10 00    8  -2
       // 10 01    9  -2
       // 10 10    10 -2

       count9[tn][8-cc]++;  
       count9[tn][9]++;  
   }
           
   double notwin = 0; 
   for (int n = 0; n < nb; n++) 
   {
//       int set[] = {1,2,3,5,6,7};
       for(int j = 0; j < 3; j++)
	   notwin += count9[n][j+1]; 
   }
   //count genoytpes that are not AA AA, AB AB and BB BB. 
   //if the count notwin === 0, then it's twin or self. 
//   clock_t clk1 = clock(); 
//	   build_matrix; 
//   if(ii == 0 && jj == 1) 
//   {
//       for (int n = 0; n < nb; n++)
//       {
//	   for (int j = 0; j < 10; j++)
//	       cout << count9[n][j] << " \t"; 
//	   cout << endl; 
//       }
//   }

   for (int i = 0; i < 9; i++)
       for (int j = 0; j < 9; j++)
	   AtA[i][j] = 0; 
   for (int j = 0; j < 9; j++)
       AtB[j] = 0; 

   for (int n = 0; n < nb; n++)
   {
       for (int j = 0; j < 10; j++)
       {
	   theta[j] = 0; 
	   for (int k = gdadu.ileft[n]; k <= gdadu.iright[n]; k++)
	       theta[j] += count9[k][j]; 
       }
       double af = 0; 
       double cc[9] = {4,3,2,3,2,1,2,1,0}; 
       if(theta[9] > 0) {
	   for (int j = 0; j < 8; j++)
	       af += theta[j] * cc[j]; 
	   af /= theta[9]; 
	   af /= 4; 
       }
       gdm(af, m1); 
       //use SNPs in a bin to re-compute allele frequency. 
       //this makes it robust to allele freq. mis-specification. 

       for (int i = 0; i < 9; i ++)
	   for (int j = 0; j < 9; j++)
	   {
	       sigma[i][j] = 0; 
	       for (int k = 0; k < 9; k++)
	       	   sigma[i][j] += m1[k*9+i] * m1[k*9+j]; 
	       AtA[i][j] += sigma[i][j] * theta[9]; 
	   }

       for (int j = 0; j < 9; j++)
	   for (int k = 0; k < 9; k++)
	   	AtB[j] += m1[k*9+j] * theta[k]; 
   }

//   if(ii == 1 && jj == 2) 
//   {
//       for (int n = 0; n < nb; n++)
//       {
//	   for (int j = 0; j < 10; j++)
//	       cout << theta[n*10+j] << " \t"; 
//	   cout << endl; 
//       }
//   }

   if(gdadu.gsl == 1) 
       gsl_solve(AtA, AtB, x, np); 
   else 
       nnls(AtA, np, np, AtB, x, &r2, w, zz, index);  
//       cout << "nnls used more than 200 iterations to converge" << endl; 
   double sumx = 0; 
   for (int i = 0; i < 9; i++)
       sumx += x[i]; 

   if(gdadu.yesk == 1) 
   {
       double phi=x[0]+ (x[2]+x[4]+x[6])/2+x[7]/4; 
       int idx = tiktok % (gdadu.nth*1000); 
       gdadu.kin[idx][0] = phi; 
       gdadu.kin[idx][1] = sumx; 
       for (int i = 0; i < 9; i++)
	    gdadu.kin[idx][2+i] = x[i]; 
   }   //pref.kin.gz has unnormalized x; 

   static int warning = 0; 
   if(sumx > 1.02 || sumx < 0.98) 
   {
       if(warning == 0) {
	   fprintf(stdout, "sumd %8.6f deviates from 1, kindred renomralize \n. use -k and examine .kin.gz file for un-renomorlized x.\n", sumx);  
	   warning = 1; 
       }
   }
   for (int i = 0; i < 9; i++)
       x[i] /= sumx; 
   double phi=x[0]+ (x[2]+x[4]+x[6])/2+x[7]/4; 
   if(notwin < 1 && phi < 0.5) 
       phi = 0.5; 
   gdadu.grm[tiktok] = 2*phi; 


   Free2DMatrix(AtA); 
   delete[] AtB; 
   delete[] m1; 
   Free2DMatrix(sigma); 
   delete[] theta; 
   Free2DMatrix(count9); 
}

void jacquard_bychr(void *par) 
{
    int nb = gdadu.nb; 
    int ns = gdadu.ns; 
    int ni = gdadu.ni; 

    long tiktok = (long) par; 
    long kk = tr2sq[tiktok]; 
    int ii = kk / ni; 
    int jj = kk % ni; 
   
   double * m1 = new double[81]; 
   double ** sigma = Allocate2DMatrix(9, 9); 
   double * theta = new double[10]; 
   double ** count9 = Allocate2DMatrix(nb,10); 
   //each allele frequency bin has 9 counts of different joint genotypes. 
   //each allele frequency bin has 1 count of number of SNPs. 
   double ** count9_all = Allocate2DMatrix(nb, 10); 
   for (int i = 0; i < nb; i++) 
       for (int j = 0; j < 10; j++) 
       	   count9_all[i][j] = 0; 

   int np = 9; 
   double ** AtA = Allocate2DMatrix(9,9); 
   double * AtB = new double[9]; 
   double w[9];
   double zz[9];
   int index[9]; 
   double r2; 
   double x[9]; 

//   clock_t clk0 = clock(); 
//	   count; 
   unsigned long ki = (long) (ns) * ii; 
   unsigned long kj = (long) (ns) * jj; 
   unsigned long qi = ki / 16;   
   uint8_t ri = ki % 16; 
   unsigned long qj = kj / 16; 
   uint8_t rj = kj % 16; 
   double total_phi = 0; 
   double total_weight = 0; 
   for (int chr = 0; chr < 22; chr++)
   {
       int beg = gdadu.beg[chr];
       int end = gdadu.end[chr]; 
       double weight = (double) (end-beg)/1000; 
       for (int n = 0; n < nb; n++) 
	   for (int j = 0; j < 10; j++) 
	       count9[n][j] = 0; 

       for (int m = beg; m < end; m++)
       {
    //       uint8_t gi = (gdadu.gt[ki/4] >> ((ki%4)*2)) & 3; 
    //       uint8_t gj = (gdadu.gt[kj/4] >> ((kj%4)*2)) & 3; 
	   uint8_t gi = (gdadu.gt[qi] >> ((ri)*2)) & 3; 
	   uint8_t gj = (gdadu.gt[qj] >> ((rj)*2)) & 3; 
	   ri++; if((ri& 0xF) == 0) {ri &= ~0xFF; qi++;} 
	   rj++; if((rj& 0xF) == 0) {rj &= ~0xFF; qj++;} 

	   if((gi^3) == 0 || (gj^3) == 0) continue;     
	   //11 is missing
	   uint8_t cc = ((gi << 2) | gj) - gi; 
	   // 00 00    0  -0
	   // 00 01    1  -0
	   // 00 10    2  -0
	   // 01 00    4  -1
	   // 01 01    5  -1
	   // 01 10    6  -1
	   // 10 00    8  -2
	   // 10 01    9  -2
	   // 10 10    10 -2

	   int tn = (int) gdadu.raf[m]; 
	   count9[tn][8-cc]++;  
	   count9[tn][9]++;  
       }
       for (int i = 0; i < nb; i++)
	   for (int j = 0; j < 10; j++)
	   count9_all[i][j] += count9[i][j]; 
	       
//       int notwin = 0; 
//       for (int n = 0; n < nb; n++) 
//       {
//	   int set[] = {1,2,3,5,6,7};
//	   for(int j = 0; j < 6; j++)
//	       notwin += count9[n*10 + set[j]]; 
//       }
       //count genoytpes that are not AA AA, AB AB and BB BB. 
       //if the count notwin === 0, then it's twin or self. 

       for (int i = 0; i < 9; i++)
	   for (int j = 0; j < 9; j++)
	       AtA[i][j] = 0; 
       for (int j = 0; j < 9; j++)
	   AtB[j] = 0; 

       for (int n = 0; n < nb; n++)
       {
	   for (int j = 0; j < 10; j++)
	   {
	       theta[j] = 0; 
	       for (int k = gdadu.ileft[n]; k <= gdadu.iright[n]; k++)
		   theta[j] += count9[k][j]; 
	   }
	   double af = 0; 
	   double cc[9] = {4,3,2,3,2,1,2,1,0}; 
	   if(theta[9] > 0) {
	       for (int j = 0; j < 8; j++)
		   af += theta[j] * cc[j]; 
	       af /= theta[9]; 
	       af /= 4; 
	   }
	   gdm(af, m1); 
	   //use SNPs in a bin to re-compute allele frequency. 
	   //this makes it robust to allele freq. mis-specification. 

	   for (int i = 0; i < 9; i ++)
	       for (int j = 0; j < 9; j++)
	       {
		   sigma[i][j] = 0; 
		   for (int k = 0; k < 9; k++)
		       sigma[i][j] += m1[k*9+i] * m1[k*9+j]; 
		   AtA[i][j] += sigma[i][j] * theta[9]; 
	       }

	   for (int j = 0; j < 9; j++)
	       for (int k = 0; k < 9; k++)
		    AtB[j] += m1[k*9+j] * theta[k]; 
       }

       if(gdadu.gsl == 1) 
	   gsl_solve(AtA, AtB, x, np); 
       else 
	   nnls(AtA, np, np, AtB, x, &r2, w, zz, index); 
       double phi=x[0]+ (x[2]+x[4]+x[6])/2+x[7]/4; 
       {
	   total_phi += phi * weight; 
	   total_weight += weight; 
       }
   }
   double phi1 = total_phi / total_weight; 

   if(gdadu.bychr == 2) {
       for (int i = 0; i < 9; i++)
	   for (int j = 0; j < 9; j++)
	       AtA[i][j] = 0; 
       for (int j = 0; j < 9; j++)
	   AtB[j] = 0; 

       for (int n = 0; n < nb; n++)
       {
	   for (int j = 0; j < 10; j++)
	   {
	       theta[j] = 0; 
	       for (int k = gdadu.ileft[n]; k <= gdadu.iright[n]; k++)
		   theta[j] += count9_all[k][j]; 
	   }
	   double af = 0; 
	   double cc[9] = {4,3,2,3,2,1,2,1,0}; 
	   if(theta[9] > 0) {
	       for (int j = 0; j < 8; j++)
		   af += theta[j] * cc[j]; 
	       af /= theta[9]; 
	       af /= 4; 
	   }
	   gdm(af, m1); 
	   //use SNPs in a bin to re-compute allele frequency. 
	   //this makes it robust to allele freq. mis-specification. 

	   for (int i = 0; i < 9; i ++)
	       for (int j = 0; j < 9; j++)
	       {
		   sigma[i][j] = 0; 
		   for (int k = 0; k < 9; k++)
		       sigma[i][j] += m1[k*9+i] * m1[k*9+j]; 
		   AtA[i][j] += sigma[i][j] * theta[9]; 
	       }

	   for (int j = 0; j < 9; j++)
	       for (int k = 0; k < 9; k++)
		    AtB[j] += m1[k*9+j] * theta[k]; 
       }

       if(gdadu.gsl == 1) 
	   gsl_solve(AtA, AtB, x, np); 
       else 
	   nnls(AtA, np, np, AtB, x, &r2, w, zz, index); 
       double phi2=x[0]+ (x[2]+x[4]+x[6])/2+x[7]/4; 
       phi1 = phi1 / 2 + phi2 / 2; 
   }

   gdadu.grm[tiktok] = phi1 * 2.0; 
   
   Free2DMatrix(AtA); 
   delete[] AtB; 
   delete[] m1; 
   Free2DMatrix(sigma); 
   delete[] theta; 
   Free2DMatrix(count9); 
   Free2DMatrix(count9_all); 
}

void jacquard3(void *par) 
{
    int nb = gdadu.nb; 
    int ns = gdadu.ns; 
    int ni = gdadu.ni; 

    long tiktok = (long) par; 
    long kk = tr2sq[tiktok]; 
    int ii = kk / ni; 
    int jj = kk % ni; 
   
   double * m1 = new double[81]; //we only use the last three columns
   double ** sigma = Allocate2DMatrix(nb*10, 3); 
   double * theta = new double[nb*10]; 
//   double ** sts = Allocate2DMatrix(100, 81); 
   double ** count9 = Allocate2DMatrix(nb,10); 

   int np = 3; 
   double ** AtA = Allocate2DMatrix(np,np); 
   double * AtB = new double[np]; 
   double w[3];
   double zz[3];
   int index[3]; 
   double r2; 
   double x[3]; 

//   clock_t clk0 = clock(); 
//	   count; 
   unsigned long ki = (long) ns * ii; 
   unsigned long kj = (long) ns * jj; 
   unsigned long qi = ki / 16;   
   uint8_t ri = ki % 16; 
   unsigned long qj = kj / 16; 
   uint8_t rj = kj % 16; 

   for (int n = 0; n < nb; n++) 
       for (int j = 0; j < 10; j++)
       	   count9[n][j] = 0; 

   for (int m = 0; m < ns; m++)
   {
//       uint8_t gi = (gdadu.gt[ki/4] >> ((ki%4)*2)) & 3; 
//       uint8_t gj = (gdadu.gt[kj/4] >> ((kj%4)*2)) & 3; 
       uint8_t gi = (gdadu.gt[qi] >> ((ri)*2)) & 3; 
       uint8_t gj = (gdadu.gt[qj] >> ((rj)*2)) & 3; 
       ri++; if((ri& 0xF) == 0) {ri &= ~0xFF; qi++;} 
       rj++; if((rj& 0xF) == 0) {rj &= ~0xFF; qj++;} 

       if((gi^3) == 0 || (gj^3) == 0) continue;     
       //11 is missing
       uint16_t tn = gdadu.raf[m]; 
       if(tn>nb) { 
	   tn = nb+nb-tn;    //turn af bin into maf bin. 
	   gi = 2 - gi;      //flip genotypes; 
	   gj = 2 - gj; 
       }
       uint8_t cc = ((gi << 2) | gj) - gi; 
       // 00 00    0  -0
       // 00 01    1  -0
       // 00 10    2  -0
       // 01 00    4  -1
       // 01 01    5  -1
       // 01 10    6  -1
       // 10 00    8  -2
       // 10 01    9  -2
       // 10 10    10 -2

       count9[tn][8-cc]++;  
       count9[tn][9]++;  
   }

//   clock_t clk1 = clock(); 
//	   build_matrix; 

   for (int n = 0; n < nb; n++)
   {
       for (int j = 0; j < 10; j++)
       {
	   theta[n*10+j] = 0; 
	   for (int k = gdadu.ileft[n]; k <= gdadu.iright[n]; k++)
	       theta[n*10+j] += count9[k][j]; 
       }
       double af = 0; 
       double cc[9] = {4,3,2,3,2,1,2,1,0}; 
       if(theta[n*10+9] > 0) {
	   for (int j = 0; j < 9; j++)
	       af += theta[n*10+j] * cc[j]; 
	   af /= theta[n*10+9]; 
	   af /= 4; 
       }
       gdm(af, m1); 

       for (int i = 0; i < 9; i ++)
	   for (int j = 0; j < np; j++)
	   {
	       sigma[n*10+i][j] = m1[i*9+6+j] * theta[n*10+9]; 
	   }   //here only use the last three columns of m1; 
       for (int j = 0; j < np; j++)
	   sigma[n*10+9][j] = theta[n*10+9]; 
   }

   for (int i = 0; i < np; i++)
       for (int j = 0; j < np; j++)
       {
	   AtA[i][j] = 0; 
	   for (int n = 0; n < nb; n++) 
	       for (int k = 0; k < 10; k++)
	       AtA[i][j] += sigma[n*10+k][i] * sigma[n*10+k][j]; // *wt[i] * wt[j]; 
       }
   for (int j = 0; j < np; j++)
   {
       AtB[j] = 0; 
       for (int n = 0; n < nb; n++) 
	   for (int k = 0; k < 10; k++)
	   	AtB[j] += sigma[n*10+k][j] * theta[n*10+k]; 
   }

   if(gdadu.gsl == 1) 
       gsl_solve(AtA, AtB, x, np); 
   else 
       nnls(AtA, np, np, AtB, x, &r2, w, zz, index); 
   double sumx = 0; 
   for (int i = 0; i < np; i++)
       sumx += x[i]; 
   static int warning = 0; 
   if(sumx > 1.02 || sumx < 0.98) 
   {
       if(warning == 0) {
	   fprintf(stdout, "sumd %8.6f deviates from 1, kindred renomralize \n. use -k and examine .kin.gz file to investigate further.\n", sumx);  
	   warning = 1; 
       }
   }
   for (int i = 0; i < np; i++)  //np=3
       x[i] /= sumx; 
   double phi=x[0]/2.0 + x[1]/4.0; 
   gdadu.grm[tiktok] = 2*phi; 

   if(gdadu.yesk == 1) 
   {
       int idx = tiktok % (gdadu.nth*1000); 
       gdadu.kin[idx][0] = phi; 
       gdadu.kin[idx][1] = sumx; 
       for (int i = 0; i < 6; i++)
	    gdadu.kin[idx][2+i] = 0;                         
       //the first 6 states has 0 weights.
       for (int i = 0; i < 3; i++)
	    gdadu.kin[idx][6+i] = x[i]; 
   }

   Free2DMatrix(AtA); 
   delete[] AtB; 
   delete[] m1; 
   Free2DMatrix(sigma); 
   delete[] theta; 
   Free2DMatrix(count9); 
}

//this convert grm matrix to gcta format. 
int make_grm(string fn, string pref, int lskip)
{
    gzFile fgz = gzopen(fn.c_str(), "rb"); 
    if(fgz == NULL) {
	    printf("can't open %s file to read. \n", fn.c_str()); 
	    exit(0); 
    }
    
   vector<double>  vv; // variable for input value
   string buf; 
    char c = gzgetc(fgz);  
    int countr = 0; 
    while (c != EOF)  {
	if(c == '\r' || c == '\n')
	    countr ++; 
	if(countr >= lskip )  
	    buf.push_back(c); 
	c=gzgetc(fgz); 
    }
    cout << buf.size() << endl; 
    cout << countr  << endl; 

    char * dup = strdup(buf.c_str()); 
    for (char * tok = strtok(dup, " \r\n\0"); tok; tok = strtok(NULL, " \r\n\0"))
    {
	double num = atof(tok); 
	vv.push_back(num); 
    }
    free(dup); 
    buf.assign(""); 
    
//    cout << vv.size() << endl; 
   int ni = sqrt(vv.size()); 
   cout << ni << endl; 
   if(ni*ni != vv.size()) 
   {
       cerr << "not a square matrix, check if header is there" << endl; 
       exit(0); 
   }

    string buf1(pref); 
    buf1.append(".grm.gz"); 
    gzFile fp1 = gzopen(buf1.c_str(), "wb"); 
    if(fp1 == NULL) 
    {
	    printf("can't open %s file to write\n", buf1.c_str()); 
	    exit(0); 
    }
   for (int i = 0; i < ni; i++)
       for (int j = 0; j <= i; j++)
       {
	   char buf[256]; 
	   sprintf(buf, "%d %d 100000 %f\n", i+1, j+1, vv.at(i*ni+j));  
	   gzprintf(fp1, "%s",  buf); 
       }
   gzclose(fp1); 

   string buf2(pref); 
   buf2.append(".grm.id"); 
   FILE * fp2 = fopen(buf2.c_str(), "w");
   if (fp2 == NULL) {
	    fprintf(stderr, "can't open file %s\n", buf2.c_str()); 
	    exit(EXIT_FAILURE);
   }    
   for (int i = 0; i < ni; i++)
       fprintf(fp2, "%d0000\t%d\n", i+1, i+1);  
   fclose(fp2); 

   return 0;
}


void gsl_solve(double ** AtA, double * AtB, double * x, int np)
{
    gsl_matrix * mU = gsl_matrix_alloc(np, np); 
    gsl_vector * vB = gsl_vector_alloc(np); 
    gsl_vector * vx = gsl_vector_alloc(np); 

    gsl_matrix * mV = gsl_matrix_alloc(np, np); 
    gsl_vector * vS = gsl_vector_alloc(np); 
    gsl_vector * work = gsl_vector_alloc(np); 

    for (int i = 0; i < np; i++)
	for (int j = 0; j < np; j++)
	    gsl_matrix_set(mU, i, j, AtA[i][j]); 
    for (int i = 0; i < np; i++)
	    gsl_vector_set(vB, i, AtB[i]); 

    gsl_linalg_SV_decomp(mU, mV, vS, work); 
    gsl_linalg_SV_solve(mU, mV, vS, vB, vx); 

    for (int i = 0; i < np; i++)
	x[i] = gsl_vector_get(vx, i); 

    gsl_matrix_free(mU); 
    gsl_vector_free(vB); 
    gsl_vector_free(vx); 

    gsl_matrix_free(mV); 
    gsl_vector_free(vS); 
    gsl_vector_free(work); 
}
