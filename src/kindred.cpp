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
//#include "gsl/gsl_linalg.h"
#include <pthread.h>
#include <stdint.h>
#include "thpool.h"


using namespace std; 

#define VERSION "0.3"

typedef struct KData {
    double * grm; 
    double ** kin; 
    uint8_t * gt; 
    double * raf; 
    int nb; 
    int ni; 
    int ns; 
    int nth; 
} KData; 

KData gdadu; 
map<long, long> tr2sq; 
map<long, long> sq2tr; 
//global data for multithreading. 

void jacquard(void *par);  

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
    double q = 1 - p; 
    double p2 = p * p; 
    double q2 = q * q; 
    double pq = p * q; 
    double p3 = p2 * p; 
    double q3 = q2 * q; 
    double p2q = p2 * q; 
    double pq2 = p * q2; 


    double tt[9][9] = { {p, p2, p2, p3, p2, p3, p2, p3, p3*p}, 
	  {0, 0,  pq, 2*p2q, 0, 0, 0, p2q, 2*p2*pq},
	  {0, pq, 0, pq2, 0, p2q, 0, 0, pq*pq}, 
	  {0, 0, 0, 0, pq, 2*p2q, 0, p2q, 2*p2*pq}, 
	  {0, 0, 0, 0, 0, 0, 2*pq, pq, 4*pq*pq}, 
	  {0, 0, 0, 0, pq, 2*pq2, 0, pq2, 2*q2*pq},
	  {0, pq, 0, p2q, 0, pq2, 0, 0, pq*pq}, 
	  {0, 0, pq, 2*pq2, 0, 0, 0, pq2, 2*q2*pq}, 
	  {q, q2, q2, q3, q2, q3, q2, q3, q3*q} }; 

    for (int i = 0, k=0; i < 9; i++)
	for (int j = 0; j < 9; j++,k++)
	    mm[k] = tt[i][j]; 
}

int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: kindred (kinship and inbreeding coefficients)\n");
	fprintf(stderr, "Version: %s\n", VERSION);
	fprintf(stderr, "Usage:   kindred [-i] in.vcf.gz [-abot]\n");
	fprintf(stderr, "Options: \n");
	fprintf(stderr, "         -a str        INFO/AF tag in vcf file [AF]\n");
	fprintf(stderr, "         -b int        number of allele freq bins [100]\n");
	fprintf(stderr, "         -i str        (indexed) bgzipped vcf file \n");
	fprintf(stderr, "         -o str        output prefix [out]\n");
	fprintf(stderr, "         -t int        number of thread [4]\n");
	fprintf(stderr, "Bug report: Yongtao Guan <ytguan@gmail.com>\n\n");
	return 1;
}


int main(int argc, char *argv[])
{
        time_t time0, time1, time2; 
	string fn_vcf; 
//	string fn_af; 
//	string fn_gt; 
	string pref("out");
	string af_tag("AF"); 
	char c;
        gdadu.nb = 100; //number of frequency bins; 
	gdadu.nth = 8; 

	while ((c = getopt(argc, argv, "a:b:i:o:t:")) >= 0) 
	{
	    switch (c) {
		case 'a': 
		    af_tag.assign(optarg); 
		    break; 
		case 'b': 
		    gdadu.nb = atoi(optarg); 
		    break;
		case 'i': 
		    fn_vcf.assign(optarg); 
		    break;
		case 'o': 
		    pref.assign(optarg); 
		    break;
		case 't': 
		    gdadu.nth = atoi(optarg); 
		    break;
		default: 
		    break; 
	    }
	}
    
//	std::cout << argc << " " << optind << std::endl; 
	if(argc > optind) 
	    fn_vcf.assign(argv[argc-1]); 
	if (argc < 2) return usage();

       time0 = time(NULL); 
       vector<float> vaf; 
        vector<uint8_t> v8; 
	//open vcf file; 
	htsFile *fpv = hts_open(fn_vcf.c_str(), "r");   
	bcf_hdr_t *hdr = bcf_hdr_read(fpv);  
	bcf1_t* line = bcf_init();   
	gdadu.ni = bcf_hdr_nsamples(hdr); 
	cout << "total number of samples = " << gdadu.ni << endl; 
//	if(hdr->samples != NULL) {
//	    for (int i = 0; i < gdadu.ni; i++)
//		cout << hdr->samples[i] << endl; 
//	}
	gdadu.ns = 0; 
	int ns0 = 0; 
	while(bcf_read1(fpv, hdr, line) == 0) {   
	  
	    bcf_get_variant_types(line); 
	    if(line->d.var_type != VCF_SNP) continue; 
	    if(line->n_allele > 2) continue;       

    	    float * dst = NULL;
	    int ndst = 0;
	    double af = 0; 
	    int dp_status = bcf_get_info_values(hdr, line, af_tag.c_str(),  (void**) &dst, &ndst, BCF_HT_REAL);
	    if(dp_status > 0) af = dst[0]; 
	    else ns0++; 
//	    cout << dst[0] << "\t " << af << endl; 
	    vaf.push_back(af); 
	    free(dst); 
			 
	    int32_t *gt_arr = NULL, ngt_arr = 0;
	    int ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
	    if ( ngt<=0 ) // GT not present 
	    {
		cout << "no genotypes at " << bcf_seqname(hdr, line) << " " << line->pos << endl; 
		continue; 
	    }
     
	    for (int i=0; i<gdadu.ni; i++)
	    {
		int gt = 3; 
		if (gt_arr[2*i+0] != bcf_gt_missing && gt_arr[2*i+1] != bcf_gt_missing)
		    gt = bcf_gt_allele(gt_arr[2*i+0]) + bcf_gt_allele(gt_arr[2*i+1]); 
		v8.push_back(gt); 
	   }
	   free(gt_arr);
	   gdadu.ns++; 
	   if(gdadu.ns % 10000 == 0) 
	       fprintf(stdout, "processed %d biallelic SNPs in vcf file \n", gdadu.ns); 
       }
       cout << "total number of SNPs = " << gdadu.ns << endl; 
       cout << "total number of SNPs without tag " << af_tag << " = " << ns0 << endl; 


   gdadu.raf = new double[gdadu.ns]; //rounded allele frequency. 
   for (int m = 0; m < gdadu.ns; m++)
       gdadu.raf[m] = vaf.at(m); 
   
   int blen = gdadu.ns*gdadu.ni/4+1; 
   gdadu.gt = new uint8_t[blen]; 
   for (int i = 0; i < blen; i++)
       gdadu.gt[i] &= ~0xFF; 
   for (int i = 0, k = 0; i < gdadu.ns; i++)
       for (int j = 0; j < gdadu.ni; j++, k++)
	   gdadu.gt[k/4] |= ((v8.at(k) & 3) << ((k%4)*2)); 
   vector<uint8_t> ().swap(v8); 

    gdadu.grm = new double[gdadu.ni*(gdadu.ni+1)/2];
    gdadu.kin = Allocate2DMatrix(gdadu.nth*1000,11); 

   for (int i=0, k=0; i<gdadu.ni; i++) 
       for (int j=i; j<gdadu.ni; j++, k++) 
       {
		   long tr = (long) k; 
		   long sq = (long) (i*gdadu.ni+j); 
		   tr2sq[tr] = sq; 
		   sq2tr[sq] = tr;
       }
   //global map for exchange indices; 

   time1 = time(NULL); 

//   double wt[9] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1, 1, 1}; 
    FILE * fp1;
    string buf; 
    buf.assign(pref); 
    buf.append(".kin");  
    fp1 = fopen(buf.c_str(), "w");
    if (fp1 == NULL) {
	    fprintf(stderr, "can't open file %s\n", buf.c_str()); 
	    exit(EXIT_FAILURE);
    }   // for SNPs. 

    FILE * fp2;
    buf.assign(pref); 
    buf.append(".grm");  
    fp2 = fopen(buf.c_str(), "w");
    if (fp2 == NULL) {
	    fprintf(stderr, "can't open file %s\n", buf.c_str()); 
	    exit(EXIT_FAILURE);
    }   // for SNPs. 

   //if i save upper triangular matrix linearly. 
   //how to efficienty recover coordinates from linear poistion?  


    fprintf(fp1, "ID1 ID2 phi sumd d1 d2 d3 d4 d5 d6 d7 d8 d9\n"); 

///////////////////////////////////////////////////////
    threadpool thpool = thpool_init(gdadu.nth);
///////////////////////////////////////////////////////

    int load = 0; 
    long tiktok = 0; 
    for (int ii = 0; ii < gdadu.ni; ii++)
       for (int jj = ii; jj < gdadu.ni; jj++)
       {
//	   long kk = ii * gdadu.ni + jj; 
///////////////////////////////////////////////////////
	   thpool_add_work(thpool, jacquard, (void*)tiktok);
///////////////////////////////////////////////////////
	   load++; 
	   tiktok++; 
	   if(load == gdadu.nth * 1000) {
///////////////////////////////////////////////////////
	       thpool_wait(thpool); 
///////////////////////////////////////////////////////
	       fprintf(stdout, "processed %ld many pair of samples \n", tiktok); 
	       for (long i = tiktok - load;  i < tiktok; i++)
	       {
		   long idx = tr2sq[i]; 
                   int i1 = idx / gdadu.ni; 
		   int i2 = idx % gdadu.ni; 
		   if(hdr->samples != NULL) 
			fprintf(fp1, "%s %s ", hdr->samples[i1], hdr->samples[i2]); 
		   else 
			fprintf(fp1, "id%d id%d ", i1, i2); 
		   for (int c = 0; c < 11; c++)
		       fprintf(fp1, "%6.5f ", gdadu.kin[i%(gdadu.nth*1000)][c]); 
		   fprintf(fp1, "\n"); 
	       }
	       load = 0; 
	   }
    //	   jacquard((void*) par); 
       }

///////////////////////////////////////////////////////
    thpool_wait(thpool);
    thpool_destroy(thpool);
///////////////////////////////////////////////////////
    if(load < gdadu.nth * 1000) {
       for (long i = tiktok - load;  i < tiktok; i++)
       {
	   long idx = tr2sq[i]; 
	   int i1 = idx / gdadu.ni; 
	   int i2 = idx % gdadu.ni; 
	   if(hdr->samples != NULL) 
		fprintf(fp1, "%s %s ", hdr->samples[i1], hdr->samples[i2]); 
	   else 
		fprintf(fp1, "id%d id%d ", i1, i2); 
	   for (int c = 0; c < 11; c++)
	       fprintf(fp1, "%6.5f ", gdadu.kin[i%(gdadu.nth*1000)][c]); 
	   fprintf(fp1, "\n"); 
       }
    }
//    for (int ii = 0, k = 0; ii < gdadu.ni; ii++)
//       for (int jj = ii; jj < gdadu.ni; jj++, k++)
//       {
//	   if(hdr->samples != NULL) 
//		fprintf(fp1, "%s %s ", hdr->samples[ii], hdr->samples[jj]); 
//	   else 
//		fprintf(fp1, "id%d id%d ", ii, jj); 
//	   for (int c = 0; c < 11; c++)
//	       fprintf(fp1, "%8.6f ", gdadu.kin[k][c]); 
//	   fprintf(fp1, "\n"); 
//       }
    fclose(fp1); 

    for (int i = 0; i < gdadu.ni; i++)
    {
	if(hdr->samples != NULL) 
		fprintf(fp2, "%s ", hdr->samples[i]); 
	else 
		fprintf(fp2, "id%d ", i); 
    }
    fprintf(fp2, "\n"); 
    for (int ii = 0; ii < gdadu.ni; ii++)
    {
       for (int jj = 0; jj < gdadu.ni; jj++)
       {
	   long idx = 0; 
	   if(ii > jj)  idx = sq2tr[(long) jj*gdadu.ni+ii]; 
	   else idx = sq2tr[(long) ii*gdadu.ni+jj]; 
	   fprintf(fp2, "%8.6f ", gdadu.grm[idx]); 
       }
       fprintf(fp2, "\n"); 
    }
    fclose(fp2); 

    time2 = time(NULL); 
    time_t dt1 = difftime(time1, time0); 
    time_t dt2 = difftime(time2, time1); 
    fprintf(stdout, "read vcf %ld seconds\n", dt1); 
    fprintf(stdout, "compute kinship %ld seconds\n", dt2); 

    Free2DMatrix(gdadu.kin); 
    delete[] gdadu.grm; 
    delete[] gdadu.gt; 
    delete[] gdadu.raf; 
    hts_close(fpv); 
    return 0;
}

void jacquard(void *par) 
{
    int nb = gdadu.nb; 
    int ns = gdadu.ns; 
    int ni = gdadu.ni; 

    long tiktok = (long) par; 
    long kk = tr2sq[tiktok]; 
    int ii = kk / ni; 
    int jj = kk % ni; 
   
   double * m1 = new double[81]; 
   double ** ste = Allocate2DMatrix(nb*10, 9); 
//   double ** sts = Allocate2DMatrix(100, 81); 
   double * count9 = new double[nb*10]; 
   double * count1 = new double[nb]; 
   //each allele frequency bin has 9 counts of different joint genotypes. 
   //each allele frequency bin has 1 count of number of SNPs. 

   int np = 9; 
   double ** AtA = Allocate2DMatrix(9,9); 
   double * AtB = new double[9]; 
   double w[9];
   double zz[9];
   int index[9]; 
   double r2; 
   double x[9]; 

//	   count; 
   for (int n = 0; n < nb; n++) 
   {
       count1[n] = 0; 
       for (int j = 0; j < 10; j++)
       count9[n*10+j] = 0; 
   }
   for (int m = 0; m < ns; m++)
   {
       if(gdadu.raf[m] < 0.01 || gdadu.raf[m] > 0.99) continue; 
       int ki = m*ni+ii; 
       int kj = m*ni+jj; 
       uint8_t gi = (gdadu.gt[ki/4] >> ((ki%4)*2)) & 3; 
       uint8_t gj = (gdadu.gt[kj/4] >> ((kj%4)*2)) & 3; 
       if(gi > 2 || gj > 2) continue; 
       int tgt = (2-gi) * 3 + 2-gj; 
//	       if(gt[m][ii] > 2 || gt[m][jj] > 2) //missing genotypes; 
//		   continue; 
//	       int tgt = (2-gt[m][ii]) * 3 + (2-gt[m][jj]); 
//               cout << (int) gi << " " << (int) gj << " " << tgt << endl; 
//	       raf[m] += (drand48()-0.5) * 0.05; 
//	       if(raf[m] < 0) raf[m] = 0; 
//	       if(raf[m] > 1) raf[m] = 0.99;
       int tn = floor(gdadu.raf[m] * nb); 
       count9[tn*10+tgt]+=1;  
       count1[tn]+=1;  
       int tn2 = nb-1-tn; 
       int tgt2 = 8-tgt; 
       count9[tn2*10+tgt2]+=1;  
       count1[tn2]+=1;  
   }
   
//	   build_matrix; 
   for (int n = 0; n < nb; n++) 
   {
       if(count1[n] == 0)  
	   continue; 
       double af = (n+0.5) / nb; 
       gdm(af, m1); 
       for (int i = 0; i < 9; i ++)
	   for (int j = 0; j < 9; j++)
	   {
	       ste[n*10+i][j] = m1[i*9+j] * count1[n]; 
	   }
       for (int j = 0; j < 9; j++)
	   ste[n*10+9][j] = count1[n]; 
       count9[n*10+9] = count1[n]; 
   }
   for (int i = 0; i < 9; i++)
       for (int j = 0; j < 9; j++)
       {
	   AtA[i][j] = 0; 
	   for (int n = 0; n < 10*nb; n++) 
	       AtA[i][j] += ste[n][i] * ste[n][j]; // *wt[i] * wt[j]; 
       }
   for (int j = 0; j < 9; j++)
   {
       AtB[j] = 0; 
       for (int n = 0; n < 10*nb; n++) 
	   AtB[j] += ste[n][j] * count9[n]; 
   }
   //solve; 

//	   for (int i = 0; i < np; i ++) 
//	   {
//	       for (int j = 0; j < np; j++)
//		   fprintf(stdout, "%5.4f ", AtA[i][j]); 
//	       fprintf(stdout, "%5.4f \n", AtB[i]); 
//	   }

   if(nnls(AtA, np, np, AtB, x, &r2, w, zz, index) != 0) 
   {
       cout << "nnls failed" << endl; 
       exit(0); 
   }  
   double phi=x[0]+ (x[2]+x[4]+x[6])/2+x[7]/4; 
   gdadu.grm[tiktok] = 2*phi; 
   double sumx = 0; 
   for (int i = 0; i < 9; i++)
       sumx += x[i]; 
   int idx = tiktok % (gdadu.nth*1000); 
   gdadu.kin[idx][0] = phi; 
   gdadu.kin[idx][1] = sumx; 
   for (int i = 0; i < 9; i++)
   	gdadu.kin[idx][2+i] = x[i]; 

   Free2DMatrix(AtA); 
   delete[] AtB; 
   delete[] m1; 
   Free2DMatrix(ste); 
   delete[] count9; 
   delete[] count1; 
}
