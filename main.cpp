#include <vector>
#include <iostream>
#include <iomanip>


#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <nag.h>
#include <nagg01.h>
#include <nagd01.h>
#include <nagg05.h>
#include <omp.h>

#include "myrand.h"


#define NSUM 100000

#define MAXMU (1.0)
#define MAXSIG (NSUM)
#define MINSIG (1.0)




std::vector<double> sum_errors(int nsamples);

int main(int argc,char** argv){


  if(argc<=1){printf("Usage: ./main nsamples\n"); return 1;}
  Integer nsamples = atol(argv[1]);

  /*Generate errors*/
  auto errs = sum_errors(nsamples);


  Integer nb=nsamples;
  double xmean;
  double xsd;
  double xskew;
  double xkurt;
  double xmin;
  double xmax;
  Integer pn=0;
  NagError fail;
  INIT_FAIL(fail);
  nag_summary_stats_onevar(nb, &errs[0], NULL,&pn, &xmean, &xsd, &xskew, &xkurt, &xmin, &xmax, NULL, &fail);


  printf("mean = %lf, std = %lf\n",xmean,xsd);
  return 0;
}





std::vector<double> sum_errors(int nsamples){
  auto xs = rand_normal(0.0,20.0,NSUM);

  std::vector<double> out(nsamples);
  double* outd = &out[0];
  double* xsd = &xs[0];
#pragma omp parallel
  {

  int nthreads=omp_get_num_threads();
  int tid = omp_get_thread_num();
  int lnsamples = nsamples/(nthreads);
  int final_nsamples=nsamples%nthreads;

  std::vector<float> xf(NSUM);
  std::vector<double> xd(NSUM);
  for(Integer i=0;i<NSUM;i++){
    xf[i]=float(xsd[i]);
    xd[i]=double(float(xsd[i]));
  }


  /* Integer scalar and array declarations */
  Integer exit_status = 0;
  Integer  lstate;
  Integer *index = 0, *state = 0;

  /* NAG structures */
  NagError fail;

  /* Number of permutations */
  Integer m = nsamples;

  /* Sample size */
  Integer n = NSUM;

  /* Choose the base generator */
  Nag_BaseRNG genid = Nag_Basic;
  Integer subid = 0;

  /* Set the seed */
  Integer seed[] = { 1762543 };
  Integer lseed = 1;

  /* Initialize the error structure */
  INIT_FAIL(fail);


  /* Get the length of the state array */
  lstate = -1;
  nag_rand_init_repeatable(genid, subid, seed, lseed, state, &lstate, &fail);
  if (fail.code != NE_NOERROR) {
    printf("Error from nag_rand_init_repeatable (g05kfc).\n%s\n",
           fail.message);
  }

  /* Allocate arrays */
  index = new Integer[n];
  state = new Integer[lstate];

  /* Initialize the generator to a repeatable sequence */
  nag_rand_init_repeatable(genid, subid, seed, lseed, state, &lstate, &fail);
  if (fail.code != NE_NOERROR) {
    printf("Error from nag_rand_init_repeatable (g05kfc).\n%s\n",
           fail.message);
  }

/*  nag_rand_skip_ahead (NSUM*tid*lnsamples, state, &fail);
  if (fail.code != NE_NOERROR) {
    printf("Error from nag_rand_skip_ahead (g05kjc).\n%s\n", fail.message);
  }*/



  double relerr=0.0;
  /* Permutate M times */
  int beg=tid*lnsamples;
  int end=beg+lnsamples;
  if(tid==nthreads-1){
    end+=final_nsamples;
  }
/*
#pragma omp critical
  {
    printf("(thread,beg,end,nsamples) = (%d,%d,%d,%d)\n",tid,beg,end,nsamples);
  }
*/

  for (Integer j = 0; j < m; j++) {
    /* Set up the index vector */
    for (Integer i = 0; i < n; i++)
      index[i] = i;
    /* Call the permutation routine */
    nag_rand_permute(index, n, state, &fail);
    if (fail.code != NE_NOERROR) {
      printf("Error from nag_rand_permute (g05ncc).\n%s\n", fail.message);
    }

    if(j>=beg && j<end){
  float tmpf=0.0;
  double tmpd=0.0;
  for(Integer i=0;i<n;i++){
    tmpf += xf[index[i]];
    tmpd += xd[index[i]];
  }
  relerr = ((double(tmpf)-tmpd)/tmpd)/1e-6;
  outd[j]=relerr;
  }}

  delete [] index;
  delete [] state;
  }
  return out;

}


