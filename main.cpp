#include <vector>
#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cmath>
#include "myrand.h"


#define NSUM 100000


std::vector<double> sum_errors(int nsamples);

int main(int argc,char** argv){

  auto errs = sum_errors(200);

  for(auto e : errs){
    std::cout<< std::scientific << std::setprecision(15) << e <<std::endl;
  }
  return 0;
}





std::vector<double> sum_errors(int nsamples){
  auto xs = rand_normal(0.0,20.0,NSUM);

  std::vector<float> xf(NSUM);
  std::vector<double> xd(NSUM);
  for(Integer i=0;i<NSUM;i++){
    xf[i]=float(xs[i]);
    xd[i]=double(float(xs[i]));
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


  std::vector<double> out;
  double relerr=0.0;
  /* Permutate M times */
  for (Integer j = 0; j < m; j++) {
    /* Set up the index vector */
    for (Integer i = 0; i < n; i++)
      index[i] = i;
    /* Call the permutation routine */
    nag_rand_permute(index, n, state, &fail);
    if (fail.code != NE_NOERROR) {
      printf("Error from nag_rand_permute (g05ncc).\n%s\n", fail.message);
    }

  float tmpf=0.0;
  double tmpd=0.0;
  for(Integer i=0;i<n;i++){
    tmpf += xf[index[i]];
    tmpd += xd[index[i]];
  }
  relerr = ((double(tmpf)-tmpd)/tmpd)/1e-6;
  out.push_back(relerr);

  }

  delete [] index;
  delete [] state;
  return out;





}
