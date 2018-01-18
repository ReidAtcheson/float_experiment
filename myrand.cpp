#include <memory>
#include <iostream>
#include <cstdio>
#include <cstdint>

#include "myrand.h"
std::vector<double> rand_normal(double mu,double sigma,Integer n){
  Integer exit_status = 0;
  Integer i, lstate;
  Integer *state = 0;
  /* Double scalar and array declarations */
  double *x = 0;

  /* Set the distribution parameters */
  double xmu = mu;
  double var = sigma;


  /* Choose the base generator */
  Nag_BaseRNG genid = Nag_Basic;
  Integer subid = 0;

  /* Set the seed */
  Integer seed[] = { 1762543 };
  Integer lseed = 1;

  NagError fail;
  /* Initialize the error structure */
  INIT_FAIL(fail);
  lstate = -1;
  nag_rand_init_repeatable(genid, subid, seed, lseed, state, &lstate, &fail);
  if (fail.code != NE_NOERROR) {
    printf("Error from nag_rand_init_repeatable (g05kfc).\n%s\n",
             fail.message);
  }


  std::unique_ptr<double> xd(new double[n]);
  std::unique_ptr<Integer> stated(new Integer[lstate]);
  x=xd.get();
  state = stated.get();


  /* Initialize the generator to a repeatable sequence */
  nag_rand_init_repeatable(genid, subid, seed, lseed, state, &lstate, &fail);

  /* Generate the variates */
  nag_rand_normal(n, xmu, var, state, x, &fail);

  std::vector<double> out;
  for(Integer i=0;i<n;i++)out.push_back(x[i]);
  return out;
}



