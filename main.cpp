#include <vector>
#include <iostream>
#include <iomanip>

#include <cstdio>
#include <cmath>
#include <nag.h>
#include <nagg01.h>
#include <nagd01.h>

#include "myrand.h"


#define NSUM 100000

#define MAXMU (1.0)
#define MAXSIG (NSUM)
#define MINSIG (10.0)


typedef struct{
  double mu_mu;
  double mu_sigma;
  double sigma_mu;
  double sigma_sigma;
  std::vector<double>* xs;
}hypothesis_and_data_t;


extern "C" double likelihood(double mu,double sigma,Nag_Comm* dat);
extern "C" double prior(double mu,double sigma,Nag_Comm* dat);
extern "C" double unnormalized_posterior(double mu,double sigma,Nag_Comm* dat);
extern "C" double phi1(double y,Nag_Comm* dat){return -10.0;};
extern "C" double phi2(double y,Nag_Comm* dat){return  10.0;};



std::vector<double> sum_errors(int nsamples);

int main(int argc,char** argv){

  /*Generate errors*/
  auto errs = sum_errors(200);

  hypothesis_and_data_t d;
  d.mu_mu=0.0;
  d.mu_sigma=10.0;
  d.sigma_mu=10.0;
  d.sigma_sigma=NSUM;
  d.xs = &errs;
  Nag_Comm dat; 
  dat.p=(void*)(&d);


  NagError err;
  INIT_FAIL(err);
  double norm=0.0;
  double absacc=1e-14;
  Integer npts;
  nag_quad_2d_fin(0.001,20.0,phi1,phi2,unnormalized_posterior,absacc,&norm,&npts,&dat,&err);


  int musamples=1024;
  int sigmasamples=1024;
  std::vector<double> outimage;
  double map_m=-1000.0;
  double map_s=-1000.0;
  double max=-1000.0;
  for(int i=0;i<musamples;i++){
    for(int j=0;j<sigmasamples;j++){
      double m=-MAXMU + (double(i)/double(musamples))*(2*MAXMU);
      double s= MINSIG +  (double(j)/double(sigmasamples))*(MAXSIG-MINSIG);
      double post=unnormalized_posterior(m,s,&dat)/norm;
      outimage.push_back(
            MAX(post,1e-6)
          );
      if(MAX(post,max)==post){
        max=post;
        map_m=m;
        map_s=s;
      }
    }
  }
  std::cout<<"map_m = "<<map_m<<std::endl;
  std::cout<<"map_s = "<<map_s<<std::endl;






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

extern "C" double likelihood(double mu,double sigma,Nag_Comm* dat){
  hypothesis_and_data_t* _dat=(hypothesis_and_data_t*)(dat->p);
  std::vector<double>xs = *((std::vector<double>*) (_dat->xs));
  NagError err;
  INIT_FAIL(err);
  double out=1.0;
  for(auto x : xs){
    out=out*nag_normal_pdf(x,mu,sigma,&err);
  }
  return out;
}

extern "C" double prior(double mu,double sigma,Nag_Comm* dat){
  hypothesis_and_data_t* _dat=(hypothesis_and_data_t*)(dat->p);
  NagError err;
  INIT_FAIL(err);
  double Pmu=nag_normal_pdf(mu,_dat->mu_mu,_dat->mu_sigma,&err);
  double Psig=nag_normal_pdf(sigma,_dat->sigma_mu,_dat->sigma_sigma,&err);
  return Pmu*Psig;
}


extern "C" double unnormalized_posterior(double mu,double sigma,Nag_Comm* dat){
  return likelihood(mu,sigma,dat)*prior(mu,sigma,dat);
}
