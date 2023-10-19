#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
/* for random numbers */
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
gsl_rng *rnd;    /* global struct for random numbers */
#include "integrators.h"
#include "bins.h"

/*****************************************************************/
/* paramters used */
/*****************************************************************/

/* some values for ODE type we solve */
enum
{
  Andy2,        /* first version */
  with_v_s      /* same as Andy2 but moving at const vel. v_s */
};

/* global parameters and their default values */
long N = 100000000;       /* number of steps */
long N0 =10000000000; 
int ODE_type = Andy2;       /* ODE we solve  */
double dt = 0.005;         /* time step we use */
double tau1 = 20.0;
double tau2 = 50.0;
double Delta =1.0;
double Noise_b =1.0; 
double Omega2 = 0.01;
double gam = 0.0;           /* coeff of friction lineaar in velocity */
double v_s =2.0;          /* velocity with which we move through fluid */
int bins_num = 10000;         /* number of bins we use */
int bins_start = -500;       /* starting value for binning */
int bins_end   = 500;        /* ending value for binning */

/* set global pars to values in par file */
void read_parfile(char *file)
{
  char str[1000];
  FILE *fp;

  printf(" opening %s:\n", file);
  fp = fopen(file, "r");
  if(!fp)
  {
    printf("cannot open parfile!\n");
    exit(1);
  }

  while(fgets(str,999, fp))
  {
    char *s = strtok(str, "= ");
    //printf("s=|%s|\n", s);

    if(strcmp(s, "N")==0)
      N = atof(strtok(NULL, "= "));
    else if(strcmp(s, "N0")==0)
      N0 = atof(strtok(NULL, "= "));
    else if(strcmp(s, "dt")==0)
      dt = atof(strtok(NULL, "= "));
    else if(strcmp(s, "Omega2")==0)
      Omega2 = atof(strtok(NULL, "= "));
    else if(strcmp(s, "gam")==0)
      gam = atof(strtok(NULL, "= "));
    else if(strcmp(s, "Delta")==0)
      Delta = atof(strtok(NULL, "= "));
    else if(strcmp(s, "v_s")==0)
      v_s = atof(strtok(NULL, "= "));
    else if(strcmp(s, "Noise_b")==0)
      Noise_b = atof(strtok(NULL, "= "));
    //else if(strcmp(s, "bin_x")==0)
      //bin_x = atoi(strtok(NULL, "= "));
    else if(strcmp(s, "bins_num")==0)
      bins_num = atof(strtok(NULL, "= "));
    else if(strcmp(s, "bins_start")==0)
      bins_start = atof(strtok(NULL, "= "));
    else if(strcmp(s, "bins_end")==0)
      bins_end = atof(strtok(NULL, "= "));
    //printf("s=|%s|\n", s);
  }

  printf(" Parameters:\n");
  printf(" N = %ld\n", N);
  printf(" N0 = %ld\n", N0);
  printf(" dt = %g\n", dt);
  printf(" Omega2 = %g\n", Omega2);
  printf(" gam = %g\n", gam);
  printf(" Delta = %g\n", Delta);
  printf(" v_s = %g\n", v_s);
  printf(" Noise_b = %g\n", Noise_b);
  printf("\n");

  fclose(fp);
}
/*****************************************************************/
/* mean and std dev */
/*****************************************************************/
void mean_stddev(int dim, double *sum_x, double *sum_x2, long N,
                 double *mean, double *stddev)
{
  int k;
  for(k=0; k<dim; k++)
  {
    mean[k]   = sum_x[k]/N;
    stddev[k] = sqrt( (sum_x2[k]/N - mean[k]) );
  }
}
/*****************************************************************/
/* start of actual program */
/*****************************************************************/
/* sign function */
double sig(double v)
{
  if(v>0.0) return 1.0;
  if(v<0.0) return -1.0;
  return 0.0;
}

double DWiener_over_dt(double x, double dt)
{
  double g = gsl_ran_gaussian(rnd, sqrt(dt));
  return g/dt;
}
/* Andy's EOMs */
/* we are using Units 2 from Notes.txt */
void Andy2_EOM(int n, double *f, double *x, double t, void *p)
{
  double *par = p;
  double dt = par[11];
  double v = x[1]; /* x = x[0] */
  double dWdt = DWiener_over_dt(x[0], dt);
  f[0] = v;
  f[1] = -sig(v+v_s) - x[0] +(0.1)*sqrt(2.0)*dWdt;
}

/* main program */
int main(int argc, char *argv[])
{
  FILE *fpbinx;
  FILE *fpbinv;
  FILE *fpbinW1;
  FILE *fpbinW2;

  int i, k, n=2;
  long it, its_inbinsx,its_inbinsv,its_inbinsWi,its_inbinsW1,its_inbinsW2;
  int noutput;  /* number of steps we output, could be less */
  int N_ini;    /* number of initial steps that we all output */
  int outit;
  double t;
  long int Ntau1,Ntau2;
  Ntau1 = tau1/dt;
  Ntau2 = tau2/dt;



  double x[2], dx[2];
  double par[20];
  tBins *binsx;
  tBins *binsv;
  tBins *binsW1;
  tBins *binsW2;


  double *W1,*W2;
  W1 = calloc((Ntau1), sizeof(double));
  W2 = calloc((Ntau2), sizeof(double));

  /*double U1=0.0, U2=0.0, DeltaU;
  double K1=0.0, K2=0.0, DeltaK;
  double DeltaS;*/
  double  sum_x[] = {0., 0.};  /* sum of all x[0..1] */
  double sum_x2[] = {0., 0.};  /* sum of all x[0..1]^2 */
  double mean[2], stddev[2];
  double binsize  = (bins_end - bins_start)/(1.0*bins_num);
  printf("bin size is %g\n",binsize);
  /* init rnd */
  rnd = gsl_rng_alloc(gsl_rng_mt19937); /* use Mersenne Twister */
  gsl_rng_set(rnd, 11);                 /* set seed */

  /* check if we started program as: andy1 parfile */
  if(argc!=2)
  {
    printf("Usage: andy1 parfile\n");
    exit(1);
    printf(" argv[1]=%s -> N=%ld\n", argv[1], N);
  }
  /* read pars */
  read_parfile(argv[1]);
  noutput = N/1000000;

    /* init bins */
  binsx = alloc_bins(bins_num, bins_start, bins_end);
  binsv = alloc_bins(bins_num, bins_start, bins_end);
  binsW1 = alloc_bins(bins_num, bins_start, bins_end);
  binsW2 = alloc_bins(bins_num, bins_start, bins_end);
  /* calculate t0 and tf */
  //t0 = 0.0;
  //tf = N*dt;

  /* set some pars */
  outit = N/noutput;
  par[10] = N;
  par[11] = dt;
  printf("\n");
  printf(" outit = %d\n", outit);
  printf("\n");
  
  /* set and show initial x */
  t = 0.0;
  x[0] = -3.80; // 0.0; //1.0;
  x[1] = -0.656;
  /* set x,v to equilibrium without random noise force */
  //x[0] = -(v_s + Delta*sig(v_s))/Omega2;
  //x[1] = 0.;
  printf(" t            x(t)            v(t)\n");
  printf(" %e, %e, %e \n", t, x[0], x[1]);
  /* output some timesteps */
  N_ini = 20./dt;
  for(it=0; it<N_ini; it++)
  {
    //RK4_vstep(n, dx, dt, Andy2_EOM, x,t, par);
    Euler_vstep(n, dx, dt, Andy2_EOM, x,t, par);
    
    for(i=0; i<n; i++) x[i] += dx[i];
    t += dt;
    if(it%outit == 0)  printf(" %e, %e, %e\n", t, x[0], x[1]);
  }
  /* take many more steps to settle down in steady state */
  for(it=N_ini; it<N0; it++)
  {
    //RK4_vstep(n, dx, dt, Andy2_EOM, x,t, par);
    Euler_vstep(n, dx, dt, Andy2_EOM, x,t, par);
    for(i=0; i<n; i++) x[i] += dx[i];
    t += dt;
    if(it%outit == 0) printf(" %e, %e, %e\n", t, x[0], x[1]);
  }

  /* take N more steps */
  for(it=0; it<N; it++)
  {
    //RK4_vstep(n, dx, dt, Andy2_EOM, x,t, par);
    Euler_vstep(n, dx, dt, Andy2_EOM, x,t, par);
    #pragma omp parallel for
    for(k=0; k<(Ntau1);k++)
    W1[k] += -(v_s)*dt*x[0];
    #pragma omp parallel for
    for(k=0; k<(Ntau2);k++)
    W2[k] += -(v_s)*dt*x[0];
    	
     for(i=0; i<n; i++) 
     { 
      x[i] += dx[i];
     }
     #pragma omp parallel for
     for(k=0; k<(Ntau1);k++)
     {
       if(it%Ntau1 == k)
       {
         addto_bins(binsW1, W1[k]/tau1);
         W1[k]=0;
       }
     }
     #pragma omp parallel for
     for(k=0; k<(Ntau2);k++)
     {
       if(it%Ntau2 == k)
       {
         addto_bins(binsW2, W2[k]/tau2);
         W2[k]=0;
       }
     }

    t += dt;
    if(it%outit == 0) printf(" %e, %e, %e\n", t, x[0], x[1]);
      addto_bins(binsx, x[0]);
      addto_bins(binsv, x[1]);

    /* gather stats */
    for(i=0; i<n; i++)
    {
      sum_x[i]  += x[i];
      sum_x2[i] += x[i]*x[i];
    }
  }
  fpbinx=fopen("x.csv","w+");
  /* just output how many were in bins */
  for(i=0;i<num_of_bins(binsx);i++)
  { if(num_in_bin(binsx, i) !=0)
    fprintf(fpbinx,"%e  ,   %ld             ,  %e \n",
           pos_in_bins(binsx, i), num_in_bin(binsx, i),num_in_bin(binsx, i)/(1.*N*binsize)  );
  }
  fpbinv=fopen("v.csv","w+");
  /* just output how many were in bins */
  for(i=0;i<num_of_bins(binsv);i++)
  {
    if(num_in_bin(binsv, i) !=0)
    fprintf(fpbinv,"%e  ,   %ld             ,%e \n  ",
            pos_in_bins(binsv, i), num_in_bin(binsv, i), num_in_bin(binsv, i)/(1.*N*binsize) );
  }


  fpbinW1=fopen("Wtau1.csv","w+");
  /* just output how many were in bins */
  for(i=0;i<num_of_bins(binsW1);i++)
  {
    if(num_in_bin(binsW1, i) !=0)
    fprintf(fpbinW1,"%e  ,   %ld             ,%e \n  ",
            pos_in_bins(binsW1, i), num_in_bin(binsW1, i), num_in_bin(binsW1, i)/(1.) );
  }

  fpbinW2=fopen("Wtau2.csv","w+");
  /* just output how many were in bins */
  for(i=0;i<num_of_bins(binsW2);i++)
  {
    if(num_in_bin(binsW2, i) !=0)
    fprintf(fpbinW2,"%e  ,   %ld             ,%e \n  ",
            pos_in_bins(binsW2, i), num_in_bin(binsW2, i), num_in_bin(binsW2, i)/(1.) );
  }


  /* output bin average as well */
  printf("\n X bin average = %e\n\n", calc_bin_average(binsx));
  printf("\n V bin average = %e\n\n", calc_bin_average(binsv));
  /* print stats */
  its_inbinsx = total_in_bins(binsx);
  its_inbinsv = total_in_bins(binsv);
  its_inbinsW1 = total_in_bins(binsW1);
  its_inbinsW2 = total_in_bins(binsW2);

  mean_stddev(n, sum_x, sum_x2, it, mean, stddev);
  printf(" stats x :  %ld in bins of %ld iterations\n", its_inbinsx, it);
  printf(" stats v :  %ld in bins of %ld iterations\n", its_inbinsv, it);
  printf(" stats Wi :  %ld in bins of %ld iterations\n", its_inbinsWi, it);
  printf(" stats W1 :  %ld in bins of %ld iterations\n", its_inbinsW1, it);
  printf(" stats W2 :  %ld in bins of %ld iterations\n", its_inbinsW2, it);
  printf(" means:    ");
  for(i=0; i<n; i++) printf("%e  ", mean[i]);
  printf("\n");
  printf(" stddevs:  ");
  for(i=0; i<n; i++) printf("%e  ", stddev[i]);
  printf("\n");
  free(W1);
  free(W2);
  free_bins(binsx);
  free_bins(binsv);
  free_bins(binsW1);
  free_bins(binsW2);
  gsl_rng_free(rnd);
}