#include <stdlib.h>
#include "integrators.h"


/* RK4 stepper */
void RK4_vstep(int n, double *Dx, double dt,
               void (*fv)(int n, double *f, double *x, double t, void *p),
               double *x, double t, void *par)
{
  double *k1 = calloc(n, sizeof(double));
  double *k2 = calloc(n, sizeof(double));
  double *k3 = calloc(n, sizeof(double));
  double *k4 = calloc(n, sizeof(double));
  double *w = calloc(n, sizeof(double));
  double *f = calloc(n, sizeof(double));
  int i;

  fv(n,f, x,t, par); /* set f = func(x,t) */
  for(i=0; i<n; i++) k1[i] = dt*f[i];

  for(i=0; i<n; i++) w[i] = x[i] + 0.5*k1[i]; 
  fv(n,f, w, t+0.5*dt, par);
  for(i=0; i<n; i++) k2[i] = dt*f[i];

  for(i=0; i<n; i++) w[i] = x[i] + 0.5*k2[i]; 
  fv(n,f, w, t+0.5*dt, par);
  for(i=0; i<n; i++) k3[i] = dt*f[i];

  for(i=0; i<n; i++) w[i] = x[i] + k3[i]; 
  fv(n,f, w, t+dt, par);
  for(i=0; i<n; i++) k4[i] = dt*f[i];

  for(i=0; i<n; i++) Dx[i] = (k1[i] + 2.0*(k2[i]+k3[i]) + k4[i])/6.0;

  free(k1);
  free(k2);
  free(k3);
  free(k4);
  free(w);
  free(f);
}

/* Simple Euler step */
void Euler_vstep(int n, double *Dx, double dt,
                 void (*fv)(int n, double *f, double *x, double t, void *p),
                 double *x, double t, void *par)
{
  double *f = calloc(n, sizeof(double));
  int i;

  fv(n,f, x,t, par); /* set f = func(x,t) */
  for(i=0; i<n; i++) Dx[i] = dt*f[i];

  free(f);
}
