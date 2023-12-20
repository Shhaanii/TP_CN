/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include <math.h>
#include <string.h>

#include "lib_poisson1D.h"

void eig_poisson1D(double *eigval, int *la)
{
  memset(eigval, 0, sizeof(double) * *la);
  for (int k = 1; k <= *la; ++k)
  {
    eigval[k - 1] = 2 - (2 * cos(k * M_PI / (*la + 1)));
  }
}

double eigmax_poisson1D(int *la)
{
  return (2 - (2 * cos(*la * M_PI / (*la + 1))));
}

double eigmin_poisson1D(int *la){
  return (2 - (2 * cos(1 * M_PI / (*la + 1))));
}

double richardson_alpha_opt(int *la){
  return (2 / (eigmax_poisson1D(la) + eigmin_poisson1D(la)));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite)
{
  // RHS : right hand side
  // tol : tolérance, epsilon
  // maxit : nombre maximal d'itération
  // resvec : stocke le résidu
  // nbites : nombre d'itération
  
  /*

  r = norm(b - Ax) / norm(b)
  -> calcule norme de b
  -> dgbmv sur b-Ax puis norme de nouveau b 
  RESVEC = RSH - AB * X
  while RESVEC[k+1] > tol && maxit > 0 
  Xk+1] = X[k] + alpha_rich(RHS - AB * X[k]) 
  RESVEC
  
  --maxit
  nbr itération ++
  */
 
  // || b ||
  const double norm_b = cblas_dnrm2(*la, RHS, 1);

  // b = b - Ax
  cblas_dgbmv(Cbla, CblasNoTrans, *la, *lab, *kl, *ku, -1.0, AB, *lab, X, 1.0, 1.0, RHS, 1.0);

  // calcul residu -> ||  b - Ax || / || b ||
  const double residu = cblas_dnrm2(*la, RHS, 1) / norm_b;




}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

