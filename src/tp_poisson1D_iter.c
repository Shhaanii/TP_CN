/******************************************/
/* tp2_poisson1D_iter.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include <time.h>
#include "lib_poisson1D.h"

#define ALPHA 0
#define JAC 1
#define GS 2

int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, lab, kv;
  int *ipiv;
  int info;
  int NRHS;
  int IMPLEM = 0;
  double T0, T1;
  double *RHS, *SOL, *EX_SOL, *X;
  double *AB;
  double *MB;
  
  double temp, relres;
  struct timespec start, end;
  double elapsed_time;
  char *file_name;

  double opt_alpha;

  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  /* Size of the problem */
  NRHS=1;
  nbpoints=12;
  la=nbpoints-2;

  /* Dirichlet Boundary conditions */
  T0=5.0;
  T1=20.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  SOL=(double *) calloc(la, sizeof(double)); 
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  /* Setup the Poisson 1D problem */
  /* General Band Storage */
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=0;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;
  
  AB = (double *) malloc(sizeof(double)*lab*la);
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  
  /* uncomment the following to check matrix A */
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");
  
  /********************************************/
  /* Solution (Richardson with optimal alpha) */

  /* Computation of optimum alpha */
  opt_alpha = richardson_alpha_opt(&la);
  printf("Optimal alpha for simple Richardson iteration is : %lf\n",opt_alpha); 

  /* Solve */
  double tol=1e-3;
  int maxit=1000;
  double *resvec;
  int nbite=0;

  resvec=(double *) calloc(maxit, sizeof(double));

  /* Solve with Richardson alpha */
  if (IMPLEM == ALPHA) {
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    richardson_alpha(AB, RHS, SOL, &opt_alpha, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Temps d'exécution pour ridchardson alpha : %f secondes\n", elapsed_time);

    write_vec(resvec, &nbite, "RESVEC_ALPHA.dat");
  }

  /* Richardson General Tridiag */

  /* get MB (:=M, D for Jacobi, (D-E) for Gauss-seidel) */
  kv = 0; // attention j'ai mit 0 mais de base c'était 1
  ku = 1;
  kl = 1;
  MB = (double *) malloc(sizeof(double)*(lab)*la);

  if (IMPLEM == JAC)
  {
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    extract_MB_jacobi_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    file_name = "RESVEC_JAC.dat";
  } else if (IMPLEM == GS) 
  {
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    extract_MB_gauss_seidel_tridiag(AB, MB, &lab, &la, &ku, &kl, &kv);
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    file_name = "RESVEC_GS.dat";
  }
  /* Solve with General Richardson */
  if (IMPLEM == JAC || IMPLEM == GS) {

    write_GB_operator_colMajor_poisson1D(MB, &lab, &la, "MB.dat");
    richardson_MB(AB, RHS, SOL, MB, &lab, &la, &ku, &kl, &tol, &maxit, resvec, &nbite);

    elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    printf("Temps d'exécution pour MB avec %s: %f secondes\n", file_name, elapsed_time);

    write_vec(resvec, &nbite, file_name);
  }

  /* Write solution */
  write_vec(SOL, &la, "SOL.dat");

  /* Write convergence history */
  write_vec(resvec, &nbite, "RESVEC.dat");

  relres = relative_forward_error(RHS, EX_SOL, &la);

  printf("\nThe relative forward error is relres = %e\n", relres);

  free(resvec);
  free(RHS);
  free(SOL);
  free(EX_SOL);
  free(X);
  free(AB);
  free(MB);
  printf("\n\n--------- End -----------\n");
}
