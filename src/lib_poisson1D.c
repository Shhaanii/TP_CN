/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <string.h>


void set_GB_operator_colMajor_poisson1D(double *AB, int *lab, int *la, int *kv) {
    // lab = nombre de colonnes
    // la = taille d'une ligne
    // kv => première ligne de la matrice bande pour la factorisation LU
    // Matrice v0 v1 .. u d l
    size_t taille = (*lab) * (*la);
    for (int indice = 0; indice < taille; indice += *lab) {

      for (int i = 0; i < *kv; ++i) {
        AB[indice + i] = 0.0;
      }
      // Assignation des valeurs U D L de la matrice poisson
      AB[indice + *kv] = -1.0;
      AB[indice + (*kv) + 1] = 2.0;
      AB[indice + (*kv) + 2] = -1.0;
    }
    // Remplace les extrémités de la matrice par 0.0
    AB[*kv] = 0.0;
    AB[taille - 1] = 0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double *AB, int *lab, int *la, int *kv)
{
  memset(AB, 0, sizeof(double) * *la * *lab);
  for (int indice = *kv + 1; indice < (*lab) * (*la); indice += *lab)
  {
    AB[indice] = 1.;
  }
}

void set_dense_RHS_DBC_1D(double *RHS, int *la, double *BC0, double *BC1)
{
  memset(RHS, 0, sizeof(double) * *la);
  RHS[0] = *BC0;
  RHS[*la - 1] = *BC1;
}


void set_analytical_solution_DBC_1D(double *EX_SOL, double *X, int *la, double *BC0, double *BC1)
{
  const double diff = *BC1 - *BC0;

  for (int i = 0; i < (*la); i++)
  {
    EX_SOL[i] = *BC0 + X[i] * diff;
  }

}

void set_grid_points_1D(double *x, int *la)
{
    // Calcul du pas entre les points sur la grille -> pour que les points soit équidistants
    const double step = 1.0 / ((*la) + 1); 

    for (int i = 0; i < *la; ++i)
    {
        x[i] = (i + 1) * step;
    }
}


double relative_forward_error(double* x, double* y, int* la){
  return 0;
}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  return *info;
}