# Exercice numéro 3:

**1. En C, comment doit on déclarer et allouer une matrice pour utiliser BLAS et LAPACK ?**

Dans https://netlib.org/clapack/readme il est conseillé de : 

```c
double *A; 
A = malloc(M*N*sizeof(double));
````

**2. Quelle est la signification de la constante LAPACK COL MAJOR ?**

    LAPACK_COL_MAJOR est une constante qui permet de spécifier l'ordre de stockage des matrices dans LAPACK. Elle indique que les matrices sont stockées en mode "column-major". Ainsi les éléments d'une colonne sont stockés de manière consécutive en mémoire et les colonnes sont placées les une après les autres. 

https://www.netlib.org/lapack/lapacke


**3. A quoi correspond la dimension principale (leading dimension) généralement notée ld ?**

    Le "leading dimension" est la distance entre les éléments d'une colonne dans un tableau 2D.

**4. Que fait la fonction dgbmv ? Quelle méthode implémente-t-elle ?**

    DGBMV - banded matrix vector multiply (en précision double : 64 bits) Matrice tridiagonale surement et stockage en general Band

**5. Que fait la fonction dgbtrf ? Quelle méthode implémente-t-elle ?**

    DGBTRF computes an LU factorization of a real m-by-n band matrix A
    using partial pivoting with row interchanges.

    This is the blocked version of the algorithm, calling Level 3 BLAS.

    https://netlib.org/lapack//explore-html/df/dd6/group__gbtrf_ga725ad1dfee52c0dff571566d64e30ae8.html

**6. Que fait la fonction dgbtrs ? Quelle méthode implémente-t-elle ?**

    - DGBTRS solves a system of linear equations
    A * X = B  or  A**T * X = B
    with a general band matrix A using the LU factorization computed by DGBTRF.

    https://netlib.org/lapack//explore-html/d1/dd7/group__gbtrs_ga21ed607226c9690cadc991f15435fa74.html

**7. Que fait la fonction dgbsv ? Quelle méthode implémente-t-elle ?**

    DGBSV computes the solution to a real system of linear equations
    A * X = B, where A is a band matrix of order N with KL subdiagonals
    and KU superdiagonals, and X and B are N-by-NRHS matrices.

    The LU decomposition with partial pivoting and row interchanges is
    used to factor A as A = L * U, where L is a product of permutation
    and unit lower triangular matrices with KL subdiagonals, and U is
    upper triangular with KL+KU superdiagonals.  The factored form of A
    is then used to solve the system of equations A * X = B.

    https://www.netlib.org/lapack/explore-html/db/df8/group__gbsv_gaff55317eb3aed2278a85919a488fec07.html#gaff55317eb3aed2278a85919a488fec07 

**8. Comment calculer la norme du résidu relatif avec des appels BLAS**

    Il y a des fonctions qui calcule la normal en Blas : exemple en double : DNRM2 - Euclidean norm DZNRM2 - Euclidean norm
    r = norm(b - Ax) / norm(b)