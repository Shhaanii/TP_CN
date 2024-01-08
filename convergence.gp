set term png size 800,600
set xlabel "Iterations"
set ylabel "Residu relatif"

set output "img/conv_richardson_alpha.png"
set title "Convergence du residu - Méthode de Richardson Alpha"
plot "bin/RESVEC_ALPHA.dat" title "Richardson Alpha" with lines

set output "img/conv_jacobi.png"
set title "Convergence du residu - Méthode de Jacobi"
plot "bin/RESVEC_JAC.dat" title "Jacobi" with lines

set output "img/conv_gauss.png"
set title "Convergence du residu - Méthode de Gauss-Seidel"
plot "bin/RESVEC_GS.dat" title "Gauss-Seidel" with lines

# pour comparaison
set output "img/conv_jacobi_gs.png"
set title "Convergence du residu - Jacobi vs Gauss-Seidel"
plot "bin/RESVEC_JAC.dat" title "Jacobi" with lines, \
     "bin/RESVEC_GS.dat" title "Gauss-Seidel" with lines