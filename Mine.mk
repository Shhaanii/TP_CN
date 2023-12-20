#######################################
# Pour mon environnement
#######################################
CC=gcc
LIBSLOCAL=-L/usr/lib/ -L/usr/lib/x86_64-linux-gnu/ -llapack -lblas -lm
INCLUDEBLASLOCAL=-I/usr/include 
OPTCLOCAL=-fPIC -march=native