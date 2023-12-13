#######################################
# ambre.mk
# Default options for ambre computer
#######################################
CC=gcc
LIBSLOCAL=-L/ -L/usr/lib/x86_64-linux-gnu -llapack -lblas -lm
INCLUDEBLASLOCAL=-I/usr/include 
OPTCLOCAL=-fPIC -march=native