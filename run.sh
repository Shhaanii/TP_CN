docker build -t conteneur_tp . && docker run -it --rm -h=$(hostname) -v $PWD:/app conteneur_tp
