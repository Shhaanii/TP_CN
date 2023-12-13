# Use an official base image -> Debian
FROM debian:latest

# Set the working directory
WORKDIR /app

# Install dependencies 
RUN apt-get update && \
    apt-get install -y make libblas-dev liblapack-dev libatlas-base-dev gcc gfortran apt-utils liblapacke-dev

# docker build -t conteneur_tp .
# docker run -it --rm -h=docker -v $PWD:/app conteneur_tp