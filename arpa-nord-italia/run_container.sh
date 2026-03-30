#!/bin/bash
cpp=/home/francesca16/Luca-Francesca/fdaPDE-cpp

docker run --rm -v $cpp:/root/fdaPDE-cpp \
    -v $(pwd):/root/arpa-nord-italia -ti aldoclemente/fdapde-docker:latest /bin/bash     