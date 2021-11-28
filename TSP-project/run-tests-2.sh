#!/bin/bash

g++ -std=c++0x create-test-cases.cpp -o tester.o &&
g++ -std=c++0x multiple-2-opt.cpp -o tsp.o &&

mkfifo fifo
for i in {1..1}
do
    ./tester.o 1 > fifo & 
    time ./tsp.o < fifo
done
rm fifo\


# do some stuff here
./graphs/create-graphs.sh

