#!/bin/bash

g++ -std=c++0x create-test-cases.cpp -o tester.o &&
g++ -std=c++0x greedy-2-opt.cpp -o tsp.o &&

mkfifo fifo
for i in {1..1}
do
    #./tester.o 1 > fifo & 
    #time ./tsp.o < fifo
    time ./tsp.o < ./test3.in
done
#rm fifo\


# do some stuff here


