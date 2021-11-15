#!/bin/bash

g++ create-test-cases.cpp -o tester.o
g++ greedy-2-opt.cpp -o tsp.o 

mkfifo fifo
for i in {1..1}
do
    ./tester.o 1 > fifo & 
    time ./tsp.o < fifo
    echo ""
done
rm fifo\


# do some stuff here


