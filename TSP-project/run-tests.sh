#!/bin/bash

g++ create-test-cases.cpp -o tester.o
g++ greedy-2-opt.cpp -o tsp.o 

mkfifo fifo
for i in {1..5}
do
    ./tester.o 1 > fifo & 
    ./tsp.o < fifo
    echo ""
done
rm fifo