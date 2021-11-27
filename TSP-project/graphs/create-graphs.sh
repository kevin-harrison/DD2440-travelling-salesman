#!/bin/bash

# Find files ending with .dot extension
for i in $(find . -name '*.dot' ); 
do
    # Using graphviz CLI tool
    neato -Tpng $i -o "${i%.dot}.png" -n 2> /dev/null # Create PNG from .dot file
done