#!/bin/bash

# Find files ending with .dot extension
for i in $(find . -name '*.dot' ); 
do
    # Using graphviz CLI tool
    neato -Tpng $i -o "${i%.dot}.png" -n # Create PNG from .dot file
done