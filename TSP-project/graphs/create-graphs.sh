#!/bin/bash

# Find files ending with .dot extension
for i in $(find . -name '*.dot' ); 
do
    neato -Tpng $i -o "${i%.dot}.png" -n # Create PNG from .dot file
done