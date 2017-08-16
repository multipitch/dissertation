#!/usr/bin/env bash

# Generate static plots
python3 ./../../src/plots.py

# Generate plots from random data 
./../../src/runmodel -s CPLEX -P ./random

# Move random data plots to this folder
mv ./random/*.pdf .

