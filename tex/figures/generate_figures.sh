#!/usr/bin/env bash
python3 ./../../src/plots.py
python3 ./../../src/model.py -s CPLEX -P ./random
mv ./random/*.pdf .

