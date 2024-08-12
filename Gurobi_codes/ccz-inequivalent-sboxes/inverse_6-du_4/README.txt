We explain how to run the codes for the verification of Table 4 in our paper.

## Prerequisites
- Gurobi C++ [version 9 or later]

## Run the code
Step 1:
- make && ./a.out -t 1 -n 64

Here t is the number of threads and n equals size of S-box

The current run write the outputs in files:
- SOL.sol_0.sol
- SOL.sol_1.sol
- .
- .

Step 2: Run python script.py which reads SOL.sol_*.sol files and returns a list of S-boxes.

Step 3: Test CCZ-inequivalence:
    - sage check_ccz-inequivalence.sage
