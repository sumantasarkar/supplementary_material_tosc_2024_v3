We give an example of how to run the codes for finding an S-box with linearity at most L.
The procedure is exactly same for other cryptographic properties. Thus, we omit those details.

## Prerequisites
- Gurobi C++ [version 9 or later]

## Run the code
- make && ./a.out -t 1 -n 16

Here t is the number of threads and n equals size of S-box

The current run will produce the following output:

S = SBox( 12, 10, 7, 13, 2, 5, 0, 11, 6, 1, 8, 14, 3, 9, 4, 15)

With SAGE command S.linearity(), one can check that the linearity of this S-box is 8.

## Varying linearity and S-box sizes
- In lat.cpp change the value in macro "#define LINEARITY 8"
- Change n to 8, 16, 32 for 3, 4, 5 bit S-boxes, respectively.
