We give an example of how to run the codes for reconstructing an S-box from given LAT.
The procedure is exactly same for other cryptographic tables. Thus, we omit those details.

## Prerequisites
- Gurobi C++ [version 9 or later]

## Run the code
- make && ./a.out -t 4 -n 32

Here t is the number of threads and n equals size of S-box

The current run takes the LAT of the inverse S-box (5-bit) and produces the following output:

S = SBox( 0, 1, 18, 28, 9, 23, 14, 12, 22, 4, 25, 16, 7, 15, 6, 13, 11, 24, 2, 29, 30, 26, 8, 5, 17, 10, 21, 31, 3, 19, 20, 27)

With SAGE command S.linear_approximation_table(), one can match the LAT of obtained S-box with the given LAT.

## Reconstructing S-boxes from other LATs
- In lat.cpp change the following:
    //model.addConstr(temp == chi_lat[la][lb] + n/2);
    //model.addConstr(temp == inv4_lat[la][lb] + n/2);
    //model.addConstr(temp == present_lat[la][lb] + n/2);
    //model.addConstr(temp == gift_lat[la][lb] + n/2);
    //model.addConstr(temp == skinny4_lat[la][lb] + n/2);
    model.addConstr(temp == inv5_lat[la][lb] + n/2);
    //model.addConstr(temp == ascon_lat[la][lb] + n/2);
    //model.addConstr(temp == keccak_lat[la][lb] + n/2);
    //model.addConstr(temp == inv6_lat[la][lb] + n/2);
    //model.addConstr(temp == apn6_lat[la][lb] + n/2);
- Change the value of n accordingly. For example, in case of inv6_lat, n is 64