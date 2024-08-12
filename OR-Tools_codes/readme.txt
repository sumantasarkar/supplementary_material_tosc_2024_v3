Install OR-Tools from https://developers.google.com/optimization/install

## DDT to S-box ##
- Change the DDT in line 557 of ddt.py. The current input is apn6_ddt.
- On terminal type: python ddt.py 64 8
    -- 64 : 2^6
    -- 8  : number of threads

## LAT to S-box ##
- Change the LAT in line 460 of lat.py. The current input is keccak_lat.
- On terminal type: python lat.py 32 8
    -- 32 : 2^5
    -- 8  : number of threads

## DLCT/ACT to S-box ##
- Change the DLCT/ACT in line 338 of dlct_act.py. The current input is inv4_act.
- On terminal type: python dlct_act.py 32 8
    -- 16 : 2^4
    -- 8  : number of threads