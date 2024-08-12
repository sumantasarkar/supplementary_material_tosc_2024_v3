def sbox(file_name):
    S = []
    with open(file_name) as file:
        for L in file:
            if(L[0]=='Y'):
                ZZ = L.split()
                S.append(int(ZZ[-1]))
    return S

for i in range(1, 8):
    file_name = "SOL.sol_"+str(i)+".sol"
    S = sbox(file_name)
    print(S)