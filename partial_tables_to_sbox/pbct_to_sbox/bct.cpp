#include "main.h"
#include "array.h"
#include <sstream>

string itos(int i) {stringstream s; s << i; return s.str(); }
GRBEnv env = GRBEnv();

void XOR_word_new(GRBModel& model, int x, int y, vector<GRBVar>& A, vector<GRBVar>& B, GRBVar Z, int N) {
    vector<GRBVar> E(N), F(N) ;
    for(int i = 0 ; i<N ; i++){
        E[i] = model.addVar(0, 1, 0, GRB_BINARY, "E_"+itos(x) + "_" + itos(y) + "_"+itos(i));
        F[i] = model.addVar(0, 1, 0, GRB_BINARY, "F_"+itos(x) + "_" + itos(y) + "_"+itos(i));
    }

    for(int i = 0 ; i< N; i++){
        model.addConstr(A[i] + B[i] + E[i] == 2*F[i]);
    }

    GRBLinExpr temp = 0 ;

    for(int i = 0 ; i< N; i++){
        temp = temp + (1<<i)* E[i];
    }

    model.addConstr(Z == temp);
}


void check_equality_b(GRBModel& model, int a, int b, int x, int y, GRBVar P, GRBVar Q, GRBVar S) {

    GRBVar T[8];
    for(int i = 0 ; i<8 ; i++){
        T[i] = model.addVar(0, 1, 0, GRB_BINARY, "T_"+ itos(a) + "_" + itos(b)+"_"+itos(x) + "_" + itos(y) + "_"+itos(i));
    }

    GRBLinExpr temp = 0;

    for(int i = 0 ; i<8; i++){
        temp += T[i];
    }
    model.addConstr(temp == 1- S);

    // S = 1 => P = b and Q = b
    model.addQConstr(S*(P-b) >= 0) ;
    model.addQConstr(S*(b-P) >= 0) ;
    model.addQConstr(S*(P-Q) >= 0) ;
    model.addQConstr(S*(Q-P) >= 0) ;

    // S = 0

    // P = b and Q >= b+1
    model.addQConstr(T[0]*(P-b) == 0);
    model.addQConstr(T[0]*(Q-b -1) >= 0);

    // P = b and Q <= b-1
    model.addQConstr(T[1]*(P-b) == 0);
    model.addQConstr(T[1]*(b -1 - Q) >= 0);

    // Q = b and P >= b+1
    model.addQConstr(T[2]*(Q-b) == 0);
    model.addQConstr(T[2]*(P-b -1) >= 0);

    // Q = b and P <= b-1
    model.addQConstr(T[3]*(Q-b) == 0);
    model.addQConstr(T[3]*(b -1 - P) >= 0);

    model.addQConstr(T[4]*(P-b-1) >= 0);
    model.addQConstr(T[4]*(Q-b-1) >= 0);

    model.addQConstr(T[5]*(P-b-1) >= 0);
    model.addQConstr(T[5]*(b -1 - Q) >= 0);

    model.addQConstr(T[6]*(b-1-P) >= 0);
    model.addQConstr(T[6]*(Q-b -1) >= 0);

    model.addQConstr(T[7]*(b-1-P) >= 0);
    model.addQConstr(T[7]*(b -1 - Q) >= 0);

}


int compute_inv(int A[][2], int lx, int ly, int n){
    int ind = 0;
    for(int i = 0 ; i<(n*(n-1))/2; i++){
        if((A[i][0] == lx && A[i][1]==ly ) || (A[i][0] == ly && A[i][1]==lx )){
            ind = i;
            break;
        }
    }
    return ind;
}

int bct_model(int threadNumber, int n, int rows);
int bct_model(int threadNumber, int n, int rows) {
    int i, j;
    try {
        env.set(GRB_IntParam_LogToConsole, 1);
        env.set(GRB_IntParam_Threads, threadNumber);

        // Gurobi search parameters
        env.set(GRB_IntParam_MIPFocus, 1);

        // Heuristics
        //env.set(GRB_DoubleParam_NoRelHeurTime, 3600*24*10);

        // All solutions
        //env.set(GRB_IntParam_PoolSearchMode, 2);
        //env.set(GRB_IntParam_PoolSolutions, 2000000);
        //env.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);

        int count = 0 ;
        int INDEXES [(n*(n-1))/2][2] ;
        for(int x = 0; x<n; x++){
            for(int y = x+1; y<n;y++){
                INDEXES[count][0] = x;
                INDEXES[count][1] = y;
                count++;
            }
        }

        int N = log2(n);
        GRBModel model = GRBModel(env);
        GRBVar X[n][n];
        GRBVar Y[n];
        GRBVar A[n][N];

        for(i = 0 ; i<n ; i++){
            for(j = 0 ; j<n; j++){
                X[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "X_"+itos(i)+"_"+itos(j));
            }
            Y[i] = model.addVar(0, n-1, 0, GRB_INTEGER, "Y_"+itos(i));
        }

        for( i = 0 ; i<n; i++){
            for(j = 0 ; j< N; j++){
                A[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "A_"+itos(i)+"_"+itos(j));
            }
        }

        GRBLinExpr temp;

        // Column sum
        for(j = 0 ; j < n ; j++){
            temp  = 0 ;
            for(i = 0 ; i<n; i++){
                temp += X[i][j];
            }
            model.addConstr(temp == 1) ;
        }

        // Row sum
        for(i = 0 ; i<n ; i++){
            temp  = 0 ;
            for(j = 0 ; j<n; j++){
                temp += X[i][j];
            }
            model.addConstr(temp == 1) ;
        }

        // Output values
        for(i = 0 ; i<n ; i++){
            temp = 0;
            for(j = 0 ; j<n; j++){
                temp += j*X[i][j];
            }
            model.addConstr(Y[i] == temp) ;
        }

        for(i = 0 ; i<n; i++){
            temp = 0;
            for(j = 0 ; j<N; j++){
                temp = temp + (1<<j)* A[i][j];
            }
            model.addConstr(temp == Y[i]);
        }


        // BCT[a, b] = | {(x, y) | S(x) + S(y) = b and S(x+a) + S(y+a) = b}|
        // Pairs (x, y); their XORs and integer representation.

        GRBVar Z[(n*(n-1))/2];
        for(i = 0 ; i<(n*(n-1))/2; i++){
            Z[i] = model.addVar(0, n-1, 0, GRB_INTEGER, "Z_"+itos(i));
        }
        vector<GRBVar> T1(N), T2(N) ;
        count = 0 ;
        for(int x = 0; x<n; x++){

            for(int k = 0 ; k<N; k++){
                T1[k] = A[x][k];
            }

            for(int y = x+1; y<n;y++){
                for(int k = 0 ; k<N; k++){
                    T2[k] = A[y][k];
                }
                XOR_word_new(model, x, y, T1, T2, Z[count], N);
                count++;
            }
        }


        int lx, ly, ind;
        for(int a = 1 ; a < rows; a++){
            for(int b = 1; b < n; b++){

                GRBVar S[(n * (n-1))/2];
                for(i = 0 ; i< (n * (n-1))/2 ; i++){
                    S[i] = model.addVar(0, 1, 0, GRB_BINARY, "S_"+ itos(a) + "_" + itos(b)+"_"+itos(i));
                }

                count = 0 ;
                for(int x = 0 ; x<n; x++){
                    for(int y = x+1; y<n; y++){
                        lx = x ^ a;
                        ly = y ^ a;
                        ind = compute_inv(INDEXES, lx, ly, n);
                        check_equality_b(model, a, b, x, y, Z[count], Z[ind], S[count]);
                        count++;
                    }
                }

                temp = 0 ;
                for(i = 0 ; i< (n * (n-1))/2 ; i++){
                    temp += S[i];
                }
                model.addConstr(temp == inv4_bct[a][b]/2);
                //model.addConstr(temp == inv4_bct[a][b]/2);
                //model.addConstr(temp == present_bct[a][b]/2);
                //model.addConstr(temp == gift_bct[a][b]/2);
                //model.addConstr(temp == skinny4_bct[a][b]/2);
                //model.addConstr(temp == inv5_bct[a][b]/2);
                //model.addConstr(temp == ascon_bct[a][b]/2);
                //model.addConstr(temp == keccak_bct[a][b]/2);

            }

        }

        model.optimize();

        // Print the solution

        /*
        int solCount = model.get(GRB_IntAttr_SolCount);
		cout << "Number of solutions:" << solCount<<endl;
		if (solCount >= 2000000000) {
		    cerr << "Number of solutions is too large" << endl;
		    exit(0);
		}

		for(i = 0; i<solCount; i++){
		    model.set(GRB_IntParam_SolutionNumber, i);
		    cout <<"[ ";
		    for(j = 0; j<n-1; j++){
		        cout<< round(Y[j].get(GRB_DoubleAttr_Xn)) <<", ";
		    }
		    cout<< round(Y[j].get(GRB_DoubleAttr_Xn)) <<"]";
		    cout << endl;
		}
        */

        cout << endl ;

        cout <<"S = SBox( ";
        for(i = 0 ; i<n-1; i++){

            cout<< round(Y[i].get(GRB_DoubleAttr_Xn)) <<", ";
        }
        cout<< round(Y[i].get(GRB_DoubleAttr_Xn)) <<")";

        cout << endl;

        if (model.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
            return -1;
        }
        else if ((model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)) {
            return 1 ;
        }
        else {
            return -2;
        }
    }
    catch (GRBException e) {
        cerr << "Error code = " << e.getErrorCode() << endl;
        cerr << e.getMessage() << endl;
    }

    catch (...) {
        cerr << "Exception during optimization" << endl;
    }
    return -1;
}

int bct_sbox(int size_sbox, int threadNumber, int rows) {
    bct_model(threadNumber, size_sbox, rows);
    return 0;
}

