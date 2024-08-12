#include "main.h"
#include <sstream>

string itos(int i) {stringstream s; s << i; return s.str(); }
GRBEnv env = GRBEnv();

#define DU 6

// Constraints for word-wise XOR
void XOR_word(GRBModel& model, vector<GRBVar>& A, vector<GRBVar>& B, vector<GRBVar>& C, int N, int R, int S) {
    vector<GRBVar> D(N);
    for(int i = 0 ; i<N ; i++){
        D[i] = model.addVar(0, 1, 0, GRB_BINARY, "D_"+itos(R)+"_"+itos(S) + "_" + itos(i));
    }
    for(int i = 0 ; i< N; i++){
        model.addConstr(A[i] + B[i] + C[i] == 2*D[i]);
    }
}

int ddt_model(int threadNumber, int n);
int ddt_model(int threadNumber, int n) {
    int i, j;
    try {
	    env.set(GRB_IntParam_LogToConsole, 1);
		env.set(GRB_IntParam_Threads, threadNumber);

        // Gurobi search parameters
		env.set(GRB_IntParam_MIPFocus, 1);
		//env.set(GRB_IntParam_Presolve, 2);
		//env.set(GRB_IntParam_Method, 2);

        // Heuristics
		//env.set(GRB_DoubleParam_NoRelHeurTime, 3600*24*10);

        // Indices (solutions) for input a given input difference
        int INDEXS [n-1][n/2][2] ;
        for(int alpha = 0 ; alpha<n-1 ; alpha++){
            for(i = 0 ; i<n/2; i++){
                INDEXS[alpha][i][0] = 0 ;
                INDEXS[alpha][i][1] = 0 ;
            }
        }

        for(int alpha = 1 ; alpha<n ; alpha++){
            int count = 0 ;
            for(i = 0 ; i<n; i++){
                for(j = i+1 ; j<n; j++){
                    if((i ^ j) == alpha){
                        INDEXS[alpha-1][count][0] = i ;
                        INDEXS[alpha-1][count][1] = j ;
                        count++;
                    }
                }
            }
        }


        int N = log2(n);
        GRBModel model = GRBModel(env);
		GRBVar X[n][n];
        GRBVar Y[n];
        GRBVar A[n][N];
        GRBVar C[n-1][n/2][N];

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

        for(i = 0 ; i<n-1; i++){
            for(j = 0 ; j<n/2; j++){
                for(int k = 0 ; k<N; k++){
                    C[i][j][k] = model.addVar(0, 1, 0, GRB_BINARY, "C_"+itos(i)+"_"+itos(j) + "_" + itos(k));
                }
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

        // Output value
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

        for(i = 0 ; i< n-1 ; i++){
             for(j = 0 ; j<n/2; j++){
                 vector<GRBVar> T1(N);
                 vector<GRBVar> T2(N);
                 vector<GRBVar> T3(N);
                 for(int k = 0 ; k<N; k++){
                     T1[k] = A[INDEXS[i][j][0]][k];
                     T2[k] = A[INDEXS[i][j][1]][k];
                     T3[k] = C[i][j][k];
                 }
                 XOR_word(model, T1, T2, T3, N, i, j);
             }
        }

        GRBLinExpr temp1;

        for(int row = 0 ; row < n-1; row++){

            GRBVar S[n/2][n-1];
            for(i = 0 ; i<n/2 ; i++){
                for(j = 0 ; j<n-1; j++){
                    S[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S_"+ itos(row) + "_" + itos(i)+"_"+itos(j));
                }
            }

            // column sum
            for(j = 0 ; j < n-1 ; j++){
                temp  = 0 ;
                for(i = 0 ; i<n/2; i++){
                    temp += S[i][j];
                }
                model.addConstr(temp <= DU/2);

            }

            // Row sum

            for(i = 0 ; i<n/2 ; i++){
                temp  = 0 ;
                for(j = 0 ; j<n-1; j++){
                    temp += S[i][j];
                }
                model.addConstr(temp == 1) ;
            }


            for(i = 0 ; i<n/2; i++) {
                temp = 0;
                for (j = 0; j < n - 1; j++) {
                    temp += (j + 1) * S[i][j];
                }
                temp1 = 0;
                for (int k = 0; k < N; k++) {
                    temp1 = temp1 + (1 << k) * C[row][i][k];
                }
                model.addConstr(temp == temp1);

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

        //cout << solCount << endl ;
        //cout << endl ;

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

int ddt_sbox(int size_sbox, int threadNumber) {
    ddt_model(threadNumber, size_sbox);
    return 0;
 }

