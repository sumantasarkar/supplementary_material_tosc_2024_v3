#include "main.h"
#include <sstream>

#define DLU 4

string itos(int i) {stringstream s; s << i; return s.str(); }
GRBEnv env = GRBEnv();


void lambda_scalar(GRBModel& model, vector<GRBVar>& A, vector<GRBVar>& B, int N, int row, int beta, GRBVar P, int lambda) {

    GRBLinExpr temp;
    temp = 0 ;
    GRBVar V, M;
    V = model.addVar(0, N, 0, GRB_INTEGER, "V_"+itos(row)+"_"+itos(beta) + "_" + itos(lambda));
    M = model.addVar(0, 1, 0, GRB_BINARY, "M_"+itos(row)+"_"+itos(beta) + "_" + itos(lambda));

    model.addConstr(M == 1 - P);
    for(int i = 0; i < N; i++){
        temp += ((lambda>>i)&1) * A[i];
        temp += ((lambda>>i)&1) * B[i];
    }
    temp += M;
    model.addConstr(temp == 2*V);

}

int dlct_model(int threadNumber, int n);
int dlct_model(int threadNumber, int n) {
    int i, j;
    try {
	    env.set(GRB_IntParam_LogToConsole, 1);
		env.set(GRB_IntParam_Threads, threadNumber);

        // Gurobi search parameters


		env.set(GRB_IntParam_MIPFocus, 1);

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
        vector<vector<GRBVar>> X(n, vector<GRBVar>(n));
        vector<GRBVar> Y(n);
        vector<vector<GRBVar>> A(n, vector<GRBVar>(N));

        for(i = 0 ; i<n ; i++){
            for(j = 0 ; j<n; j++){
                X[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "X_"+itos(i)+"_"+itos(j));
            }
            Y[i] = model.addVar(0, n-1, 0, GRB_INTEGER, "Y_"+itos(i));
        }

        // variables for binary representation of Y
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

        // Output value
        for(i = 0 ; i<n ; i++){
            temp = 0;
            for(j = 0 ; j<n; j++){
                temp += j*X[i][j];
            }
            model.addConstr(Y[i] == temp) ;
        }


        // Binary representations of Y
        for(i = 0 ; i<n; i++){
            temp = 0;
            for(j = 0 ; j<N; j++){
                temp = temp + (1<<j)* A[i][j];
            }
            model.addConstr(temp == Y[i]);
        }

        // Adding DLCT constraints

        for(int row = 0 ; row < n-1; row++){

            for(int la = 0; la < n-1; la++){

                GRBVar W[n/2];
                for(i = 0 ; i<n/2 ; i++){
                    W[i] = model.addVar(0, 1, 0, GRB_BINARY, "W_"+ itos(row) + "_" + itos(la)+"_"+itos(i));
                }


                for(i = 0 ; i<n/2; i++){
                    vector<GRBVar> T1(N);
                    vector<GRBVar> T2(N);
                    for(int k = 0 ; k<N; k++){
                        T1[k] = A[INDEXS[row][i][0]][k];
                        T2[k] = A[INDEXS[row][i][1]][k];
                    }
                    lambda_scalar(model, T1, T2, N, row, i, W[i], la + 1);
                }

                temp  = 0 ;
                for(i = 0 ; i<n/2; i++){
                    temp += W[i];
                }
                model.addConstr(temp <= n/4 + DLU/2);
                model.addConstr(temp >= n/4 - DLU/2);


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

int dlct_sbox(int size_sbox, int threadNumber) {
    dlct_model(threadNumber, size_sbox);
    return 0;
 }


