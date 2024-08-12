#include "main.h"
#include <sstream>

#define LINEARITY 8

string itos(int i) {stringstream s; s << i; return s.str(); }
GRBEnv env = GRBEnv();

void scalar(GRBModel& model, int bit, int LB[], vector<GRBVar>& A, GRBVar P, int x, int N, int la, int lb) {

    GRBLinExpr temp;
    GRBVar V, M;
    V = model.addVar(0, ceil(N+1/2), 0, GRB_INTEGER, "V_"+itos(la)+"_"+itos(lb) + "_" + itos(x));
    M = model.addVar(0, 1, 0, GRB_BINARY, "M_"+itos(la)+"_"+itos(lb) + "_" + itos(x));

    model.addConstr(M == 1 - P);
    temp = 0 ;
    for(int i = 0; i < N; i++){
        temp += LB[i] * A[i];
    }
    temp += bit;
    temp += M;
    model.addConstr(temp == 2*V);

}

void get_bin(int A[], int x, int N){
    for(int i = 0; i<N;i++){
        A[i] = ( x >> i ) & 1 ;
    }
}

int lat_model(int threadNumber, int n);
int lat_model(int threadNumber, int n) {
    int i, j;
    try {
	    env.set(GRB_IntParam_LogToConsole, 1);
		env.set(GRB_IntParam_Threads, threadNumber);

        // Gurobi search parameters

		env.set(GRB_IntParam_MIPFocus, 1);
        // Heuristics
		//env.set(GRB_DoubleParam_NoRelHeurTime, 3600*24*10);

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

        // Adding LAT constraints
        int LA[N], LB[N], LX[N], bit;

        for(int la = 1 ; la < n; la++){

            // Binary representation of la
            get_bin(LA, la, N);

            for(int lb = 1; lb < n; lb++){

                // Binary representation of lb
                get_bin(LB, lb, N);

                GRBVar W[n];
                for(i = 0 ; i<n ; i++){
                    W[i] = model.addVar(0, 1, 0, GRB_BINARY, "W_"+ itos(la) + "_" + itos(lb)+"_"+itos(i));
                }
                temp = 0;
                for(i = 0 ; i<n; i++){
                    temp += W[i];
                }

                model.addConstr(temp <= n/2 + LINEARITY/2);
                model.addConstr(temp >= n/2 - LINEARITY/2);


                for(int x = 0 ; x <n; x++){

                    // Binary representation of x
                    get_bin(LX,x,N);

                    bit = 0 ;
                    for(i = 0 ; i<N; i++){
                        bit ^= LA[i] & LX[i];
                    }

                    vector<GRBVar> T1(N);
                    for(int k = 0 ; k<N; k++){
                        T1[k] = A[x][k];
                    }
                    scalar(model, bit, LB, T1, W[x], x, N, la, lb);
                }


            }
        }


        model.optimize();

		// Print the solutions

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
        //cout << solCount<<endl;

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

int lat_sbox(int size_sbox, int threadNumber) {
    lat_model(threadNumber, size_sbox);
    return 0;
 }


