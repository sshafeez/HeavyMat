#include <vector>
#include <utility>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <algorithm>
using namespace std;

#define epsilon 0.000001


///////////////////// NAIVE MATRIX //////////////////////////////

class matrix{
    public:
    vector<vector<double>> grid;
	
    //Create N rows of M columns
    matrix(int rows, int cols, bool initRandom=false){
        grid.resize(rows, vector<double>(cols,0));
        if(initRandom){
            for(int i=0; i<rows; ++i){
                for(int j=0; j<cols; ++j){
                    grid[i][j] = rand() % 1000;
                }
            }
        }
    }

	//copy constructor
	matrix(matrix& mat) {
		grid =  mat.grid;
	}

	//append col vector
	void append(vector<double>& col){
		if (grid.empty()) {
			grid.resize(col.size());
		}
		for(int i=0; i<grid.size(); ++i){
			grid[i].push_back(col[i]);
		}
	}

    //individual element access
    double& at(int row, int col){
        return grid[row][col];
    }
    
    //print contents
    void print(){
        for(int i=0; i<grid.size(); ++i){
            for(int j=0; j<grid[0].size(); ++j){
                cout<<grid[i][j]<<" ";
                //cout<<std::setprecision(2)<<std::fixed<<grid[i][j]<<" ";
            }
            cout<<"\n";
        }
        cout<<endl;
    }

    //get size in <rows,columns>
    pair<int,int> size(){
        return pair<int,int>(grid.size(),grid[0].size());
    }

	//invert the matrix
	void invert() {
		if (size().first != size().second) throw "not invertible";
		//initialize I matrix on the right
		vector<vector<double>> identity;
		identity.resize(size().first, vector<double>(size().second, 0));
		for (int i = 0; i < size().first; ++i) {
			identity[i][i] = 1;
		}
		//start rref procedure
		int lead = 0;
		int rows = size().first;
		int cols = size().second;

		for (int r = 0; r < rows; ++r) {
			//if (rows == 5 && cols == 5) print();
			if (lead >= cols) return;
			int i = r;
			//find the next pivot
			while (abs(grid[i][lead]) < epsilon) {
				++i;
				if (i == rows) {
					i = r;
					++lead;
					if (lead == cols) {
						return;
					}
				}
			}
			//cout << "pivoting @ " << i << " " << lead << endl;

			//swap row i and r to put leading nonzero up high
			for (int k = 0; k < cols; ++k) {
				double tmp = grid[i][k];
				grid[i][k] = grid[r][k];
				grid[r][k] = tmp;

				//apply to I as well
				tmp = identity[i][k];
				identity[i][k] = identity[r][k];
				identity[r][k] = tmp;
			}

			//divide it so that leading nonzero is leading 1
			if (abs(grid[r][lead]) > epsilon) {
				double divisor = grid[r][lead];
				for (int k = 0; k < cols; ++k) grid[r][k] /= divisor;
				//apply to I as well
				for (int k = 0; k < cols; ++k) identity[r][k] /= divisor;
			}

			for (i = 0; i < rows; ++i) {
				if (i != r) {
					double mlpr = grid[i][lead];
					for (int k = 0; k < cols; ++k) {
						grid[i][k] -= grid[r][k] * mlpr;
					}
					//apply to I as well
					for (int k = 0; k < cols; ++k) {
						identity[i][k] -= identity[r][k] * mlpr;
					}
				}
			}
			++lead;
		}
		grid = identity;
	}

	//transpose the matrix
	void transpose() {
		vector<vector<double>> transposed;
		transposed.resize(size().second, vector<double>(size().first, 0));
		for (int i = 0; i < size().first; ++i) {
			for (int j = 0; j < size().second; ++j) {
				transposed[j][i] = grid[i][j];
			}
		}
		grid = transposed;
	}


};


//basic multiplication

void multiply(matrix& left, matrix& right, matrix& dest){
    if(left.size().second != right.size().first) throw "inner dimension mismatch";
	dest.grid.resize(left.size().first, vector<double>(right.size().second,0));

    /// (mxn) * (pxq) = (mxq)
    int rows = left.size().first;
    int n = left.size().second;
    int cols = right.size().second;
    for(int i=0; i<rows; ++i){
        for(int j=0; j<cols; ++j){
            dest.at(i,j) = 0;
            for(int k=0; k<n; ++k){
				dest.at(i,j) += left.at(i,k) * right.at(k,j);
            }
        }
    }
}


///////////////////// MATRIX WITH LINEAR DEPENDENCIES CACHED //////////////////////////////

class heavy_matrix : public matrix{
    public:
    struct dep{
        int index;
        double scalar;
    };
    unordered_map<int,vector<dep>> rowDeps;
    unordered_map<int,vector<dep>> colDeps;

    
    
    heavy_matrix(int rows, int cols, bool initRandom = false): matrix(rows,cols,initRandom){}

	//find linear dependencies
    void cache(double error=epsilon){
        rowDeps.clear(); 
        colDeps.clear(); 

		for (int iter = 0; iter < 2; ++iter) {
			//get row deps
			matrix A(0, 0); //holds basis
			vector<int> rowsInBasis = { 0 };
			A.append(grid[0]);
			for (int i = 1; i < size().first; ++i) {
				//cout << "A: \n"; A.print();
				//compute least squares B vector
				matrix b(0, 0); b.append(grid[i]);
				//cout << "b: \n"; b.print();
				matrix At = A; At.transpose();
				//cout << "At: \n"; At.print();
				matrix AtA(0, 0); multiply(At, A, AtA);
				//cout << "AtA: \n"; AtA.print();
				AtA.invert();
				//cout << "AtA inverted: \n"; AtA.print();
				matrix AtAAt(0, 0); multiply(AtA, At, AtAAt);
				//cout << "AtAAt : \n"; AtAAt.print();
				matrix AtAAtb(0, 0); multiply(AtAAt, b, AtAAtb);
				//cout << "AtAAtb : \n"; AtAAtb.print();

				//compute recreated vector
				matrix proj(0, 0); multiply(A, AtAAtb, proj);
				//cout << "proj: \n"; proj.print();

				//compute max deviation
				double cost = 0;
				for (int k = 0; k < b.size().first; ++k) {
					cost = max(cost, abs(b.grid[k][0] - proj.grid[k][0]));
				}

				//add linear dependence or add to basis
				if (cost <= error) {
					vector<dep> deps;
					for (int k = 0; k < rowsInBasis.size(); ++k) {
						if(abs(AtAAtb.grid[k][0]) > epsilon) deps.push_back({ rowsInBasis[k], AtAAtb.grid[k][0] });
					}
					if(iter==0) rowDeps[i] = deps; 
					else colDeps[i] = deps;

				}
				else {
					A.append(grid[i]);
					rowsInBasis.push_back(i);
				}
			}
			transpose();
		}
    }

};

void multiply(heavy_matrix& left, heavy_matrix& right, heavy_matrix& dest){
    if(left.size().second != right.size().first) throw "inner dimension mismatch";

    /// (mxn) * (nxq) = (mxq)
    int rows = left.size().first;
    int n = left.size().second;
    int cols = right.size().second;
    for(int i=0; i<rows; ++i){
        for(int j=0; j<cols; ++j){
			dest.at(i,j) = 0;
            if(left.rowDeps.count(i)){
                for(heavy_matrix::dep& comb : left.rowDeps[i]){
					dest.at(i,j) += comb.scalar * dest.at(comb.index,j);
                }
            }
            else if(right.colDeps.count(j)){
                for(heavy_matrix::dep& comb : right.colDeps[i]){
					dest.at(i,j) += comb.scalar * dest.at(i,comb.index);
                }
            }
            else{
                for(int k=0; k<n; ++k){
					dest.at(i,j) += left.at(i,k) * right.at(k,j);
                }
            }
            
        }
    }
}