#include <vector>
#include <utility>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <algorithm>
using namespace std;


///////////////////// NAIVE MATRIX //////////////////////////////

class matrix{
    protected:
    vector<vector<double>> grid;
	public:

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
			while (abs(grid[i][lead]) < 0.0001) {
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
			if (abs(grid[r][lead]) > 0.0001) {
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
		print();

	}

    //Create N rows of M columns
    matrix(int rows, int cols, bool initRandom=false){
        if(rows==0 || cols==0){
            throw "Matrix cannot have 0 rows or columns";
        }
        grid.resize(rows, vector<double>(cols,0));
        if(initRandom){
            for(int i=0; i<rows; ++i){
                for(int j=0; j<cols; ++j){
                    grid[i][j] = rand() % 100;
                }
            }
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
    friend matrix* multiply(matrix& left, matrix& right);

};


//basic multiplication

matrix* multiply(matrix& left, matrix& right){
    if(left.size().second != right.size().first) throw "inner dimension mismatch";

    /// (mxn) * (pxq) = (mxq)
    int rows = left.size().first;
    int n = left.size().second;
    int cols = right.size().second;
    matrix* result = new matrix(rows,cols);
    for(int i=0; i<rows; ++i){
        for(int j=0; j<cols; ++j){
            result->at(i,j) = 0;
            for(int k=0; k<n; ++k){
                result->at(i,j) += left.at(i,k) * right.at(k,j);
            }
        }
    }
    
    return result;
}


///////////////////// MATRIX WITH LINEAR DEPENDENCIES CACHED //////////////////////////////

class heavy_matrix : public matrix{
    private:
    struct dep{
        int index;
        double scalar;
    };
    unordered_map<int,vector<dep>> rowDeps;
    unordered_map<int,vector<dep>> colDeps;

    public:
    
    heavy_matrix(int rows, int cols, bool initRandom = false): matrix(rows,cols,initRandom){}
    void cache(){
        rowDeps.clear(); 
        colDeps.clear(); 
    }

    friend heavy_matrix* multiply(heavy_matrix& left, heavy_matrix& right);
};

heavy_matrix* multiply(heavy_matrix& left, heavy_matrix& right){
    if(left.size().second != right.size().first) throw "inner dimension mismatch";

    /// (mxn) * (nxq) = (mxq)
    int rows = left.size().first;
    int n = left.size().second;
    int cols = right.size().second;
    heavy_matrix* result = new heavy_matrix(rows,cols);
    for(int i=0; i<rows; ++i){
        for(int j=0; j<cols; ++j){
            result->at(i,j) = 0;
            if(left.rowDeps.count(i)){
                for(heavy_matrix::dep& comb : left.rowDeps[i]){
                    result->at(i,j) += comb.scalar * result->at(comb.index,j);
                }
            }
            else if(right.colDeps.count(j)){
                for(heavy_matrix::dep& comb : right.colDeps[i]){
                    result->at(i,j) += comb.scalar * result->at(i,comb.index);
                }
            }
            else{
                for(int k=0; k<n; ++k){
                    result->at(i,j) += left.at(i,k) * right.at(k,j);
                }
            }
            
        }
    }
    
    return result;
}