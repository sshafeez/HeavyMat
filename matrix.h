#include <vector>
#include <utility>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <unordered_map>
#include <algorithm>
#include <math.h>  
#include <mutex>
using namespace std;

#define epsilon 0.00001
#define max_num 256

long int floating_mults = 0;
long int floating_adds = 0;

mutex printLock;

///////////////////// NAIVE MATRIX //////////////////////////////

class matrix{
    public:
    vector<vector<double>> grid;
	
    //Create N rows of M columns
	matrix(int rows, int cols, bool initRandom = false) {
		grid.resize(rows, vector<double>(cols, 0));
		if (initRandom) {
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < cols; ++j) {
					grid[i][j] = rand() % max_num;
				}
			}
		}
	}

	//copy constructors
	matrix(const matrix& mat) {
		grid =  mat.grid;
	}
	matrix(const vector<vector<double>>& data) {
		grid = data;
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
	void set(int row, int col, double val){
		grid[row][col]=val;
	}

    //get size in <rows,columns>
    pair<int,int> size(){
        return pair<int,int>(grid.size(),grid[0].size());
    }
	
    //matrix addition
    matrix operator+(matrix mat) {
        matrix result(mat.size().first,mat.size().first);
        for (int row = 0; row < mat.size().first; ++row) {
            for (int col = 0; col < mat.size().first; ++col) {
                result.grid[row][col] = grid[row][col] + mat.grid[row][col];
            }
        }
		floating_adds += mat.size().first * mat.size().first;
        return result;
    }
	matrix operator-(matrix mat) {
        matrix result(mat.size().first,mat.size().first);
        for (int row = 0; row < mat.size().first; ++row) {
            for (int col = 0; col < mat.size().first; ++col) {
                result.grid[row][col] = grid[row][col] - mat.grid[row][col];
            }
        }
		floating_adds += mat.size().first * mat.size().first;
        return result;
    }
	
    //return n/2 by n/2 matrix consisting of one of 4 quadrants of current matrix
    matrix split(int pos) {
        matrix mat = matrix(grid.size()/2,grid.size()/2);
        // if 0 or 1, start in row 0, end in row grid.size()/2-1;
        // if 2 or 3, start in row grid.size()/2, end in row grid.size();
        int rowStart = ((pos == 0 || pos == 1) ? 0 : grid.size()/2);
        int rowEnd = ((pos == 0 || pos == 1) ? grid.size()/2 : grid.size());
        // if 0 or 2, start in col 0, end in col grid.size()/2-1;
        // if 1 or 3, start in col grid.size()/2, end in col grid.size();
        int colStart = ((pos == 0 || pos == 2) ? 0 : grid.size()/2);
        int colEnd = ((pos == 0 || pos == 2) ? grid.size()/2 : grid.size());
        int matRow = 0;
        for (int row = rowStart; row < rowEnd; ++row) {
            int matCol = 0;
            for (int col = colStart; col < colEnd; ++col) {
                mat.grid[matRow][matCol] = grid[row][col];
                ++matCol;
            }
            ++matRow;
        }
        return mat;
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
	
	void saturate() {
        for (int i = 0; i < size().first; ++i) {
                for (int j = 0; j < size().second; ++j) {
                    if (grid[i][j] < 0) {
                            grid[i][j] = 0;
                    }	
                    if (grid[i][j] > 255) {
                            grid[i][j] = 255;
                    }
                }
        }
    }

	void print(){
		cout<<grid.size()<<" : "<<grid[0].size()<<endl;
		for(vector<double>& row : grid){
			for(double val : row){
				cout<<val<<" ";			
			}
			cout<<endl;
		}
	}

};

// Joins 4 matrices into one large matrix (2n by 2n) where each matrix is a quadrant of the larger matrix
matrix join(matrix &first, matrix &second, matrix &third, matrix &fourth) {
    matrix result(first.size().first*2, first.size().first*2);
    int gridRow = 0;
    for (int row = 0; row < first.size().first; ++row) {
        int gridCol = 0;
        for (int col = 0; col < first.size().first; ++col) {
            result.grid[gridRow][gridCol] = first.grid[gridRow][gridCol];
            result.grid[gridRow][gridCol+first.size().first] = second.grid[gridRow][gridCol];
            result.grid[gridRow+first.size().first][gridCol] = third.grid[gridRow][gridCol];
            result.grid[gridRow+first.size().first][gridCol+first.size().first] = fourth.grid[gridRow][gridCol];
            ++gridCol;
        }
        ++gridRow;
    }
    return result;
}


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
				floating_mults++; floating_adds++;
				dest.at(i,j) += left.at(i,k) * right.at(k,j);
            }
	    floating_adds--;
        }
    }
}

//Strassen multiplication

matrix strassen_multiply_helper(matrix &left, matrix &right) {
    //Base case
    if (left.size().first == 2) {
        int a = left.at(0,0);
        int b = left.at(0,1);
        int c = left.at(1,0);
        int d = left.at(1,1);
        int e = right.at(0,0);
        int f = right.at(0,1);
        int g = right.at(1,0);
        int h = right.at(1,1);
        int m1 = a*(f-h);
        int m2 = (a+b)*h;
        int m3 = (c+d)*e;
        int m4 = d*(g-e);
        int m5 = (a+d)*(e+h);
        int m6 = (b-d)*(g+h);
        int m7 = (a-c)*(e+f);
        matrix result(2,2);
        result.at(0,0) = m5+m4-m2+m6;
        result.at(0,1) = m1+m2;
        result.at(1,0) = m3+m4;
        result.at(1,1) = m1+m5-m3-m7;
	floating_mults += 7; floating_adds += 18;
        return result;
    }
    //Inductive step
    matrix a11 = left.split(0);
    matrix a12 = left.split(1);
    matrix a21 = left.split(2);
    matrix a22 = left.split(3);
    matrix b11 = right.split(0);
    matrix b12 = right.split(1);
    matrix b21 = right.split(2);
    matrix b22 = right.split(3);
    matrix c11 = strassen_multiply_helper(a11,b11)+strassen_multiply_helper(a12,b21);
    matrix c12 = strassen_multiply_helper(a11,b12)+strassen_multiply_helper(a12,b22);
    matrix c21 = strassen_multiply_helper(a21,b11)+strassen_multiply_helper(a22,b21);
    matrix c22 = strassen_multiply_helper(a21,b12)+strassen_multiply_helper(a22,b22);
    return join(c11,c12,c21,c22);
}

void strassen_multiply(matrix &left, matrix &right, matrix &dest) {
    dest = strassen_multiply_helper(left,right);
}


