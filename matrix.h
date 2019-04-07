#include <vector>
#include <utility>
#include <stdlib.h>
#include <iostream>
using namespace std;

class Matrix{
    private:
    vector<vector<double>> grid;

    public:
    //Create N rows of M columns
    Matrix(int rows, int cols, bool initRandom=false){
        if(rows==0 || cols==0){
            throw "Matrix cannot have 0 rows or columns";
        }
        grid.resize(rows, vector<double>(cols,0));
        if(initRandom){
            for(int i=0; i<rows; ++i){
                for(int j=0; j<cols; ++j){
                    grid[i][j] = rand() % 1000;
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
            }
            cout<<"\n";
        }
    }

    //get size in <rows,columns>
    pair<int,int> size(){
        return pair<int,int>(grid.size(),grid[0].size());
    }
    friend Matrix* multiply(Matrix& left, Matrix& right);

};


//basic multiplication

Matrix* multiply(Matrix& left, Matrix& right){
    if(left.size().second != right.size().first) throw "inner dimension mismatch";

    /// (mxn) * (pxq) = (mxq)
    int rows = left.size().first;
    int n = left.size().second;
    int cols = right.size().second;
    Matrix* result = new Matrix(rows,cols);

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