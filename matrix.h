#include <vector>
#include <utility>
#include <stdlib.h>
#include <iostream>
using namespace std;

template <typename T>
class Matrix{
    private:
    vector<vector<T>> grid;

    public:
    //Create N rows of M columns
    Matrix(int rows, int cols, bool initRandom=false){
        if(rows==0 || cols==0){
            throw "Matrix cannot have 0 rows or columns";
        }
        grid.resize(rows, vector<T>(cols,0));
        if(initRandom){
            for(int i=0; i<rows; ++i){
                for(int j=0; j<cols; ++j){
                    grid[i][j] = rand() % 1000;
                }
            }
        }
    }

    //individual element access
    T& at(int row, int col){
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
    friend Matrix<T>* multiply(Matrix<T>& left, Matrix<T>& right);

};


//basic multiplication
template <typename T>
Matrix<T>* multiply(Matrix<T>& left, Matrix<T>& right){
    if(left.size().second != right.size().first) throw "inner dimension mismatch";

    /// (mxn) * (pxq) = (mxq)
    int rows = left.size().first;
    int n = left.size().second;
    int cols = right.size().second;
    Matrix<T>* result = new Matrix<T>(rows,cols);

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