#include <vector>
#include <utility>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <unordered_map>
using namespace std;


///////////////////// NAIVE MATRIX //////////////////////////////

class matrix{
    private:
    vector<vector<double>> grid;

    public:
    //Create N rows of M columns
    matrix(int rows, int cols, bool initRandom=false){
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
                cout<<std::setprecision(0)<<std::fixed<<grid[i][j]<<" ";
            }
            cout<<"\n";
        }
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
    public:
    struct dep{
        int index;
        int scalar;
    };
    unordered_map<int,vector<dep>> rowDeps;
    unordered_map<int,vector<dep>> colDeps;

    heavy_matrix(int rows, int cols, bool initRandom = false): matrix(rows,cols,initRandom){}
};

heavy_matrix* multiply(heavy_matrix& left, heavy_matrix& right){
    if(left.size().second != right.size().first) throw "inner dimension mismatch";

    /// (mxn) * (pxq) = (mxq)
    int rows = left.size().first;
    int n = left.size().second;
    int cols = right.size().second;
    heavy_matrix* result = new heavy_matrix(rows,cols);
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