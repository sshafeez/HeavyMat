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
                cout<<std::setprecision(2)<<std::fixed<<grid[i][j]<<" ";
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
    private:
    struct dep{
        int index;
        double scalar;
    };
    unordered_map<int,vector<dep>> rowDeps;
    unordered_map<int,vector<dep>> colDeps;

    public:
    //////RREF CODE//////
    
    void rref(){
        int lead = 0;
        int rows = size().first;
        int cols = size().second;

        for(int r=0; r<rows; ++r){
            if(lead >= cols) break;
            int i = r;
            while( abs(grid[i][lead]) < 0.0001){
                ++i;
                if(i == rows){
                    i=r;
                    ++lead;
                    if(lead == cols){
                        break;
                    }
                }
            }
            iter_swap(grid[i].begin(), grid[r].begin());
            if(abs(grid[r][lead]) > 0.0001){
                double divisor = grid[r][lead];
                for(int k=0; k<cols; ++k) grid[r][k] /= divisor;
            } 
            //for_each(grid[r].begin(),grid[r].end(), [divisor](double& val){ val /= divisor;});
            for(i=0; i<rows; ++i){
                if(i != r){
                    for(int k=0; k<cols; ++k) grid[i][k] -= grid[r][k] * grid[i][lead];
                }
            }
            ++lead;
        }

    }
    
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