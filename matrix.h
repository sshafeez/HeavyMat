#include <vector>
#include <utility>
#include <iostream>
using namespace std;

template <class T>
class Matrix{
    private:
    vector<vector<T>> grid;


    public:
    //Create N rows of M columns
    Matrix(int N, int M){
        if(M==0 || N==0){
            throw "Matrix cannot have 0 rows or columns";
        }
        grid.resize(N, vector<T>(M,0));
    }

    T& at(int row, int col){
        return grid[row][col];
    }
    
    Matrix* operator*(const Matrix<T> right){
        Matrix* result = new Matrix();
        return result;
    }
    
    void print(){
        for(int i=0; i<grid.size(); ++i){
            for(int j=0; j<grid[0].size(); ++j){
                cout<<grid[i][j]<<" ";
            }
            cout<<"\n";
        }
    }
    pair<int,int> size(){
        return pair<int,int>(grid.size(),grid[0].size());
    }

};

