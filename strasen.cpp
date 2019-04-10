//
//  main.cpp
//  whatever
//
//  Created by Ozan Aktay on 3/30/19.
//  Copyright © 2019 Ozan Aktay. All rights reserved.
//

#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

class matrix {
public:
    matrix operator+(matrix mat) {
        matrix result(mat.size(),mat.size());
        for (int row = 0; row < mat.size(); ++row) {
            for (int col = 0; col < mat.size(); ++col) {
                result.grid[row][col] = grid[row][col] + mat.grid[row][col];
            }
        }
        return result;
    }
    int size() {
        return grid.size();
    }
    matrix(int rows, int cols, bool initRandom=false) {
        grid.resize(rows, vector<int>(cols,0));
    }
    void fill(vector<int> inputs) {
        int i = 0;
        for (int row = 0; row < grid.size(); ++row) {
            for (int col = 0; col < grid.size(); ++col) {
                grid[row][col] = inputs[i];
                ++i;
            }
        }
    }
    int& at(int row, int col){
        return grid[row][col];
    }
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
    vector<vector<int>> grid;
};

matrix join(matrix first, matrix second, matrix third, matrix fourth) {
    matrix result(first.size()*2, first.size()*2);
    int gridRow = 0;
    for (int row = 0; row < first.size(); ++row) {
        int gridCol = 0;
        for (int col = 0; col < first.size(); ++col) {
            result.grid[gridRow][gridCol] = first.grid[gridRow][gridCol];
            result.grid[gridRow][gridCol+first.size()] = second.grid[gridRow][gridCol];
            result.grid[gridRow+first.size()][gridCol] = third.grid[gridRow][gridCol];
            result.grid[gridRow+first.size()][gridCol+first.size()] = fourth.grid[gridRow][gridCol];
            ++gridCol;
        }
        ++gridRow;
    }
    return result;
}

matrix multiply(matrix left, matrix right) {
    //Base case
    if (left.size() == 2) {
        int a = left.grid[0][0];
        int b = left.grid[0][1];
        int c = left.grid[1][0];
        int d = left.grid[1][1];
        int e = right.grid[0][0];
        int f = right.grid[0][1];
        int g = right.grid[1][0];
        int h = right.grid[1][1];
        int m1 = a*(f-h);
        int m2 = (a+b)*h;
        int m3 = (c+d)*e;
        int m4 = d*(g-e);
        int m5 = (a+d)*(e+h);
        int m6 = (b-d)*(g+h);
        int m7 = (a-c)*(e+f);
        matrix result(2,2);
        result.grid[0][0] = m5+m4-m2+m6;
        result.grid[0][1] = m1+m2;
        result.grid[1][0] = m3+m4;
        result.grid[1][1] = m1+m5-m3-m7;
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
    matrix c11 = multiply(a11,b11)+multiply(a12,b21);
    matrix c12 = multiply(a11,b12)+multiply(a12,b22);
    matrix c21 = multiply(a21,b11)+multiply(a22,b21);
    matrix c22 = multiply(a21,b12)+multiply(a22,b22);
    return join(c11,c12,c21,c22);
}


int main(int argc, const char * argv[]) {
    matrix left(4,4);
    matrix right(4,4);
    left.fill({1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16});
    right.fill({17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32});
    multiply(left,right).print();
    return 0;
}
