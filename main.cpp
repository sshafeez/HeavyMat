#include <iostream>
#include "SVD.h"
#include "SVD.cpp"
#include "Matrix.h"
#include "image.h"
#include <vector>

using namespace std;


int main(int argc, const char * argv[]) {
    int rows = 250;
    int cols = 250;
    
    float* R = (float *) malloc(sizeof(float)*rows*cols);
    float* G = (float *) malloc(sizeof(float)*rows*cols);
    float* B = (float *) malloc(sizeof(float)*rows*cols);
    
    image_read("Mario_Luigi.png", rows, cols, R, G, B);
    
    Matrix<float> R1(rows,cols,R);
    Matrix<float> R2(rows,cols);
    Matrix<float> R3(rows,cols);
    
    Matrix<float> G1(rows,cols,G);
    Matrix<float> G2(rows,cols);
    Matrix<float> G3(rows,cols);
    
    Matrix<float> B1(rows,cols,B);
    Matrix<float> B2(rows,cols);
    Matrix<float> B3(rows,cols);
    
    SVD(R1,R2,R3);
    eliminateBasis(R2, 50); // Eliminate all sigma values < 50
    Matrix<float> resR = R1*R2*R3.Transposed();
    SVD(G1,G2,G3);
    eliminateBasis(G2, 50);
    Matrix<float> resG = G1*G2*G3.Transposed();
    SVD(B1,B2,B3);
    eliminateBasis(B2, 50);
    Matrix<float> resB = B1*B2*B3.Transposed();
    
    image_write("image.png", rows, cols, resR.GetRawData(), resG.GetRawData(), resB.GetRawData());
    
   /* SVD(G1);
    cout << "done\n";
    SVD(B1);
    cout << "done\n";*/
    //free(R); free(G); free(B);
    /*int rows = 5;
    int cols = 3;
    float* data = (float *) malloc(sizeof(float)*rows*cols);
    float data2[15] = {1,2,3,4,1,2,3,4,1,2,3,4,1,2,3};
    //vector<float> data2 = {1,2,3,4,1,2,3,4,1,2,3,4};
    for (int i = 0; i < rows*cols; ++i) {
        data[i] = data2[i];
    }
    Matrix<float> a(rows,cols,data);
    SVD(a);
    cout << '\n';
    a.Print();
    a.Print();
    int max = std::max(cols,rows);
    int min = std::min(cols,rows);
    Matrix<float> v(max,max);
    float* w;
    dsvd(a,rows,cols,w,v);
    cout << "u:\n";
    a.Print();
    cout << "v:\n";
    v.Print();
    cout << "vecs in Sigma:\n";
    for (int i = 0; i < min; ++i) {
        cout << w[i] << '\n';
    }*/
    return 0;
}
