#include "matrix.h"
#include <iostream>
using namespace std;

int main(){
    heavy_matrix A(3,3,true);
    A.at(0,0) = 1; A.at(0,1) = 0; A.at(0,2) = 1;
    A.at(1,0) = 1; A.at(1,1) = 1; A.at(1,2) = 0;
    A.at(2,0) = 0; A.at(2,1) = 1; A.at(2,2) = 0;

	A.invert();

   
    

    return 0;

}