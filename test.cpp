#include "matrix.h"
#include <iostream>
using namespace std;

int main(){
    heavy_matrix A(4,5,true);
    A.print();
    A.rref();
    A.print();

    heavy_matrix B(5,5,true);
    B.print();
    B.rref();
    B.print();

    heavy_matrix C(5,4,true);
    C.print();
    C.rref();
    C.print();

    return 0;

}