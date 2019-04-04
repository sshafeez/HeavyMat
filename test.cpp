#include "matrix.h"
#include <iostream>
using namespace std;

int main(){
    Matrix<int> A(4,5,true);
    Matrix<int> B(5,4,true);
    A.print();
    cout<<endl;
    B.print();
    cout<<endl;
    Matrix<int>* C = A*B;
    C->print();
    return 0;

}