#include "matrix.h"
#include <iostream>
using namespace std;

int main(){
    Matrix A(4,5,true);
    Matrix B(5,4,true);
    A.print();
    cout<<endl;
    B.print();
    cout<<endl;
    Matrix* C = multiply(A,B);
    C->print();
    return 0;

}