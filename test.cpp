#include "matrix.h"
#include <iostream>
using namespace std;

int main(){
    matrix A(4,5,true);
    matrix B(5,4,true);
    A.print();
    cout<<endl;
    B.print();
    cout<<endl;
    matrix* C = multiply(A,B);
    C->print();
    return 0;

}