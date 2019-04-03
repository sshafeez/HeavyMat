#include "matrix.h"

int main(){
    Matrix<int> mat(4,5);
    cout<<mat.size().first<<" rows\n"<< mat.size().second<<" cols\n"<<endl;
    mat.at(1,1) = 9;
    mat.print();
    return 0;

}