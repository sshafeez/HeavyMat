#include "matrix.h"
#include <iostream>
using namespace std;

int main(){
    heavy_matrix A(5,3,true);
	A.grid[0] = { 1,2,3 };
	A.grid[1] = { 3,5,7 };
	A.grid[2] = { 5,9,14 };
	A.grid[3] = { 3,2,1 };
	A.grid[4] = { 1,1,1 };

	A.cache(0.001);

   
    

    return 0;

}