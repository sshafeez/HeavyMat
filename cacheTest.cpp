#include <vector>
#include "image.h"
#include "matrix.h"

using namespace std;

int main() {
    heavy_matrix A(10,3,true);
	A.print(res);
    A.cache(50);
	A.writeback();
	cout << endl;
	A.print(res);
    return 0;
}
