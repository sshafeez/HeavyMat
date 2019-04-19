#include <vector>
#include "image.h"
#include "matrix.h"

using namespace std;

int main() {
	heavy_matrix R(0, 0);
	heavy_matrix G(0, 0);
	heavy_matrix B(0, 0);
	
	image_read("Mario_Luigi.png", R.grid, G.grid, B.grid);
	R.cache(100); R.writeback(255);
	cout << "R \n";
	G.cache(100); G.writeback(255);
	cout << "G \n";
	B.cache(100); B.writeback(255);
	cout << "B \n";
	  
	image_write("Mario_Luigi2.png", R.grid, G.grid, B.grid);

    return 0;
}
