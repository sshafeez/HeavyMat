#include <vector>
#include <thread>
#include "PNG_Interface/image.h"
#include "matrix.h"

using namespace std;

void process(heavy_matrix& mat, double error){
	mat.cache(error);
	mat.writeback();

}

int main(int argc, char** argv) {
	heavy_matrix R(0, 0);
	heavy_matrix G(0, 0);
	heavy_matrix B(0, 0);
	

	image_read(argv[1], R.grid, G.grid, B.grid);
	printLock.lock();
	thread t1(process,ref(R),atof(argv[2]));
	thread t2(process,ref(G),atof(argv[2]));
	thread t3(process,ref(B),atof(argv[2]));
	cout<<"R: "<<t1.get_id()<< " G: "<<t2.get_id()<<" B: "<<t3.get_id()<<endl; printLock.unlock();
	t1.join(); t2.join(); t3.join();
	

	heavy_matrix R_orig(0, 0);
	heavy_matrix G_orig(0, 0);
	heavy_matrix B_orig(0, 0);
	image_read(argv[1], R_orig.grid, G_orig.grid, B_orig.grid);
	double totalError=0;
	totalError += calcError(R,R_orig);
	totalError += calcError(G,G_orig);
	totalError += calcError(B,B_orig);
	cout<<"total error: "<<totalError<<endl;

	string file = to_string(atof(argv[2]))+ ".png";
	image_write(file.c_str(), R.grid, G.grid, B.grid);

	
    return 0;
}
