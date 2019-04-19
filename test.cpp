#include <vector>
#include <thread>
#include "image.h"
#include "matrix.h"

using namespace std;

void process(heavy_matrix& mat){
	mat.cache(51);
	mat.writeback();

}

int main() {
	heavy_matrix R(0, 0);
	heavy_matrix G(0, 0);
	heavy_matrix B(0, 0);

	image_read("M.png", R.grid, G.grid, B.grid);
	printLock.lock();
	thread t1(process,ref(R));
	thread t2(process,ref(G));
	thread t3(process,ref(B));
	cout<<"R: "<<t1.get_id()<< " G: "<<t2.get_id()<<" B: "<<t3.get_id()<<endl; printLock.unlock();
	t1.join(); t2.join(); t3.join();
	

	heavy_matrix R_orig(0, 0);
	heavy_matrix G_orig(0, 0);
	heavy_matrix B_orig(0, 0);
	image_read("M.png", R_orig.grid, G_orig.grid, B_orig.grid);
	double totalError=0;
	totalError += calcError(R,R_orig);
	totalError += calcError(G,G_orig);
	totalError += calcError(B,B_orig);
	cout<<"total error: "<<totalError<<endl;


	image_write("M_51.png", R.grid, G.grid, B.grid);

	
    return 0;
}
