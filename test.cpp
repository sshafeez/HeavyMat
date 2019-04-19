#include <vector>
#include <thread>
#include "image.h"
#include "matrix.h"

using namespace std;

void process(heavy_matrix& mat){
	mat.cache(epsilon);
	mat.writeback();

}

int main() {
	heavy_matrix R(0, 0);
	heavy_matrix G(0, 0);
	heavy_matrix B(0, 0);

	image_read("docs.png", R.grid, G.grid, B.grid);
	thread t1(process,ref(R));
	thread t2(process,ref(G));
	thread t3(process,ref(B));
	cout<<"R: "<<t1.get_id()<< " G: "<<t2.get_id()<<" B: "<<t3.get_id()<<endl;
	t1.join(); t2.join(); t3.join();
	image_write("docs2.png", R.grid, G.grid, B.grid);

    return 0;
}
