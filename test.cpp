#include <vector>
#include "image.h"

using namespace std;

int main() {
    vector<vector<double> > vec;
    image_read("usFlag.png", 425, 276, vec);
    image_write("image.png", vec);
    return 0;
}
