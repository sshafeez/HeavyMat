#include "matrix.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

int main(){
	int width, height, bpp;

	uint8_t* rgb_image = stbi_load("image.png", &width, &height, &bpp, 0);
	cout << bpp;

	stbi_write_jpg("image2.jpg", width, height, bpp, rgb_image, 100);

	stbi_image_free(rgb_image);

	

   
    

    return 0;

}