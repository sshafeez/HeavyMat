#include <stdint.h>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image.h"
#include "stb_image_write.h"
#include <vector>

// Resources: https://stackoverflow.com/questions/2076475/reading-an-image-file-in-c-c, https://github.com/nothings/stb

// Specify filepath of image to be read, width/height of image, and a 2D vec of unsigned ints (any size).  Modifies provided 2D vector to contain pixels represented as unsigned ints.
// Note: This is designed only for colored, RGB PNG images with 3 bytes per pixel (hence width of 2D vec is width*3).
void image_read(char const *filepath, vector<vector<double> > &R, vector<vector<double> > &G, vector<vector<double> > &B) {
    int bpp,height,width;
    uint8_t* rgb_image = stbi_load(filepath, &width, &height, &bpp, 3);
    bpp = 3;
    R.resize(height); G.resize(height); B.resize(height);
    for (int i = 0; i < height; ++i) {
        R[i] = vector<double>(width,0);
        G[i] = vector<double>(width,0);
        B[i] = vector<double>(width,0);
    }
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width*bpp; j+=bpp) {
            R[i][j/bpp] = double(rgb_image[width*bpp*i+j]);
            G[i][j/bpp] = double(rgb_image[width*bpp*i+j+1]);
            B[i][j/bpp] = double(rgb_image[width*bpp*i+j+2]);
        }
    }
    stbi_image_free(rgb_image);
}

// Specify filepath where image should be written, and provide 2D vec of unsigned ints containing pixels of image.
// Note: This is designed only for colored, RGB PNG images with 3 bytes per pixel
void image_write(char const *filepath, vector<vector<double> > &R, vector<vector<double> > &G, vector<vector<double> > &B) {
    int width = R[0].size();
    int height = R.size();
    int bpp = 3;
    uint8_t* rgb_image;
    rgb_image = (uint8_t*) malloc(width*height*bpp);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            rgb_image[width*bpp*i+bpp*j] = uint8_t(R[i][j]);
            rgb_image[width*bpp*i+bpp*j+1] = uint8_t(G[i][j]);
            rgb_image[width*bpp*i+bpp*j+2] = uint8_t(B[i][j]);
        }
    }
    stbi_write_png(filepath, width, height, bpp, rgb_image, width*bpp);
    free(rgb_image);
}
