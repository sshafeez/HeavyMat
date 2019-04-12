#include <stdint.h>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image.h"
#include "stb_image_write.h"
#include <vector>

// Resources: https://stackoverflow.com/questions/2076475/reading-an-image-file-in-c-c, https://github.com/nothings/stb

// Specify filepath of image to be read, width/height of image, and a 2D vec of unsigned ints (any size).  Modifies provided 2D vector to contain pixels represented as unsigned ints.
// Note: This is designed only for colored, RGB PNG images with 3 bytes per pixel (hence width of 2D vec is width*3).
void image_read(char const *filepath, int width, int height, vector<vector<unsigned> > &vec) {
    int bpp = 3;
    uint8_t* rgb_image = stbi_load(filepath, &width, &height, &bpp, 3);
    vec.resize(height);
    for (int i = 0; i < vec.size(); ++i) {
        vec[i] = vector<unsigned>(width*bpp,0);
    }
    cout << vec.size() << " " << vec[0].size() << endl;
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width*bpp; ++j) {
            vec[i][j] = unsigned(rgb_image[width*bpp*i+j]);
        }
    }
    stbi_image_free(rgb_image);
}

// Specify filepath where image should be written, and provide 2D vec of unsigned ints containing pixels of image.
// Note: This is designed only for colored, RGB PNG images with 3 bytes per pixel
void image_write(char const *filepath, vector<vector<unsigned> > &vec) {
    int width = vec[0].size();
    int height = vec.size();
    cout << width << " " << height << endl;
    int bpp = 3;
    uint8_t* rgb_image;
    rgb_image = (uint8_t*) malloc(width*height*bpp);
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width*bpp; ++j) {
            rgb_image[width*bpp*i+j] = uint8_t(vec[i][j]);
        }
    }
    stbi_write_png(filepath, width/3, height, bpp, rgb_image, width*bpp);
    free(rgb_image);
}
