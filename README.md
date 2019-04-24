## HEAVY MAT
A memory intensive matrix class with fast matrix multiplication

## Files
Samples folder contains images with flexible basis technique applied at various maximum error levels <br />
Samples_SVD folder contains image with SVD technique applied to save same # of mults as corresponding images in samples folder <br />
matrix.h contains the main matrix class, and fundamental matrix operations, along with naive multiplication, Strassen multiplication, and flexible basis multiplication methods implemented <br />
SVD_Matrix.h, SVD_Matrix.cpp, SVD.h, SVD.cpp, and Vector.h are adaptations of code originating online (source listed below) to do singular value decomposition on square matrices <br />
M.png, arches.png, north.png are the three images used for the poster.  We also analyzed sign.png but it was not used in our poster <br />
stb_image.h and stb_image_write.h convert images into an array of bytes <br />
image.h contains functions that turn image data (in form of array of bytes) into a 2D vector of double, or an array of floats, for both grayscale and RGB images, and vice versa <br />
 
All resources used:
Strassen Multiplication Method: https://www.youtube.com/watch?v=1AIvlizGo7Y <br />
Strassen Multiplication: https://www.geeksforgeeks.org/strassens-matrix-multiplication/ <br />
SVD code: http://svn.lirec.eu/libs/magicsquares/src/SVD.cpp <br />
Mlogo image: https://etaluma.com/company/our-customers/umich-logo/ <br />
Other images: https://www.pinterest.com/uofmichigan/tour-of-umich/ <br />


