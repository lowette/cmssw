#include "RunSteps/Clusterizer/interface/SiPixelArrayBuffer.h"

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include <vector>

// Create a 2D matrix of zeros of size rows x cols
SiPixelArrayBuffer::SiPixelArrayBuffer(int rows, int cols) :
    nrows(rows),
    ncols(cols) {
    pixel_vec.resize(rows * cols);
    for (std::vector<int>::iterator it = pixel_vec.begin(); it != pixel_vec.end(); ++it) *it = 0;
}

// Change the size of the matrix and reset it
void SiPixelArrayBuffer::setSize(int rows, int cols) {
    nrows = rows;
    ncols = cols;
    pixel_vec.resize(rows * cols);
    for (std::vector<int>::iterator it = pixel_vec.begin(); it != pixel_vec.end(); ++it) *it = 0;
}

// Check if an element is inside the matrix (no overflow)
bool SiPixelArrayBuffer::inside(int row, int col) const {
    return (row >= 0 && row < nrows && col >= 0 && col < ncols);
}

// Return an element of the matrix
int SiPixelArrayBuffer::operator()(int row, int col) const {
    if (inside(row,col)) return pixel_vec[index(row,col)];
    else return 0;
}

// Return an element of the matrix (get the row and column from the pixel data)
int SiPixelArrayBuffer::operator()(const SiPixelCluster::PixelPos& pix) const {
    if (inside(pix.row(), pix.col())) return pixel_vec[index(pix)];
    else return 0;
}

// Set an element of the matrix
void SiPixelArrayBuffer::set(int row, int col, int adc) {
    pixel_vec[index(row,col)] = adc;
}

// Set an element of the matrix (get the row and column from the pixel data)
void SiPixelArrayBuffer::set(const SiPixelCluster::PixelPos& pix, int adc) {
    pixel_vec[index(pix)] = adc;
}
