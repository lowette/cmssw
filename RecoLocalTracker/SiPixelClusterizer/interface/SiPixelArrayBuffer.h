#ifndef RecoLocalTracker_SiPixelClusterizer_SiPixelArrayBuffer_H
#define RecoLocalTracker_SiPixelClusterizer_SiPixelArrayBuffer_H

#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"

#include <vector>
#include <iostream>

class SiPixelArrayBuffer {

public:
  	inline SiPixelArrayBuffer(int rows, int cols);
  	inline SiPixelArrayBuffer() { }
  
  	inline void setSize(int rows, int cols);
 	inline int operator()(int row, int col) const;
  	inline int operator()(const SiPixelCluster::PixelPos&) const;
  	inline int rows() const { return nrows; }
  	inline int columns() const { return ncols; }

  	inline bool inside(int row, int col) const;
  	inline void set_adc(int row, int col, int adc);
  	inline void set_adc(const SiPixelCluster::PixelPos&, int adc);
  	int size() const { return pixel_vec.size(); }

  	int index(int row, int col) const { return col * nrows + row; }
  	int index(const SiPixelCluster::PixelPos& pix) const { return index(pix.row(), pix.col()); }

private:
  	int nrows;
  	int ncols;
  	std::vector<int> pixel_vec;
};

SiPixelArrayBuffer::SiPixelArrayBuffer(int rows, int cols) : 
	nrows(rows), 
	ncols(cols) {
  	pixel_vec.resize(rows * cols);
  	for (std::vector<int>::iterator it = pixel_vec.begin(); it != pixel_vec.end(); ++it) *it = 0;
}

void SiPixelArrayBuffer::setSize(int rows, int cols) {
  	nrows = rows;
  	ncols = cols;
  	pixel_vec.resize(rows * cols);
	for (std::vector<int>::iterator it = pixel_vec.begin(); it != pixel_vec.end(); ++it) *it = 0;
}

bool SiPixelArrayBuffer::inside(int row, int col) const {
	return (row >= 0 && row < nrows && col >= 0 && col < ncols);
}

int SiPixelArrayBuffer::operator()(int row, int col) const {
  	if (inside(row,col)) return pixel_vec[index(row,col)];
  	else return 0;
}

int SiPixelArrayBuffer::operator()(const SiPixelCluster::PixelPos& pix) const {
  	if (inside(pix.row(), pix.col())) return pixel_vec[index(pix)];
  	else return 0;
}

void SiPixelArrayBuffer::set_adc(int row, int col, int adc) {
	pixel_vec[index(row,col)] = adc;
}

void SiPixelArrayBuffer::set_adc(const SiPixelCluster::PixelPos& pix, int adc) {
  	pixel_vec[index(pix)] = adc;
}

#endif
