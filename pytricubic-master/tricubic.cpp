#include "tricubic.hpp"

//This code is adapted from https://github.com/deepzot/likely
TriCubicInterpolator::TriCubicInterpolator(list data, list nkpoints){

  
  _initialized = false;
  _spacing = 1.0;
  _n1 = extract<int>(nkpoints[0]);
  _n2 = extract<int>(nkpoints[1]);
  _n3 = extract<int>(nkpoints[2]);
  _data.resize(boost::extents[_n1*_n2*_n3]);
  for(int k=0;k<_n3;k++){
    for(int j=0;j<_n2;j++){
      for(int i=0;i<_n1;i++)
	_data[_index(i,j,k)] = extract<fptype>(data[i][j][k]);
    }
  }
  
  //temporary array is necessary, otherwise compiler has problems with Eigen and takes very long to compile
  int temp[64][64] = {
    { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 9,-9,-9, 9, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 6,-6, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 6,-6, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 4,-4,-4, 4, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0},
    {-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 9,-9, 0, 0,-9, 9, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 0, 0, 6,-6, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9, 0, 0,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0},
    { 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0},
    {-27,27,27,-27,27,-27,-27,27,-18,-9,18, 9,18, 9,-18,-9,-18,18,-9, 9,18,-18, 9,-9,-18,18,18,-18,-9, 9, 9,-9,-12,-6,-6,-3,12, 6, 6, 3,-12,-6,12, 6,-6,-3, 6, 3,-12,12,-6, 6,-6, 6,-3, 3,-8,-4,-4,-2,-4,-2,-2,-1},
    {18,-18,-18,18,-18,18,18,-18, 9, 9,-9,-9,-9,-9, 9, 9,12,-12, 6,-6,-12,12,-6, 6,12,-12,-12,12, 6,-6,-6, 6, 6, 6, 3, 3,-6,-6,-3,-3, 6, 6,-6,-6, 3, 3,-3,-3, 8,-8, 4,-4, 4,-4, 2,-2, 4, 4, 2, 2, 2, 2, 1, 1},
    {-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0},
    {18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6, 9,-9, 9,-9,-9, 9,-9, 9,12,-12,-12,12, 6,-6,-6, 6, 6, 3, 6, 3,-6,-3,-6,-3, 8, 4,-8,-4, 4, 2,-4,-2, 6,-6, 6,-6, 3,-3, 3,-3, 4, 2, 4, 2, 2, 1, 2, 1},
    {-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-6, 6,-6, 6, 6,-6, 6,-6,-8, 8, 8,-8,-4, 4, 4,-4,-3,-3,-3,-3, 3, 3, 3, 3,-4,-4, 4, 4,-2,-2, 2, 2,-4, 4,-4, 4,-2, 2,-2, 2,-2,-2,-2,-2,-1,-1,-1,-1},
    { 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {-6, 6, 0, 0, 6,-6, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 4,-4, 0, 0,-4, 4, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4, 0, 0,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0},
    {-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0},
    {18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6,12,-12, 6,-6,-12,12,-6, 6, 9,-9,-9, 9, 9,-9,-9, 9, 8, 4, 4, 2,-8,-4,-4,-2, 6, 3,-6,-3, 6, 3,-6,-3, 6,-6, 3,-3, 6,-6, 3,-3, 4, 2, 2, 1, 4, 2, 2, 1},
    {-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-8, 8,-4, 4, 8,-8, 4,-4,-6, 6, 6,-6,-6, 6, 6,-6,-4,-4,-2,-2, 4, 4, 2, 2,-3,-3, 3, 3,-3,-3, 3, 3,-4, 4,-2, 2,-4, 4,-2, 2,-2,-2,-1,-1,-2,-2,-1,-1},
    { 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    { 0, 0, 0, 0, 0, 0, 0, 0, 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0},
    {-12,12,12,-12,12,-12,-12,12,-8,-4, 8, 4, 8, 4,-8,-4,-6, 6,-6, 6, 6,-6, 6,-6,-6, 6, 6,-6,-6, 6, 6,-6,-4,-2,-4,-2, 4, 2, 4, 2,-4,-2, 4, 2,-4,-2, 4, 2,-3, 3,-3, 3,-3, 3,-3, 3,-2,-1,-2,-1,-2,-1,-2,-1},
    { 8,-8,-8, 8,-8, 8, 8,-8, 4, 4,-4,-4,-4,-4, 4, 4, 4,-4, 4,-4,-4, 4,-4, 4, 4,-4,-4, 4, 4,-4,-4, 4, 2, 2, 2, 2,-2,-2,-2,-2, 2, 2,-2,-2, 2, 2,-2,-2, 2,-2, 2,-2, 2,-2, 2,-2, 1, 1, 1, 1, 1, 1, 1, 1}
  };
  
  for(int i=0;i<64;i++){
    for(int j=0;j<64;j++){
      _C(i,j) = temp[i][j];
    }
  }
}

fptype TriCubicInterpolator::ip(list xyz){
  
  fptype x = extract<fptype>(xyz[0]);
  fptype y = extract<fptype>(xyz[1]);
  fptype z = extract<fptype>(xyz[2]);
  
  fptype dx = fmod(x/_spacing, _n1), dy = fmod(y/_spacing, _n2), dz = fmod(z/_spacing, _n3); //determine the relative position in the box enclosed by nearest data points
  
  if(dx < 0) dx += _n1; //periodicity is built in
  if(dy < 0) dy += _n2;
  if(dz < 0) dz += _n3;
  
  int xi = (int)floor(dx); //calculate lower-bound grid indices
  int yi = (int)floor(dy);
  int zi = (int)floor(dz);
  
  // Check if we can re-use coefficients from the last interpolation.
  if(!_initialized || xi != _i1 || yi != _i2 || zi != _i3) {
  // Extract the local vocal values and calculate partial derivatives.
  Eigen::Matrix<fptype,64,1> x;
  x << 
      // values of f(x,y,z) at each corner.
      _data[_index(xi,yi,zi)],_data[_index(xi+1,yi,zi)],_data[_index(xi,yi+1,zi)],
      _data[_index(xi+1,yi+1,zi)],_data[_index(xi,yi,zi+1)],_data[_index(xi+1,yi,zi+1)],
      _data[_index(xi,yi+1,zi+1)],_data[_index(xi+1,yi+1,zi+1)],
      // values of df/dx at each corner.
      0.5*(_data[_index(xi+1,yi,zi)]-_data[_index(xi-1,yi,zi)]),
      0.5*(_data[_index(xi+2,yi,zi)]-_data[_index(xi,yi,zi)]),
      0.5*(_data[_index(xi+1,yi+1,zi)]-_data[_index(xi-1,yi+1,zi)]),
      0.5*(_data[_index(xi+2,yi+1,zi)]-_data[_index(xi,yi+1,zi)]),
      0.5*(_data[_index(xi+1,yi,zi+1)]-_data[_index(xi-1,yi,zi+1)]),
      0.5*(_data[_index(xi+2,yi,zi+1)]-_data[_index(xi,yi,zi+1)]),
      0.5*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi-1,yi+1,zi+1)]),
      0.5*(_data[_index(xi+2,yi+1,zi+1)]-_data[_index(xi,yi+1,zi+1)]),
      // values of df/dy at each corner.
      0.5*(_data[_index(xi,yi+1,zi)]-_data[_index(xi,yi-1,zi)]),
      0.5*(_data[_index(xi+1,yi+1,zi)]-_data[_index(xi+1,yi-1,zi)]),
      0.5*(_data[_index(xi,yi+2,zi)]-_data[_index(xi,yi,zi)]),
      0.5*(_data[_index(xi+1,yi+2,zi)]-_data[_index(xi+1,yi,zi)]),
      0.5*(_data[_index(xi,yi+1,zi+1)]-_data[_index(xi,yi-1,zi+1)]),
      0.5*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi+1,yi-1,zi+1)]),
      0.5*(_data[_index(xi,yi+2,zi+1)]-_data[_index(xi,yi,zi+1)]),
      0.5*(_data[_index(xi+1,yi+2,zi+1)]-_data[_index(xi+1,yi,zi+1)]),
      // values of df/dz at each corner.
      0.5*(_data[_index(xi,yi,zi+1)]-_data[_index(xi,yi,zi-1)]),
      0.5*(_data[_index(xi+1,yi,zi+1)]-_data[_index(xi+1,yi,zi-1)]),
      0.5*(_data[_index(xi,yi+1,zi+1)]-_data[_index(xi,yi+1,zi-1)]),
      0.5*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi+1,yi+1,zi-1)]),
      0.5*(_data[_index(xi,yi,zi+2)]-_data[_index(xi,yi,zi)]),
      0.5*(_data[_index(xi+1,yi,zi+2)]-_data[_index(xi+1,yi,zi)]),
      0.5*(_data[_index(xi,yi+1,zi+2)]-_data[_index(xi,yi+1,zi)]),
      0.5*(_data[_index(xi+1,yi+1,zi+2)]-_data[_index(xi+1,yi+1,zi)]),
      // values of d2f/dxdy at each corner.
      0.25*(_data[_index(xi+1,yi+1,zi)]-_data[_index(xi-1,yi+1,zi)]-_data[_index(xi+1,yi-1,zi)]+_data[_index(xi-1,yi-1,zi)]),
      0.25*(_data[_index(xi+2,yi+1,zi)]-_data[_index(xi,yi+1,zi)]-_data[_index(xi+2,yi-1,zi)]+_data[_index(xi,yi-1,zi)]),
      0.25*(_data[_index(xi+1,yi+2,zi)]-_data[_index(xi-1,yi+2,zi)]-_data[_index(xi+1,yi,zi)]+_data[_index(xi-1,yi,zi)]),
      0.25*(_data[_index(xi+2,yi+2,zi)]-_data[_index(xi,yi+2,zi)]-_data[_index(xi+2,yi,zi)]+_data[_index(xi,yi,zi)]),
      0.25*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi-1,yi+1,zi+1)]-_data[_index(xi+1,yi-1,zi+1)]+_data[_index(xi-1,yi-1,zi+1)]),
      0.25*(_data[_index(xi+2,yi+1,zi+1)]-_data[_index(xi,yi+1,zi+1)]-_data[_index(xi+2,yi-1,zi+1)]+_data[_index(xi,yi-1,zi+1)]),
      0.25*(_data[_index(xi+1,yi+2,zi+1)]-_data[_index(xi-1,yi+2,zi+1)]-_data[_index(xi+1,yi,zi+1)]+_data[_index(xi-1,yi,zi+1)]),
      0.25*(_data[_index(xi+2,yi+2,zi+1)]-_data[_index(xi,yi+2,zi+1)]-_data[_index(xi+2,yi,zi+1)]+_data[_index(xi,yi,zi+1)]),
      // values of d2f/dxdz at each corner.
      0.25*(_data[_index(xi+1,yi,zi+1)]-_data[_index(xi-1,yi,zi+1)]-_data[_index(xi+1,yi,zi-1)]+_data[_index(xi-1,yi,zi-1)]),
      0.25*(_data[_index(xi+2,yi,zi+1)]-_data[_index(xi,yi,zi+1)]-_data[_index(xi+2,yi,zi-1)]+_data[_index(xi,yi,zi-1)]),
      0.25*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi-1,yi+1,zi+1)]-_data[_index(xi+1,yi+1,zi-1)]+_data[_index(xi-1,yi+1,zi-1)]),
      0.25*(_data[_index(xi+2,yi+1,zi+1)]-_data[_index(xi,yi+1,zi+1)]-_data[_index(xi+2,yi+1,zi-1)]+_data[_index(xi,yi+1,zi-1)]),
      0.25*(_data[_index(xi+1,yi,zi+2)]-_data[_index(xi-1,yi,zi+2)]-_data[_index(xi+1,yi,zi)]+_data[_index(xi-1,yi,zi)]),
      0.25*(_data[_index(xi+2,yi,zi+2)]-_data[_index(xi,yi,zi+2)]-_data[_index(xi+2,yi,zi)]+_data[_index(xi,yi,zi)]),
      0.25*(_data[_index(xi+1,yi+1,zi+2)]-_data[_index(xi-1,yi+1,zi+2)]-_data[_index(xi+1,yi+1,zi)]+_data[_index(xi-1,yi+1,zi)]),
      0.25*(_data[_index(xi+2,yi+1,zi+2)]-_data[_index(xi,yi+1,zi+2)]-_data[_index(xi+2,yi+1,zi)]+_data[_index(xi,yi+1,zi)]),
      // values of d2f/dydz at each corner.
      0.25*(_data[_index(xi,yi+1,zi+1)]-_data[_index(xi,yi-1,zi+1)]-_data[_index(xi,yi+1,zi-1)]+_data[_index(xi,yi-1,zi-1)]),
      0.25*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi+1,yi-1,zi+1)]-_data[_index(xi+1,yi+1,zi-1)]+_data[_index(xi+1,yi-1,zi-1)]),
      0.25*(_data[_index(xi,yi+2,zi+1)]-_data[_index(xi,yi,zi+1)]-_data[_index(xi,yi+2,zi-1)]+_data[_index(xi,yi,zi-1)]),
      0.25*(_data[_index(xi+1,yi+2,zi+1)]-_data[_index(xi+1,yi,zi+1)]-_data[_index(xi+1,yi+2,zi-1)]+_data[_index(xi+1,yi,zi-1)]),
      0.25*(_data[_index(xi,yi+1,zi+2)]-_data[_index(xi,yi-1,zi+2)]-_data[_index(xi,yi+1,zi)]+_data[_index(xi,yi-1,zi)]),
      0.25*(_data[_index(xi+1,yi+1,zi+2)]-_data[_index(xi+1,yi-1,zi+2)]-_data[_index(xi+1,yi+1,zi)]+_data[_index(xi+1,yi-1,zi)]),
      0.25*(_data[_index(xi,yi+2,zi+2)]-_data[_index(xi,yi,zi+2)]-_data[_index(xi,yi+2,zi)]+_data[_index(xi,yi,zi)]),
      0.25*(_data[_index(xi+1,yi+2,zi+2)]-_data[_index(xi+1,yi,zi+2)]-_data[_index(xi+1,yi+2,zi)]+_data[_index(xi+1,yi,zi)]),
      // values of d3f/dxdydz at each corner.
      0.125*(_data[_index(xi+1,yi+1,zi+1)]-_data[_index(xi-1,yi+1,zi+1)]-_data[_index(xi+1,yi-1,zi+1)]+_data[_index(xi-1,yi-1,zi+1)]-_data[_index(xi+1,yi+1,zi-1)]+_data[_index(xi-1,yi+1,zi-1)]+_data[_index(xi+1,yi-1,zi-1)]-_data[_index(xi-1,yi-1,zi-1)]),
      0.125*(_data[_index(xi+2,yi+1,zi+1)]-_data[_index(xi,yi+1,zi+1)]-_data[_index(xi+2,yi-1,zi+1)]+_data[_index(xi,yi-1,zi+1)]-_data[_index(xi+2,yi+1,zi-1)]+_data[_index(xi,yi+1,zi-1)]+_data[_index(xi+2,yi-1,zi-1)]-_data[_index(xi,yi-1,zi-1)]),
      0.125*(_data[_index(xi+1,yi+2,zi+1)]-_data[_index(xi-1,yi+2,zi+1)]-_data[_index(xi+1,yi,zi+1)]+_data[_index(xi-1,yi,zi+1)]-_data[_index(xi+1,yi+2,zi-1)]+_data[_index(xi-1,yi+2,zi-1)]+_data[_index(xi+1,yi,zi-1)]-_data[_index(xi-1,yi,zi-1)]),
      0.125*(_data[_index(xi+2,yi+2,zi+1)]-_data[_index(xi,yi+2,zi+1)]-_data[_index(xi+2,yi,zi+1)]+_data[_index(xi,yi,zi+1)]-_data[_index(xi+2,yi+2,zi-1)]+_data[_index(xi,yi+2,zi-1)]+_data[_index(xi+2,yi,zi-1)]-_data[_index(xi,yi,zi-1)]),
      0.125*(_data[_index(xi+1,yi+1,zi+2)]-_data[_index(xi-1,yi+1,zi+2)]-_data[_index(xi+1,yi-1,zi+2)]+_data[_index(xi-1,yi-1,zi+2)]-_data[_index(xi+1,yi+1,zi)]+_data[_index(xi-1,yi+1,zi)]+_data[_index(xi+1,yi-1,zi)]-_data[_index(xi-1,yi-1,zi)]),
      0.125*(_data[_index(xi+2,yi+1,zi+2)]-_data[_index(xi,yi+1,zi+2)]-_data[_index(xi+2,yi-1,zi+2)]+_data[_index(xi,yi-1,zi+2)]-_data[_index(xi+2,yi+1,zi)]+_data[_index(xi,yi+1,zi)]+_data[_index(xi+2,yi-1,zi)]-_data[_index(xi,yi-1,zi)]),
      0.125*(_data[_index(xi+1,yi+2,zi+2)]-_data[_index(xi-1,yi+2,zi+2)]-_data[_index(xi+1,yi,zi+2)]+_data[_index(xi-1,yi,zi+2)]-_data[_index(xi+1,yi+2,zi)]+_data[_index(xi-1,yi+2,zi)]+_data[_index(xi+1,yi,zi)]-_data[_index(xi-1,yi,zi)]),
      0.125*(_data[_index(xi+2,yi+2,zi+2)]-_data[_index(xi,yi+2,zi+2)]-_data[_index(xi+2,yi,zi+2)]+_data[_index(xi,yi,zi+2)]-_data[_index(xi+2,yi+2,zi)]+_data[_index(xi,yi+2,zi)]+_data[_index(xi+2,yi,zi)]-_data[_index(xi,yi,zi)])
    ;
    // Convert voxel values and partial derivatives to interpolation coefficients.
    _coefs = _C * x;
    // Remember this voxel for next time.
    _i1 = xi;
    _i2 = yi;
    _i3 = zi;
    _initialized = true;
  }
  // Evaluate the interpolation within this grid voxel.
  dx -= xi;
  dy -= yi;
  dz -= zi;
  int ijkn(0);
  fptype dzpow(1);
  fptype result(0);
  for(int k = 0; k < 4; ++k) {
    fptype dypow(1);
    for(int j = 0; j < 4; ++j) {
      result += dypow*dzpow*(_coefs[ijkn] + dx*(_coefs[ijkn+1] + dx*(_coefs[ijkn+2] + dx*_coefs[ijkn+3])));
      ijkn += 4;
      dypow *= dy;
    }
    dzpow *= dz;
  }
  return result;
}

BOOST_PYTHON_MODULE(tricubic)
{
  class_<TriCubicInterpolator>("tricubic", init<list, list>())
    .def("ip", &TriCubicInterpolator::ip)
  ;
}
