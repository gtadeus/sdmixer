#include "sdmixer.h"

void sortrows(Matrix& m, int column);
void setUpperMax(Matrix& m, int column, double max);
Matrix unique_elements(const Matrix& m, int column);
void find_position(const Matrix& m, int column, Matrix& un, Matrix& j);
double cvsum(const ColumnVector& cv);
void get_statistics(const ColumnVector& cv);
Matrix accum_array(const Matrix& j);
Matrix coord2indx(dim_vector& dv, Matrix& coord);
int load_file(const std::string& filename, Matrix& m);
