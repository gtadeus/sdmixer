#include "sdmixer.h"

void sortrows(Matrix& m, int column)
{
    Array<double> t = m.column(column);
    Array<octave_idx_type> idx;
    idx=t.sort_rows_idx();

    idx_vector i(idx);
    idx_vector j(octave_idx_type(0), octave_idx_type(m.cols()));

    m = m.index(i,j);

}

void setUpperMax(Matrix& m, int column, double max)
{

    Array<octave_idx_type> nz_idx =  mx_el_gt(octave_value(m.column(column)).array_value(), max).find();
    idx_vector i_nz(nz_idx);

    m.assign(i_nz, octave_value(max).array_value());

}
void deleteRowsThreshold(Matrix& m, int column, double min)
{

    Array<octave_idx_type> nz_idx =  mx_el_lt(octave_value(m.column(column)).array_value(), min).find();
    idx_vector i_nz(nz_idx);

    m.delete_elements(0, i_nz);

}
void deleteRowsMaxThreshold(Matrix& m, int column, double max)
{

    Array<octave_idx_type> nz_idx =  mx_el_gt(octave_value(m.column(column)).array_value(), max).find();
    idx_vector i_nz(nz_idx);

    m.delete_elements(0, i_nz);

}

Matrix unique_elements(const Matrix& m, int column)
{

    NDArray arr=m.column(column);
    double * pix_p=(double *)arr.fortran_vec();

    std::vector<double> tmp_vec (pix_p,pix_p+arr.length());
    std::sort( tmp_vec.begin(), tmp_vec.end() );

    tmp_vec.erase( std::unique( tmp_vec.begin(), tmp_vec.end() ), tmp_vec.end() );

    double* a = &tmp_vec[0];
    Matrix un(dim_vector (tmp_vec.size(), 1));

    double *pm = un.fortran_vec();
    memcpy(pm, a, tmp_vec.size()*sizeof(double));

    return un;

}

void find_position(const Matrix& m, int column, Matrix& un, Matrix& j)
{
    int n = m.rows();
    j.resize(n, 1, 0);

    un = m.column(column);
    Array<octave_idx_type> idx_sort;
    un=un.sort(idx_sort);
    idx_vector i_sort(idx_sort);

    Matrix match = mx_el_eq(un.extract(0,0,n-2,0), un.extract(1,0,n-1,0));
    Array<octave_idx_type> nz_idx = match.all(2).find();

    idx_vector i_nz(nz_idx);
    un.delete_elements(i_nz);

    j.assign(i_sort, octave_value(1).matrix_value().stack(!match).cumsum());

}


double cvsum(const ColumnVector& cv)
{
  return octave_value(do_mx_red_op<double, double> (cv, -1, mx_inline_sum)).double_value();
}

void get_statistics(const ColumnVector& cv)
{

    int n = cv.length();

    double mean = cvsum(cv) / double(n);
    double accum = 0.0;
    double skewness =  0.0;
    double kurtosis = 0.0;


    for (int i = 0; i < cv.length(); ++i)
    {
        accum += (cv(i) - mean) * (cv(i) - mean);
        skewness += (cv(i) - mean) * (cv(i) - mean)* (cv(i) - mean);
        kurtosis += (cv(i) - mean) * (cv(i) - mean) * (cv(i) - mean) * (cv(i) - mean);
    }

    double stdev =  sqrt(accum / (n-1));
    // Sort x
    //std::sort(x.begin(), x.end());
    ColumnVector cv2=cv.sort();


    double minimum = cv2(0);
    // double firstquartile = (x[(n - 1) / 4] + x[n / 4]) / 2.;
    double median = (cv2((n - 1) / 2) + cv2(n / 2))/ 2.;
    //double thirdquartile = (x[n - 1 - (n - 1) / 4] + x[n - 1 - n / 4]) / 2.;
    double maximum = cv2(n - 1);

    skewness*=(1/(double)n)*pow(stdev,-3);
    kurtosis*=(1/(double)n)*pow(stdev,-4);
    kurtosis-=3;

    std::cout << "N:   " << n << std::endl;
    std::cout << "min:   " << minimum << std::endl;
    std::cout << "median:   " << median << std::endl;
    std::cout << "max:   " << maximum << std::endl;
    std::cout << "mean:   " << mean << std::endl;
    std::cout << "stdev:   " << stdev << std::endl;
    std::cout << "skewness:   " << skewness << std::endl;
    std::cout << "kurtosis:   " << kurtosis << std::endl;
}


Matrix accum_array(const Matrix& j)
{

    idx_vector idx = octave_value(j).index_vector();
    octave_idx_type n = idx.extent (0);
    Matrix out(n, 1, 0);
    out.idx_add(idx, 1);

    return out;
}

int load_file(const std::string& filename, Matrix& m)
{
    std::filebuf fb;
    octave_value retval;
    fb.open (filename.c_str(), std::ios::in);
    std::istream is(&fb);
    read_mat_ascii_data(is, filename, retval);
    m = retval.matrix_value();

    fb.close();

    return 0;
}

Matrix coord2indx(dim_vector& dv, Matrix& coord)
{

    int dim = coord.cols();

    Matrix retval;

    Array<idx_vector> idxa (dim_vector (dim, 1));

    for (int j = 0; j < dim; j++)
        idxa(j)=octave_value(coord.column(j)).index_vector();

    idx_vector idx = sub2ind (dv, idxa);

    return octave_value(idx).matrix_value();


}
