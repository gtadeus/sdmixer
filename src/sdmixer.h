//=================================
// include guard
#ifndef __SDMIXER_H_INCLUDED__
#define __SDMIXER_H_INCLUDED__

#include <iostream>
#include <cstdlib>
#include <sys/types.h>
#include <iomanip>
#include <dirent.h>
#include <errno.h>
#include <vector>
#include <string>
#include <sstream>
#include <stdio.h>
#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/regex.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/string_path.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
//#include "LogFile.hpp"

#include <time.h>
#include <omp.h>
#include <cmath>

#include <parallel/algorithm>
#include <parallel/settings.h>

//#include "iniReader.h"	//  * Author: benmaynard, under GPL
			//  * Created on August 26, 2010, 2:49 PM
			// http://code.google.com/p/xplugins/source/browse/

#include <octave/oct.h>
#include <octave/load-save.h>     // Load/Save
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/dMatrix.h>       // Real double precision matrices

extern std::ofstream logFile;

#include <Magick++.h>
#include <clocale>
#include <cmath>
#include "OctaveTools.hpp"
#include <tiffio.h>

template <class T> static void
encode_uint_image (std::vector<Magick::Image>& imvec, const octave_value& img, bool has_map)
{
    unsigned int bitdepth = 8;
    T m;
    m = img.array_value();
    dim_vector dsizes = m.dims ();
    unsigned int nframes = 1;
    /*if (dsizes.length () == 4)
        nframes = dsizes(3);
*/
    if ( dsizes.length() == 3)
        nframes = dsizes(2);

    Array<octave_idx_type> idx (dim_vector (dsizes.length (), 1));
    octave_idx_type rows = m.rows ();
    octave_idx_type columns = m.columns ();

    unsigned int div_factor = (1 << bitdepth) - 1;
   /* logFile << "rows: " << rows << std::endl;
    logFile << "columns: " << columns << std::endl;
    logFile << "init ready, bitdepth = " << bitdepth << std::endl;
    logFile << "nframes = " << nframes << std::endl;*/

    for (unsigned int ii = 0; ii < nframes; ii++)
    {
        if(ii > 1)
            if( (ii % (int)(nframes/5)) == 0)
                std::cout << "writing frame " << ii << "/" << nframes << std::endl;

        Magick::Image im(Magick::Geometry (columns, rows), "black");
        im.depth (bitdepth);
        //im.density(Magick::Geometry (1, 1));
        im.compressType(Magick::LZWCompression);
        im.endian(Magick::MSBEndian);
        //im.resolutionUnits(Magick::UndefinedResolution);
        im.classType (Magick::PseudoClass);
        im.type (Magick::GrayscaleType);

        Magick::ColorGray c;
        for (int y = 0; y < columns; y++)
        {
            idx(1) = y;
            for (int x=0; x < rows; x++)
            {
                idx(0) = x;
                if (nframes > 1)
                {
                    //idx(2) = 0;
                    //idx(3) = ii;
                    idx(2) = ii;
                }
                c.shade (static_cast<double>(m(idx)));
                try
                {
                    im.pixelColor (y, x, c);

                }
                catch (Magick::Exception& e)
                {
                        error ("Magick++ exception: %s", e.what ());
                }

            }
        }
    imvec.push_back (im);
    }
}

template <class T>
octave_value_list
read_images (const std::vector<Magick::Image>& imvec,
             const Array<int>& frameidx, unsigned int depth)
{
  typedef typename T::element_type P;

  octave_value_list retval (3, Matrix ());

  T im;

  int rows = imvec[0].baseRows ();
  int columns = imvec[0].baseColumns ();
  int nframes = frameidx.length ();

  dim_vector idim = dim_vector ();
  idim.resize (4);
  idim(0) = rows;
  idim(1) = columns;
  idim(2) = 1;
  idim(3) = nframes;

  Magick::ImageType type = imvec[0].type ();
  const int divisor = ((uint64_t (1) << QuantumDepth) - 1) /
                      ((uint64_t (1) << depth) - 1);

  switch (type)
    {
    case Magick::BilevelType:
    case Magick::GrayscaleType:
      {
        im = T (idim);
        P *vec = im.fortran_vec ();

        for (int frame = 0; frame < nframes; frame++)
          {
            const Magick::PixelPacket *pix
              = imvec[frameidx(frame)].getConstPixels (0, 0, columns, rows);

            P *rbuf = vec;
            for (int y = 0; y < rows; y++)
              {
                for (int x = 0; x < columns; x++)
                  {
                    *rbuf = pix->red / divisor;
                    pix++;
                    rbuf += rows;
                  }
                rbuf -= rows * columns - 1;
              }

            // Next frame.
            vec += rows * columns * idim(2);
          }
        }
      break;

    case Magick::GrayscaleMatteType:
      {
        idim(2) = 2;
        im = T (idim);
        P *vec = im.fortran_vec ();

        for (int frame = 0; frame < nframes; frame++)
          {
            const Magick::PixelPacket *pix
              = imvec[frameidx(frame)].getConstPixels (0, 0, columns, rows);

            P *rbuf = vec;
            P *obuf = vec + rows * columns;
            for (int y = 0; y < rows; y++)
              {
                for (int x = 0; x < columns; x++)
                  {
                    *rbuf = pix->red / divisor;
                    *obuf = pix->opacity / divisor;
                    pix++;
                    rbuf += rows;
                    obuf += rows;
                  }
                rbuf -= rows * columns - 1;
                obuf -= rows * columns - 1;
              }

            // Next frame.
            vec += rows * columns * idim(2);
          }
        }
      break;

    case Magick::PaletteType:
    case Magick::TrueColorType:
      {
        idim(2) = 3;
        im = T (idim);
        P *vec = im.fortran_vec ();

        for (int frame = 0; frame < nframes; frame++)
          {
            const Magick::PixelPacket *pix
              = imvec[frameidx(frame)].getConstPixels (0, 0, columns, rows);

            P *rbuf = vec;
            P *gbuf = vec + rows * columns;
            P *bbuf = vec + rows * columns * 2;
            for (int y = 0; y < rows; y++)
              {
                for (int x = 0; x < columns; x++)
                  {
                    *rbuf = pix->red / divisor;
                    *gbuf = pix->green / divisor;
                    *bbuf = pix->blue / divisor;
                    pix++;
                    rbuf += rows;
                    gbuf += rows;
                    bbuf += rows;
                  }
                rbuf -= rows * columns - 1;
                gbuf -= rows * columns - 1;
                bbuf -= rows * columns - 1;
              }

            // Next frame.
            vec += rows * columns * idim(2);
          }
        }
      break;

    case Magick::PaletteMatteType:
    case Magick::TrueColorMatteType:
    case Magick::ColorSeparationType:
      {
        idim(2) = 4;
        im = T (idim);
        P *vec = im.fortran_vec ();

        for (int frame = 0; frame < nframes; frame++)
          {
            const Magick::PixelPacket *pix
              = imvec[frameidx(frame)].getConstPixels (0, 0, columns, rows);

            P *rbuf = vec;
            P *gbuf = vec + rows * columns;
            P *bbuf = vec + rows * columns * 2;
            P *obuf = vec + rows * columns * 3;
            for (int y = 0; y < rows; y++)
              {
                for (int x = 0; x < columns; x++)
                  {
                    *rbuf = pix->red / divisor;
                    *gbuf = pix->green / divisor;
                    *bbuf = pix->blue / divisor;
                    *obuf = pix->opacity / divisor;
                    pix++;
                    rbuf += rows;
                    gbuf += rows;
                    bbuf += rows;
                    obuf += rows;
                  }
                rbuf -= rows * columns - 1;
                gbuf -= rows * columns - 1;
                bbuf -= rows * columns - 1;
                obuf -= rows * columns - 1;
              }

            // Next frame.
            vec += rows * columns * idim(2);
          }
        }
      break;

    default:
      error ("__magick_read__: undefined ImageMagick image type");
      return retval;
    }

  retval(0) = im;

  return retval;
}
#include <octave/f77-fcn.h>

#include <octave/oct-convn.h>
#include <octave/oct-locbuf.h>
template <class T, class R>
static void
convolve_2d (const T *a, octave_idx_type ma, octave_idx_type na,
             const R *b, octave_idx_type mb, octave_idx_type nb,
             T *c, bool inner);

// Forward instances to our Fortran implementations.
#define FORWARD_IMPL(T,R,f,F) \
extern "C" \
F77_RET_T \
F77_FUNC (f##conv2o, F##CONV2O) (const octave_idx_type&, \
                                 const octave_idx_type&, \
                                 const T*, const octave_idx_type&, \
                                 const octave_idx_type&, const R*, T *); \
\
extern "C" \
F77_RET_T \
F77_FUNC (f##conv2i, F##CONV2I) (const octave_idx_type&, \
                                 const octave_idx_type&, \
                                 const T*, const octave_idx_type&, \
                                 const octave_idx_type&, const R*, T *); \
\
template <> void \
convolve_2d<T, R> (const T *a, octave_idx_type ma, octave_idx_type na, \
                   const R *b, octave_idx_type mb, octave_idx_type nb, \
                   T *c, bool inner) \
{ \
  if (inner) \
    F77_XFCN (f##conv2i, F##CONV2I, (ma, na, a, mb, nb, b, c)); \
  else \
    F77_XFCN (f##conv2o, F##CONV2O, (ma, na, a, mb, nb, b, c)); \
}
//FORWARD_IMPL (uint16, uint16, i, I)
FORWARD_IMPL (double, double, d, D)
FORWARD_IMPL (float, float, s, S)
FORWARD_IMPL (Complex, Complex, z, Z)
FORWARD_IMPL (FloatComplex, FloatComplex, c, C)
FORWARD_IMPL (Complex, double, zd, ZD)
FORWARD_IMPL (FloatComplex, float, cs, CS)
template <class T, class R>
void convolve_nd (const T *a, const dim_vector& ad, const dim_vector& acd,
                  const R *b, const dim_vector& bd, const dim_vector& bcd,
                  T *c, const dim_vector& ccd, int nd, bool inner)
{
  if (nd == 2)
    convolve_2d<T, R> (a, ad(0), ad(1), b, bd(0), bd(1), c, inner);
    //convole
  else
    {
      octave_idx_type ma = acd(nd-2), na = ad(nd-1), mb = bcd(nd-2), nb = bd(nd-1);
      octave_idx_type ldc = ccd(nd-2);
      if (inner)
        {
          for (octave_idx_type ja = 0; ja < na - nb + 1; ja++)
            for (octave_idx_type jb = 0; jb < nb; jb++)
              convolve_nd<T, R> (a + ma*(ja + jb), ad, acd, b + mb*jb, bd, bcd,
                                 c + ldc*ja, ccd, nd-1, inner);
        }
      else
        {
          for (octave_idx_type ja = 0; ja < na; ja++)
          {
              if(ja > 1)
                 if( (ja % (int)(na/5)) == 0)
                    std::cout << "running 3D convolution  " << ((float)ja/(float)na)*100 << "% .." << std::endl;

            for (octave_idx_type jb = 0; jb < nb; jb++)
              convolve_nd<T, R> (a + ma*ja, ad, acd, b + mb*jb, bd, bcd,
                                 c + ldc*(ja+jb), ccd, nd-1, inner);
          }
        }
    }
}

// Arbitrary convolutor.
// The 2nd array is assumed to be the smaller one.
template <class T, class R>
static MArray<T>
convolve (const MArray<T>& a, const MArray<R>& b,
          convn_type ct)
{
  if (a.is_empty () || b.is_empty ())
    return MArray<T> ();

  int nd = std::max (a.ndims (), b.ndims ());
  const dim_vector adims = a.dims ().redim (nd), bdims = b.dims ().redim (nd);
  dim_vector cdims = dim_vector::alloc (nd);

  for (int i = 0; i < nd; i++)
    {
      if (ct == convn_valid)
        cdims(i) = std::max (adims(i) - bdims(i) + 1,
                             static_cast<octave_idx_type> (0));
      else
        cdims(i) = std::max (adims(i) + bdims(i) - 1,
                             static_cast<octave_idx_type> (0));
    }

  MArray<T> c (cdims, T());

  convolve_nd<T, R> (a.fortran_vec (), adims, adims.cumulative (),
                     b.fortran_vec (), bdims, bdims.cumulative (),
                     c.fortran_vec (), cdims.cumulative (), nd, ct == convn_valid);

  if (ct == convn_same)
    {
      // Pick the relevant part.
      Array<idx_vector> sidx (dim_vector (nd, 1));

      for (int i = 0; i < nd; i++)
        sidx(i) = idx_vector::make_range (bdims(i)/2, 1, adims(i));
      c = c.index (sidx);
    }

  return c;
}


//template <class T> static void encode_uint_image(std::vector<Magick::Image>& imvec, const octave_value& img, bool has_map);
namespace fs = boost::filesystem;
octave_value ov_sub2ind(octave_value_list& args);
octave_value ov_linspace(octave_value_list args);
Matrix round_col_vec(Matrix& m);
void read_mat_ascii_data (std::istream& is, const std::string& filename, octave_value& tc);

struct Localizations
{
    std::string  identifier;
    std::string  syntax;
    std::string  semantic;
    std::string  unit;
    std::string  min;
    std::string  max;
};

typedef std::vector<Localizations> loc;
std::vector<double> readMaxValues( std::string xmlConfigString);

//template <class T> void encode_uint_image (std::vector<Magick::Image>& imvec,// const NDArray& img);
       //                                             const octave_value& img);

//void write_uintimage (std::string filename, const octave_value& img); //const NDArray& img);

//template <class T> octave_value_list read_images (std::vector<Magick::Image>& imvec, unsigned int depth);

int AnalyseDirs(boost::property_tree::ptree pt, fs::path config_file_full_path, bool force2D, bool one_color, std::vector<std::string> inputdirs,
                int skip_stage, std::string  working_dir, int end_stage);
int AnalyseFiles(boost::property_tree::ptree pt, fs::path config_file_full_path, bool force2D, bool one_color, std::vector<std::string> inputfiles,
                 int skip_stage, std::string  working_dir, int end_stage, std::vector<std::string> inputdirs, int& nr_of_dir);

int PairFinderMain(boost::property_tree::ptree *pt, fs::path *config_file_full_path,  std::string  inputfile);
int FilterMain(boost::property_tree::ptree *pt, fs::path *config_file_full_path,  std::string  inputfile, std::string  working_dir,
               std::string raw_input_filename, Matrix& FullData, bool one_color);

int ReconstructorMain(boost::property_tree::ptree *pt, fs::path *config_file_full_path, std::string inputfile,
                      std::string working_dir, std::string raw_input_filename, bool force2D, Matrix& FullData);

extern int WriteTiff(octave_value& index, FloatNDArray& kernel, std::string raw_input_filename, dim_vector dv, int max_x, int max_y, int max_z, int channel_nr,double histo_correct[]);
int ShowIntensitiesMain(Matrix& xInt, Matrix& yInt, Matrix Filter[], int nr_of_filters, std::string inputfile, int max_x, int max_y);

FloatNDArray gaussFilter(double sigma_xy);
FloatNDArray gaussFilter3D(double sigma_xy, double sigma_z, double z_scale);

std::string get_date(void);
int ReturnDimensions( const std::string & pattern, const std::string & firstLine ) ;

int appendMatrix(Matrix& m, Matrix& subMatrix);
int statistics1D(Matrix& row_vec, Matrix& statisticsOut, std::string& formatedOutput);
int CentroidPoints(Matrix& coordinates, Matrix& intensities, Matrix& centroid);
int CalculateFrequency(Matrix& raw_data, Matrix& Frequency);
void deleteAllRows(Matrix& m);
void deleteRow(Matrix& m, const int& row);
Matrix vector2matrix(std::vector<double>& vec);
void deleteRowsThreshold(Matrix& m, int column, double min);
void deleteRowsMaxThreshold(Matrix& m, int column, double max);


//void exit_sdmixer(void);
#endif // __SDMIXER_H_INCLUDED__

