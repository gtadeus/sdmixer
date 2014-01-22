#include "sdmixer.h"

// Helper functions when reading from ascii files.

// These function take care of CR/LF issues when files are opened in
// text-mode for reading.

// Skip characters from stream IS until a newline is reached.
// Depending on KEEP_NEWLINE, either eat newline from stream or
// keep it unread.

template <class MT>
static octave_value do_linspace (const octave_value& base, const octave_value& limit,
                                octave_idx_type n)
{
  typedef typename MT::column_vector_type CVT;
  typedef typename MT::element_type T;

  octave_value retval;

  if (base.is_scalar_type ())
    {
      T bs = octave_value_extract<T> (base);
      if (limit.is_scalar_type ())
        {
          T ls = octave_value_extract<T> (limit);
          retval = linspace (bs, ls, n);
        }
      else
        {
          CVT lv = octave_value_extract<CVT> (limit);
          CVT bv (lv.length (), bs);
          retval = linspace (bv, lv, n);
        }
    }
  else
    {
      CVT bv = octave_value_extract<CVT> (base);
      if (limit.is_scalar_type ())
        {
          T ls = octave_value_extract<T> (limit);
          CVT lv (bv.length (), ls);
          retval = linspace (bv, lv, n);
        }
      else
        {
          CVT lv = octave_value_extract<CVT> (limit);
          retval = linspace (bv, lv, n);
        }
    }

  return retval;
}
octave_value ov_linspace(octave_value_list args)
{
  octave_value retval;

  int nargin = args.length ();

  octave_idx_type npoints = 100;

  if (nargin == 3)
    npoints = args(2).idx_type_value ();

    octave_value arg_1 = args(0);
    octave_value arg_2 = args(1);

    retval = do_linspace<Matrix> (arg_1, arg_2, npoints);


  return retval;
}

static dim_vector get_dim_vector (const octave_value& val, const char *name)
{
  RowVector dimsv = val.row_vector_value (false, true);
  dim_vector dv;
  octave_idx_type n = dimsv.length ();

  if (n < 1)
    error ("%s: dimension vector DIMS must not be empty", name);
  else
    {
      dv.resize (std::max (n, static_cast<octave_idx_type> (2)));
      dv(1) = 1;
      for (octave_idx_type i = 0; i < n; i++)
        {
          octave_idx_type ii = dimsv(i);
          if (ii == dimsv(i) && ii >= 0)
            dv(i) = ii;
          else
            {
              error ("%s: dimension vector DIMS must contain integers", name);
              break;
            }
        }
    }

  return dv;
}

octave_value ov_sub2ind(octave_value_list& args)
{
  int nargin = args.length ();
  //std::cout << nargin << std::endl;
octave_value retval;

      dim_vector dv = get_dim_vector (args(0), "sub2ind");
      Array<idx_vector> idxa (dim_vector (nargin-1, 1));

          dv = dv.redim (nargin -1);
          for (int j = 0; j < nargin -1 ; j++)
            {
//std::cout << idxa.as_matrix()<< std::endl;
                  idxa(j) = args(j+1).index_vector ();
                    if (j > 0 && args(j+1).dims () != args(1).dims ())
                    std::cout << "\nERROR: sub2ind: all subscripts must be of the same size" << std::endl;


            }

          idx_vector idx = sub2ind (dv, idxa);
          retval = idx;

  return retval;
}

Matrix round_col_vec(Matrix& m)
{
    Matrix return_val(m.rows(), m.cols());
    for (int i = 0; i < m.rows(); ++i)
    {
        return_val(i,0) = round(m(i,0));
    }
    return return_val;
}

void skip_until_newline (std::istream& is, bool keep_newline)
{
  if (! is)
    return;

  while (is)
    {
      char c = is.peek ();

      if (c == '\n' || c == '\r')
        {
          // Reached newline.
          if (! keep_newline)
            {
              // Eat the CR or LF character.
              char d;
              is.get (d);

              // Make sure that for binary-mode opened ascii files
              // containing CRLF line endings we skip the LF after CR.
              if (c == '\r' && is.peek () == '\n')
                {
                  // Yes, LF following CR, eat it.
                  is.get (d);
                }
            }

          // Newline was found, and read from stream if
          // keep_newline == true, so exit loop.
          break;
        }
      else
        {
          // No newline charater peeked, so read it and proceed to next
          // character.
          char d;
          is.get (d);
        }
    }
}


// If stream IS currently points to a newline (a leftover from a
// previous read) then eat newline(s) until a non-newline character is
// found.

void
skip_preceeding_newline (std::istream& is)
{
  if (! is)
    return;

  // Check whether IS currently points to newline character.
  char c = is.peek ();

  if (c == '\n' || c == '\r')
    {
      // Yes, at newline.
      do
        {
          // Eat the CR or LF character.
          char d;
          is.get (d);

          // Make sure that for binary-mode opened ascii files
          // containing CRLF line endings we skip the LF after CR.
          if (c == '\r' && is.peek () == '\n')
            {
              // Yes, LF following CR, eat it.
              is.get (d);
          }

          // Peek into next character.
          c = is.peek ();

          // Loop while still a newline ahead.
        }
      while (c == '\n' || c == '\r');
    }
}

// Read charaters from stream IS until a newline is reached.
// Depending on KEEP_NEWLINE, either eat newline from stream or keep
// it unread.  Characters read are stored and returned as
// std::string.

std::string
read_until_newline (std::istream& is, bool keep_newline)
{
  if (! is)
    return std::string ();

  std::ostringstream buf;

  while (is)
    {
      char c = is.peek ();

      if (c == '\n' || c == '\r')
        {
          // Reached newline.
          if (! keep_newline)
            {
              // Eat the CR or LF character.
              char d;
              is.get (d);

              // Make sure that for binary-mode opened ascii files
              // containing CRLF line endings we skip the LF after
              // CR.

              if (c == '\r' && is.peek () == '\n')
                {
                  // Yes, LF following CR, eat it.
                  is.get (d);
                }
            }

          // Newline was found, and read from stream if
          // keep_newline == true, so exit loop.
          break;
        }
      else
        {
          // No newline charater peeked, so read it, store it, and
          // proceed to next.
          char d;
          is.get (d);
          buf << d;
        }
    }

  return buf.str ();
}

static std::string get_mat_data_input_line (std::istream& is)
{
  std::string retval;

  bool have_data = false;

  do
    {
      retval = "";

      char c;
      while (is.get (c))
        {
          if (c == '\n' || c == '\r')
            {
              is.putback (c);
              skip_preceeding_newline (is);
              break;
            }

          if (c == '%' || c == '#')
            {
              skip_until_newline (is, false);
              break;
            }

          if (! is.eof ())
            {
              if (! have_data && c != ' ' && c != '\t')
                have_data = true;

              retval += c;
            }
        }
    }
  while (! (have_data || is.eof ()));

  return retval;
}
static void get_lines_and_columns (std::istream& is, const std::string& filename, octave_idx_type& nr, octave_idx_type& nc)
{
  std::streampos pos = is.tellg ();

  int file_line_number = 0;

  nr = 0;
  nc = 0;

  while (is && ! error_state)
    {
      octave_quit ();

      std::string buf = get_mat_data_input_line (is);

      file_line_number++;

      size_t beg = buf.find_first_not_of (", \t");

      // If we see a CR as the last character in the buffer, we had a
      // CRLF pair as the line separator.  Any other CR in the text
      // will not be considered as whitespace.

      if (beg != std::string::npos && buf[beg] == '\r' && beg == buf.length () - 1)
        {
          // We had a blank line ending with a CRLF.  Handle it the
          // same as an empty line.
          beg = std::string::npos;
        }

      octave_idx_type tmp_nc = 0;

      while (beg != std::string::npos)
        {
          tmp_nc++;

          size_t end = buf.find_first_of (", \t", beg);

          if (end != std::string::npos)
            {
              beg = buf.find_first_not_of (", \t", end);

              if (beg == std::string::npos || (buf[beg] == '\r' &&
                                  beg == buf.length () - 1))
                {
                  // We had a line with trailing spaces and
                  // ending with a CRLF, so this should look like EOL,
                  // not a new colum.
                  break;
                }
            }
          else
            break;
        }

      if (tmp_nc > 0)
        {
          if (nc == 0)
            {
              nc = tmp_nc;
              nr++;
            }
          else if (nc == tmp_nc)
            nr++;
          else
            error ("load: %s: inconsistent number of columns near line %d",
                   filename.c_str (), file_line_number);
        }
    }

  if (nr == 0 || nc == 0)
    error ("load: file `%s' seems to be empty!", filename.c_str ());

  is.clear ();
  is.seekg (pos);
}

// Extract a matrix from a file of numbers only.
//
// Comments are not allowed.  The file should only have numeric values.
//
// Reads the file twice.  Once to find the number of rows and columns,
// and once to extract the matrix.
//
// FILENAME is used for error messages.
//
// This format provides no way to tag the data as global.

void read_mat_ascii_data (std::istream& is, const std::string& filename, octave_value& tc)
{


      octave_idx_type nr = 0;
      octave_idx_type nc = 0;

      int total_count = 0;

      get_lines_and_columns (is, filename, nr, nc);

      octave_quit ();

      if (! error_state && nr > 0 && nc > 0)
        {
          Matrix tmp (nr, nc);

          if (nr < 1 || nc < 1)
            is.clear (std::ios::badbit);
          else
            {
              double d;
              for (octave_idx_type i = 0; i < nr; i++)
                {
                  std::string buf = get_mat_data_input_line (is);

                  std::istringstream tmp_stream (buf);

                  for (octave_idx_type j = 0; j < nc; j++)
                    {
                      octave_quit ();

                      d = octave_read_value<double> (tmp_stream);

                      if (tmp_stream || tmp_stream.eof ())
                        {
                          tmp.elem (i, j) = d;
                          total_count++;

                          // Skip whitespace and commas.
                          char c;
                          while (1)
                            {
                              tmp_stream >> c;

                              if (! tmp_stream)
                                break;

                              if (! (c == ' ' || c == '\t' || c == ','))
                                {
                                  tmp_stream.putback (c);
                                  break;
                                }
                            }

                          if (tmp_stream.eof ())
                            break;
                        }
                      else
                        {
                          error ("load: failed to read matrix from file `%s'",
                                 filename.c_str ());

                         // return retval;
                        }

                    }
                }
            }

          if (is || is.eof ())
            {
              // FIXME -- not sure this is best, but it works.

              if (is.eof ())
                is.clear ();

              octave_idx_type expected = nr * nc;

              if (expected == total_count)
                {
                  tc = tmp;
                  //retval = varname;
                }
              else
                error ("load: expected %d elements, found %d",
                       expected, total_count);
            }
          else
            error ("load: failed to read matrix from file `%s'",
                   filename.c_str ());
        }
      else
        error ("load: unable to extract matrix size from file `%s'",
               filename.c_str ());


  //return retval;
}

/*
void write_uintimage (std::string filename, const octave_value& img)
                      //const NDArray& img)
{
    Magick::InitializeMagick("");
    //std::cout << "writing " << filename << std::endl;
    std::vector<Magick::Image> imvec;
    Matrix m;
    m = img.uint8_array_value();

    //encode_uint_image<uint8NDArray> (imvec, img);
    encode_uint_image<uint8NDArray> (imvec, m);

    try
    {
      Magick::writeImages (imvec.begin (), imvec.end (), filename);
    }
    catch (Magick::WarningCoder& w)
    {
        std::cout << "Magick++ CoderWarning: " << w.what () << std::endl;
    }
    catch (Magick::Warning& w)
    {
        std::cout << "Magick++ warning: " << w.what () << std::endl;
    }
    catch (Magick::ErrorCoder& e)
    {
        std::cout << "Magick++ coder error: " << e.what () << std::endl;
    }
    catch (Magick::Exception& e)
    {
        std::cout << "Magick++ exception: " << e.what () << std::endl;
    }

}

template <class T> void encode_uint_image (std::vector<Magick::Image>& imvec,// const NDArray& img)
                                                 const octave_value& img)
{
  unsigned int bitdepth = 0;
  T m;

 //if (img.is_uint8_type ())
 //   {
      bitdepth = 8;
      m = img.uint8_array_value ();
 //  }
 //   else{
 //   std::cout << "\nERROR: encode_uint_image, unsupported output format" << std::endl;
//    }

   // m = img;
  dim_vector dsizes = m.dims ();
  unsigned int nframes = 1;
  if (dsizes.length () == 4)
    nframes = dsizes(4);

  Array<octave_idx_type> idx (dim_vector (dsizes.length (), 1));
  octave_idx_type rows = m.rows ();
  octave_idx_type columns = m.columns ();

  unsigned int div_factor = (1 << bitdepth) - 1;

  for (unsigned int ii = 0; ii < nframes; ii++)
    {
      Magick::Image im (Magick::Geometry (columns, rows), "black");

      im.depth (bitdepth);

      //im.classType (Magick::DirectClass);
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
                      idx(3) = ii;
                    }

                  c.shade (static_cast<double>(m(idx)) / div_factor);

                  im.pixelColor (y, x, c);
                }
            }

          im.quantizeColorSpace (Magick::GRAYColorspace);
          im.quantizeColors (1 << bitdepth);
          im.quantize ();


      imvec.push_back (im);
    }
}
*/
