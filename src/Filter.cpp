
#include "sdmixer.h"
/*


template <class T> octave_value_list read_images (std::vector<Magick::Image>& imvec, unsigned int depth)
{
  typedef typename T::element_type P;

  octave_value_list retval (3, Matrix ());

  T im;

  int rows = imvec[0].baseRows ();
  int columns = imvec[0].baseColumns ();

  dim_vector idim = dim_vector ();
  idim.resize (4);
  idim(0) = rows;
  idim(1) = columns;
  idim(2) = 1;
  idim(3) = 1;


  const int divisor = ((uint64_t (1) << 8) - 1) /
                      ((uint64_t (1) << depth) - 1);


        im = T (idim);
        P *vec = im.fortran_vec ();
           const Magick::PixelPacket *pix
              = imvec[0].getConstPixels (0, 0, columns, rows);
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


  retval(0) = im;

  return retval;
}
*/

int FilterMain(boost::property_tree::ptree *pt, fs::path *config_file_full_path, std::string inputfile, std::string working_dir,
               std::string raw_input_filename, Matrix& FullData, bool one_color)
{
    //std::cout << "\n** starting filter **" << std::endl;
    // we start with the definitions from the config file
    if(!std::ifstream(inputfile.c_str()))
    {
        std::cout << "ERROR: Filter failed to open "<< inputfile << std::endl;
        return 1;
    }

    std::string str_xOffset;
    std::string str_yOffset;

    std::string str_max_xInt;
    std::string str_max_yInt;
    std::string str_precision;
    std::string str_nr_of_filters;
    std::string path_to_filters;
    std::string filter_file_pattern;
    std::string filter_FileStamp;

    try{
        str_xOffset = pt->get<std::string>("PairFinder.xOffset");
        str_yOffset = pt->get<std::string>("PairFinder.yOffset");

        str_max_xInt = pt->get<std::string>("Filter.max_xInt");
        str_max_yInt = pt->get<std::string>("Filter.max_yInt");
        str_precision = pt->get<std::string>("Filter.precision");
        str_nr_of_filters = pt->get<std::string>("Filter.nr_of_filters");
        path_to_filters = pt->get<std::string>("Filter.path_to_filters");
        filter_file_pattern = pt->get<std::string>("Filter.filter_file_pattern");
        filter_FileStamp = pt->get<std::string>("Filter.FileStamp");
    }
    catch( const std::exception & ex ){

        std::cout << "\nERROR: Filter found syntax error in config file: "  << std::endl << config_file_full_path->string() << std::endl << ex.what() << std::endl;

        return 1;

    }
    double xOffset = atof(str_xOffset.c_str());
    double yOffset = atof(str_yOffset.c_str());
    //std::cout << xOffset << "  "  << yOffset << std::endl;

    int max_xInt = atoi(str_max_xInt.c_str());
    int max_yInt = atoi(str_max_yInt.c_str());
    double precision = atof(str_precision.c_str());
    int nr_of_filters = 1;

    if (!one_color)
    {
      nr_of_filters  = atoi(str_nr_of_filters.c_str());
    }
    int dimensions, dyeColOctave, xColIntOctave, yColIntOctave, full_size;

    int counter_crosstalk = 0;
    octave_value_list imread, rot90, result_imread, result_rot90, load, sort, result;

    octave_value_list warn;
    warn(0) = "off" ;
    result = feval("warning", warn, 1);     // turns off this annoying GraphicsMagick warning

    Matrix FilterIm[nr_of_filters], Filter_raw[nr_of_filters];

    if (!one_color)
    {
      nr_of_filters  = atoi(str_nr_of_filters.c_str());

    // we need this pattern matching to find our filter files.

        boost::regex pattern_filter(filter_file_pattern.c_str(), boost::regex_constants::icase);
        fs::path full_path( fs::initial_path<fs::path>() );

        if (path_to_filters.empty())            // no path to filters specified in the config.. search in working dir
        {
            full_path = fs::system_complete( fs::path( working_dir.c_str() ) );
        }
        else                                    // path to filters specified in the config.. fuck the working dir
        {
            full_path = fs::system_complete( fs::path( path_to_filters.c_str() ) );
        }
        std::vector<std::string> filter_list;

        fs::directory_iterator end_iter;
        int file_match_count = 0;


        for ( fs::directory_iterator dir_itr( full_path ); dir_itr != end_iter; ++dir_itr ) //iterate the  path to the filters
        {
            if ( fs::is_regular_file( dir_itr->status() ) )                 // found a regular file, could this be a filter?
            {

                std::stringstream ssfile;
                ssfile << dir_itr->path().filename();

                if (boost::regex_match(ssfile.str(), pattern_filter))           // yes! it matches the filter pattern from config.file
                {

                    std::string match_file = ssfile.str();

                    fs::path full_path_matched_file = operator/(full_path.string(), dir_itr->path().filename());
                    filter_list.push_back(full_path_matched_file.string());

                    ++file_match_count;
                }
            }

        }
        if (file_match_count == 0)
        {
            std::cout << "ERROR: found no filter files in " << full_path.string() << "\ncheck directory and regex pattern in config file!" << std::endl;
            return 1;
        }
        if(nr_of_filters < file_match_count)
        {
            std::cout << "WARNING: found " << file_match_count << " filter files in " << full_path.string() << "\nbut in configfile nr_of_filters is set to " <<
            nr_of_filters << "... only "  << nr_of_filters << " filter(s) will be used.\n" << std::endl;
        }
        if(nr_of_filters > file_match_count)
        {
            std::cout << "ERROR: found not enough filter files in " << full_path.string() <<
            "\nin configfile nr_of_filters is set to " <<
            nr_of_filters << "... only "  << file_match_count << " filter(s) were found." << std::endl;
        return 1;
        }

        std::sort(filter_list.begin(), filter_list.end());              // we want to apply the filters starting from 1


        for (int i = 0; i < nr_of_filters; ++i)
        {
           /* std::vector<Magick::Image> imvec;
            Magick::Image img;

            img.read(filter_list[i]);
            imvec.push_back(img);
            imread = read_images<class Matrix>(imvec,(unsigned int) 8);
            Filter_raw[i] = imread(0).matrix_value();
    */
           // Filter_raw[i]=imread(0).matrix_value();
            //std::cout << filter_list[i].c_str() << std::endl;
            Magick::InitializeMagick(NULL);

            std::vector<Magick::Image> imvec;
            Magick::readImages (&imvec, filter_list[i].c_str());
            int nframes = imvec.size ();
            Array<int> frameidx;
            frameidx = Array<int> (dim_vector (1, nframes));
            for (int j = 0; j < frameidx.length (); j++)
                frameidx(j) = j;

            imread = read_images<uint8NDArray> (imvec, frameidx, 8);
            Filter_raw[i] = imread(0).matrix_value();

            /*imread(0) = filter_list[i].c_str();
            result_imread = feval("imread", imread, 1);                 //load the filters
            Filter_raw[i] = result_imread(0).matrix_value();*/

            //std::cout << Filter_raw[i].cols() << "  " << Filter_raw[i].rows() << std::endl;

            rot90(0) = imread(0).matrix_value();
           // rot90(0) = imread(0);
            rot90(1) = 3;
            result_rot90 = feval("rot90", rot90, 2);                    // rotate three times = 270 °
            FilterIm[i] = result_rot90(0).matrix_value();

            /*rot90(1) = 2;
            result_rot90 = feval("rot90", rot90, 2);                    // rotate
            Filter_raw[i] = result_rot90(0).matrix_value();*/
            if( FilterIm[i].cols() != max_yInt * precision && FilterIm[i].rows() != max_xInt * precision)
            {
                std::cout << "ERROR: dimensions of filter file " << filter_list[i].c_str() << " ( " << FilterIm[i].rows() << " x " << FilterIm[i].cols()
                << " ) do not match definitions (max_yInt, max_xInt, precision) in config file (" <<
                max_xInt * precision << " x " << max_yInt * precision << " )"
                << std::endl;
                return 1;
            }
        }
    }
    else            //one-color switch set
    {

        for (int i = 0; i < nr_of_filters; ++i)
        {
            Filter_raw[i].resize(max_yInt * precision, max_xInt * precision, 0);
            //std::cout << Filter_raw[i].cols() << "  " << Filter_raw[i].rows() << std::endl;
            //Filter_raw[i].fill(0);
            rot90(0) = octave_value(Filter_raw[i]);
           // rot90(0) = imread(0);
            rot90(1) = 3;
            result_rot90 = feval("rot90", rot90, 2);                    // rotate three times = 270 °
            FilterIm[i] = result_rot90(0).matrix_value();
        }
    }
    long int counter_ch[nr_of_filters];

    for (int i = 0; i<nr_of_filters; ++i)
    {
        counter_ch[i] = 0;                                              //initialize
    }


    //load(0) = inputfile.c_str();                        // load file.. _pairs.out
    //result = feval("load", load, 1);

    //FullData = result(0).matrix_value();

    int totalColNr = FullData.cols();

    full_size = FullData.rows();
    //dimensions = (FullData.cols() - 6)/2;               // force2D is futile here

    std::ifstream ifs(raw_input_filename.c_str());
    std::string firstLine;
    getline (ifs, firstLine);                           // load first line (header) and

    dimensions = ReturnDimensions("Position", firstLine);       //count the occurence of "Position" --> this is the nr of dimensions!
    //std::cout << dimensions << std::endl;
    /*
    dyeColOctave = 2*dimensions + 7;
    xColIntOctave = 2*dimensions + 5;                       // 11  x y z frame yint dint x y z frame xint dint
    yColIntOctave = dimensions + 2;
    */

    dyeColOctave  = totalColNr + 1;
    if ( xOffset > yOffset)
    {
        //std::cout << "1" << std::endl;
        xColIntOctave = (totalColNr/2) + dimensions + 2;
        yColIntOctave = dimensions + 2;
    }
    else
    {
        //std::cout << "2" << std::endl;
        yColIntOctave = (totalColNr/2) + dimensions + 2;
        xColIntOctave = dimensions + 2;
    }



    logFile << "Max. xInt/yInt from config:\t" << max_xInt << " x " << max_yInt << " " << std::endl;
    logFile << "xInt/yInt binning (precision):\t" << precision << std::endl;
    logFile << "Max. xInt/yInt from File:\t" << FullData.column(xColIntOctave-1).max()  << " x " << FullData.column(yColIntOctave-1).max() << std::endl;
    logFile << "Filter image size from config:\t" << max_xInt * precision << " x " << max_yInt * precision << " px" << std::endl;



    // Apply x/y cut-off to FullData
    sort(0) = octave_value(FullData);
    sort(1) = xColIntOctave;

    result = feval("sortrows", sort, 2);

    FullData = result(0).matrix_value();
                                                        // the next one is tricky, it finds all entries in FullData, xColumns wich are > max_xInt
                                                        /* in Octave it looks like: max_xInt_coord=find(FullData(:, xcolInt) > max_xInt);*/
    result = feval ("find", do_binary_op (octave_value::op_gt, FullData.column(xColIntOctave - 1), max_xInt), 1);

    Matrix max_x_Int_coord = result(0).matrix_value();

    if (!max_x_Int_coord.rows() == 0 )                  // are the entries that match the find condition?
    {
            for (int i = max_x_Int_coord(0, 0); i < full_size + 1; ++i)
            {
                FullData(i - 1, xColIntOctave - 1) = max_xInt;          // set them to max_xInt
            }
    }

    sort(0) = octave_value(FullData);
    sort(1) = yColIntOctave;
    result = feval("sortrows", sort, 2);
    FullData = result(0).matrix_value();
                                                    // do the same to max_yInt

    result = feval ("find", do_binary_op (octave_value::op_gt, FullData.column(yColIntOctave - 1), max_yInt), 1);
    Matrix max_y_Int_coord = result(0).matrix_value();

    if (!max_y_Int_coord.rows() == 0 )
    {
            for (int i = max_y_Int_coord(0, 0); i < full_size + 1; ++i)
            {
                FullData(i - 1, yColIntOctave - 1) = max_yInt;
            }
    }

/*
    if(FullData.column(xColIntOctave-1).max() < max_xInt ||  FullData.column(yColIntOctave-1).max() < max_yInt)
    {
        std::cout << "\nWARNING: max_xInt (" << max_xInt << ") or max_yInt (" << max_yInt << ") from config file are larger" <<
        " than your actual Intensity limits: " << FullData.column(xColIntOctave-1).max() << " and " <<
        FullData.column(yColIntOctave-1).max() << " your Filter(s) will be cropped!" << std::endl;

        logFile << "Filter(s) cropped to:\t\t" << ceil(FullData.column(xColIntOctave-1).max()*precision) << " x "
        << ceil(FullData.column(yColIntOctave-1).max()*precision) << " px" << std::endl;

        for(int i = 0; i < nr_of_filters; ++i)
        {
            //std::cout << round(FullData.column(yColIntOctave-1).max()*precision) + 1 << std::endl;
            //FilterIm[i].resize(round(FullData.column(xColIntOctave-1).max()*precision), round(FullData.column(yColIntOctave-1).max()*precision));
            Filter_raw[i].resize(ceil(FullData.column(yColIntOctave-1).max()*precision), ceil(FullData.column(xColIntOctave-1).max()*precision));
        }
    }
*/
    //std::cout << "Filter initialized\n" << std::endl;


    octave_value_list round_x, round_y, result_round_x, result_round_y;
    Matrix mx_round_x, mx_round_y;

    round_x(0) = (FullData.column(xColIntOctave-1) * precision);      // now we round to int, and multiply with precision
    result_round_x = feval("ceil", round_x, 1);                    // (usually that means division by 10)
    mx_round_x = result_round_x(0).matrix_value();

    round_y(0) = (FullData.column(yColIntOctave-1) * precision);      // first we had Intensity range (>80000x>30000)
    result_round_y = feval("ceil", round_y, 1);                    // and now nice litte pics with (8000x3000)
    mx_round_y = result_round_y(0).matrix_value();

    //std::cout << "  " << std::endl;

    //std::cin.sync();

    if (ShowIntensitiesMain(mx_round_x, mx_round_y, Filter_raw, nr_of_filters, inputfile, max_xInt*precision, max_yInt*precision) == 1)
    {
        std::cout << "\nERROR: ShowIntensitiesMain failed" << std::endl;
        return 1;
    }




    int temp;

    /*Matrix dyeCol(full_size, 1);                                    // redimension the FullData matrix, because we want to add an other
    dyeCol.fill(0);                                                 // column for the Color

    octave_value_list horzcat;                                      // This could be done much much easier with FullData.reshape(..) !!
    horzcat(0) = (octave_value) FullData.extract(0, 0, full_size - 1, yColIntOctave - 2);
    horzcat(1) = result_round_y(0).matrix_value();
    horzcat(2) = (octave_value) FullData.extract(0, yColIntOctave, full_size - 1, xColIntOctave - 2);
    horzcat(3) = result_round_x(0).matrix_value();
    horzcat(4) = (octave_value) FullData.extract(0, xColIntOctave, full_size - 1, xColIntOctave);
    horzcat(5) = (octave_value) dyeCol;*/

    FullData.resize(full_size, totalColNr + 1, 0);

    FullData.insert(result_round_y(0).matrix_value().column(0),0,yColIntOctave - 1);
    FullData.insert(result_round_x(0).matrix_value().column(0),0,xColIntOctave - 1);

    //result = feval("horzcat", horzcat, 6);                          // redimension here. no we have 13 columns

    //Matrix FullData_rounded = FullData;

    for (int i = 0; i < full_size; ++i)                             // could this be done better?
    {
        temp = 0;                                                   // we loop through the whole data set and look pixel wise
                                                                    // if the pixel is in a certain filter
        for (int j = 0; j < nr_of_filters; ++j)
        {

            //if(FilterIm[j](FullData_rounded(i, xColIntOctave - 1) -1, FullData_rounded(i, yColIntOctave - 1) -1) == 0)
           // std::cout << FullData(i, xColIntOctave - 1) -1 << "  " << FullData(i, yColIntOctave - 1) -1 << std::endl;
            if(FilterIm[j](FullData(i, xColIntOctave - 1) -1, FullData(i, yColIntOctave - 1) -1) == 0)
            {
                FullData(i, dyeColOctave - 1) = j+1;

                counter_ch[j]++;
                temp++;                                             // jeah this pixel was in a certain channel
            }
        }

        if (temp == 0)                                              // no, this pixel was not assigned to a certain channel
        {
            counter_crosstalk++;
        }

    }


    std::ofstream ofsFile;
    ofsFile.open(inputfile.append(filter_FileStamp).c_str());      // write outputfile _pairs.out_filter.out

    ofsFile << FullData;
    ofsFile.close();

    logFile << "\n" << full_size << " pairs total." << std::endl;

    for( int i=0; i < nr_of_filters; ++i)
    {
         logFile << "\nfound " << counter_ch[i] << " pairs in channel " << i+1 << std::endl;
    }

    logFile << counter_crosstalk << " pairs were sorted out (crosstalk = " << ((double)counter_crosstalk/(double)full_size) *100 << " %)" << std::endl;


    return 0;

}
