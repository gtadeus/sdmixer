
#include "sdmixer.h"

#include <octave/lo-mappers.h>

int ReconstructorMain(boost::property_tree::ptree *pt, fs::path *config_file_full_path, std::string inputfile,
                      std::string working_dir, std::string raw_input_filename, bool force2D, Matrix& FullData)
{
    //std::cout << "*** Reconstructor *** "  << std::endl;
    // lets start with the definitions
    if(!std::ifstream(inputfile.c_str()))
    {
        std::cout << "ERROR: Reconstructor failed to open "<< inputfile << std::endl;
        return 1;
    }
    std::string str_xy_binning;
    std::string str_z_binning;
    std::string str_FWHM_xy;
    std::string str_FWHM_z;
    std::string str_NM_PER_PX;
    std::string str_histo_correct;
    std::string str_histo_threshold;
    std::string str_sqrt_cummulation;

    try{


        str_NM_PER_PX = pt->get<std::string>("General.NM_PER_PX");
        str_xy_binning = pt->get<std::string>("Reconstructor.xy_binning");
        str_z_binning = pt->get<std::string>("Reconstructor.z_binning");
        str_FWHM_xy = pt->get<std::string>("Reconstructor.FWHM_xy");
        str_FWHM_z = pt->get<std::string>("Reconstructor.FWHM_z");
        str_histo_correct = pt->get<std::string>("Reconstructor.histo_correct");
    }
    catch( const std::exception & ex ){

        std::cout << "\nERROR: Reconstructor found syntax error in config file: "  << std::endl << config_file_full_path->string() << std::endl << ex.what() << std::endl;

        return 1;

    }
    double histo_correct[3];
    try{

        str_histo_threshold = pt->get<std::string>("Reconstructor.histo_threshold");
        str_sqrt_cummulation = pt->get<std::string>("Reconstructor.sqrt_cummulation");
        histo_correct[1] = atof(str_histo_threshold.c_str());
        histo_correct[2] = atof(str_sqrt_cummulation.c_str());
    }
    catch( const std::exception & ex ){

        logFile << "did not find optional histogram equalisation parameter (histo_threshold, sqrt_cummulation) in config file, loading default\n" <<std::endl;
        histo_correct[1] = 1;
        histo_correct[2] = 0;
    }

    int NM_PER_PX = atoi(str_NM_PER_PX.c_str());
    double xy_binning = atof(str_xy_binning.c_str());
    double z_binning = atof(str_z_binning.c_str());
    double sigma_xy = ((atof(str_FWHM_xy.c_str())) / ( 2 * sqrt(2*log(2)))) / ((double) NM_PER_PX/xy_binning);
    double sigma_z = ((atof(str_FWHM_z.c_str())) / ( 2 * sqrt(2*log(2)))) / ((double) NM_PER_PX/xy_binning);

    histo_correct[0] = atof(str_histo_correct.c_str());

//    double threshold = atof(str_threshold.c_str());

    logFile << "xy binning:\t\t" << xy_binning << "\nz_binning:\t\t" << z_binning << "\nFWHM xy:\t\t" << str_FWHM_xy << "\nFWHM z:\t\t\t"
    << str_FWHM_z << "\nhisto_correct:\t\t" << histo_correct[0] << "\nhisto_threshold:\t" << histo_correct[1] << "\nsqrt_cummulation:\t" << histo_correct[2]<<std::endl;

    octave_value_list load, result, sort, unique,  ind_result, histc, histc_result, imwrite, uint, conv;

    octave_value_list warn;
    warn(0) = "off" ;
    result = feval("warning", warn, 1);     // turns off this annoying GraphicsMagick warning

    int xColOctave = 2;                     // weird, isnt it? otherwise images would be rotated90
    int yColOctave = 1;                     // should we put this into configfile?
    int zColOctave = 3;


    octave_idx_type max_x, max_y;
    octave_idx_type max_z;
   // int dimensions = (FullData.cols()-7)/2;     //the real dims from file
    //int dyeColOctave = 2*dimensions + 7;
//    int frameColOctave=dimensions+1;
    std::ifstream ifs(raw_input_filename.c_str());
    std::string firstLine;
    getline (ifs, firstLine);                           // load first line (header) and
    int dimensions = ReturnDimensions("Position", firstLine);       //count the occurence of "Position"
                                                                //  --> this is the nr of dimensions!

    std::vector<double> maxValues = readMaxValues(firstLine);
    max_x = round(maxValues[0]*1e9/xy_binning);
    max_y = round(maxValues[1]*1e9/xy_binning);
    max_z = 1;
    logFile << "Max from header: " << maxValues[0] << "  " << maxValues[1];

    if (!force2D && dimensions == 3)
    {
        deleteRowsThreshold(FullData, zColOctave - 1, 0);
        //deleteRowsMaxThreshold(FullData, zColOctave - 1, maxValues[2]*1e9);
        max_z = round(maxValues[2]*1e9 / z_binning);
        logFile << "  " << maxValues[2];
    }
    logFile << std::endl;

    logFile << "image size (x, y, z): " << max_x << "  " << max_y << "  " << max_z << std::endl;


    int old = FullData.rows();
 /*   deleteRowsThreshold(FullData, xColOctave - 1, 0);       // delete all x, y < 0 for later reconstruction
    logFile << "deleted " << old-FullData.rows() << " rows. x < 0" <<std::endl;
    old = FullData.rows();
    deleteRowsThreshold(FullData, yColOctave - 1, 0);
    logFile << "deleted " << old-FullData.rows() << " rows. y < 0" <<std::endl;

    old = FullData.rows();
    deleteRowsMaxThreshold(FullData, xColOctave - 1, maxValues[0]*1e9);
    logFile << "deleted " << old-FullData.rows() << " rows. x > maxValues[0] " <<std::endl;

    old = FullData.rows();
    deleteRowsMaxThreshold(FullData, yColOctave - 1, maxValues[1]*1e9);
    logFile << "deleted " << old-FullData.rows() << " rows. y > maxValues[1]" <<std::endl;
*/

    int totalColNr = FullData.cols();
    int dyeColOctave  = totalColNr;

    int xColOctave2 = ((totalColNr-1)/2)+xColOctave;
    int yColOctave2 = ((totalColNr-1)/2)+yColOctave;
    int zColOctave2 = ((totalColNr-1)/2)+zColOctave;

    sort(0) = octave_value(FullData);
    sort(1) = dyeColOctave;

    result = feval("sortrows", sort, 2);
    FullData = result(0).matrix_value();

    unique(0) = FullData.column(dyeColOctave - 1);
    result = feval ("unique", unique, 2);           // how many dyes do we have? uniqe(dye_col)

    Matrix dyes, Nd;

    dyes = result(0).matrix_value();
    Nd = result(1).matrix_value();

    int nr_of_dyes = dyes.rows();                      // now we know how many dyes

    Matrix x, y, z;

    Matrix size_of_img2D(1,2);                           // matrix for the size of our output tiff
    Matrix size_of_img3D(1,3);                                                  // matrix for the size of our output tiff
    std::string str_ch_nr;
    FloatNDArray kernel2D;
    FloatNDArray kernel3D;                                   // there are no 3d matrixes in octave c++ API


    kernel2D = gaussFilter(sigma_xy);          // matrix gauss kernel
    kernel3D = gaussFilter3D(sigma_xy, sigma_z, z_binning/xy_binning ); // cube gauss kernel for 3D convolution

    octave_value ov_result;
    octave_value_list args;

    double img_x_origin = 0, img_y_origin = 0, img_z_origin = 0;
    /*Matrix glob_x, glob_y, glob_z;

    glob_x = FullData.extract(0, xColOctave - 1, FullData.rows() - 1, xColOctave - 1) / xy_binning;
    glob_y  = FullData.extract(0, yColOctave - 1, FullData.rows() - 1, yColOctave - 1) / xy_binning;

    //std::cout << glob_x.column_max().max() << "  " << glob_y.column_max().max()  << std::endl;
    //std::cout << glob_x.column_min().min() << "  " << glob_y.column_min().min()  << std::endl;

    if (!force2D && dimensions == 3)
    {
        glob_z = FullData.extract(0, zColOctave - 1, FullData.rows() - 1, zColOctave - 1) / z_binning;
        //img_z_origin = glob_z.column_min().min();
        //max_z = ceil(glob_z.column_max().max()) - floor(img_z_origin) + 1;
    }*/

    //img_x_origin = glob_x.column_min().min();
    //img_y_origin = glob_y.column_min().min();
    //max_x =  ceil(glob_x.column_max().max()) - floor(img_x_origin) + 1;
    //max_y =  ceil(glob_y.column_max().max()) - floor(img_y_origin) + 1;

    //std::cout << "dimensions: " << dimensions << std::endl;
    //std::cout << max_x << "  " << max_y << "  " << max_z << std::endl;


    for (int i = 0; i < nr_of_dyes; i++)
    {
       // std::cout << "i: " << i << std::endl;
       // std::cout << "dye: "<< dyes(i, 0) << std::endl;
        if (dyes(i, 0) != 0)        //do not reconstruct from cross-talk (channel 0)
        {
            if (i == 0)                 // extract data for each channel
            {
                //std::cout << "0 - " << Nd(i, 0) - 1 << std::endl;
                x = FullData.extract(0, xColOctave - 1, Nd(i, 0) - 1, xColOctave - 1) / xy_binning;
                //std::cout << x.rows() << std::endl;
                x = x.stack(FullData.extract(0, xColOctave2 - 1, Nd(i, 0) - 1, xColOctave2 - 1) / xy_binning);
               // std::cout << x.rows() << std::endl;
                y = FullData.extract(0, yColOctave - 1, Nd(i, 0) - 1, yColOctave - 1) / xy_binning;
                y = y.stack(FullData.extract(0, yColOctave2 - 1, Nd(i, 0) - 1, yColOctave2 - 1) / xy_binning);
                if (!force2D && dimensions == 3)
                {
                    z = FullData.extract(0, zColOctave - 1, Nd(i, 0) - 1, zColOctave - 1) / z_binning;
                    z = z.stack(FullData.extract(0, zColOctave2 - 1, Nd(i, 0) - 1, zColOctave2 - 1) / z_binning);
                }
            }
            else
            {
                //std::cout <<  Nd(i - 1, 0) << " - " <<  Nd(i, 0) - 1  << std::endl;
                x = FullData.extract(Nd(i - 1, 0), xColOctave - 1, Nd(i, 0) - 1, xColOctave - 1) / xy_binning;
                x = x.stack(FullData.extract(Nd(i - 1, 0), xColOctave2 - 1, Nd(i, 0) - 1, xColOctave2 - 1) / xy_binning);
                y = FullData.extract(Nd(i - 1, 0), yColOctave - 1, Nd(i, 0) - 1, yColOctave - 1) / xy_binning;
                y = y.stack(FullData.extract(Nd(i - 1, 0), yColOctave2 - 1, Nd(i, 0) - 1, yColOctave2 - 1) / xy_binning);
                if (!force2D && dimensions == 3)
                {
                    z = FullData.extract(Nd(i - 1, 0), zColOctave - 1, Nd(i, 0) - 1, zColOctave - 1) / z_binning;
                    z = z.stack(FullData.extract(Nd(i - 1, 0), zColOctave2 - 1, Nd(i, 0) - 1, zColOctave2 - 1) / z_binning);
                }
            }

            //deleteRowsThreshold(x, 0, 0);
            //deleteRowsThreshold(y, 0, 0);
            //std::cout << x.column_min().min() << "  " << y.column_min().min() << std::endl;

            //x = x - img_x_origin;       // subtract minimum to get rid of negative numbers
             //round to integer...
            //std::cout << x.column_min().min() << std::endl;
            x = round_col_vec(x);//+1;           // no zeros are allowed

            //y = y - img_y_origin;
            y = round_col_vec(y);//+1;

           //std::cout << result.matrix_value() << std::endl;
           //std::cout << dimensions << std::endl;
            if (!force2D && dimensions == 3)
            {

               // z = z - img_z_origin;
                z = round_col_vec(z); //+ 1;

                dim_vector dv(max_y, max_x, max_z);         // dimension vector

                Matrix ind = max_y * y + x + 1 + z * max_x * max_y ;
                ov_result = octave_value(ind);

                if(nr_of_dyes == 1)
                    WriteTiff(ov_result, kernel3D, raw_input_filename, dv, max_y, max_x, max_z, i + 1, histo_correct);
                else
                    WriteTiff(ov_result, kernel3D, raw_input_filename, dv, max_y, max_x, max_z, i, histo_correct);
            }

            else{
              //  size_of_img2D(0,0) =  max_y;                       //size of output image
               // size_of_img2D(1,0) =  max_x;

                Matrix ind = 1 + x + y*max_y;

                ov_result = octave_value(ind);

                dim_vector dv(max_y, max_x);

                if(nr_of_dyes == 1)
                    WriteTiff(ov_result, kernel2D, raw_input_filename, dv, max_x, max_y, max_z, i + 1, histo_correct);
                else
                    WriteTiff(ov_result, kernel2D, raw_input_filename, dv, max_x, max_y, max_z, i, histo_correct);
            }
        }
    }
    return 0;
}


int WriteTiff(octave_value& index, FloatNDArray& kernel, std::string raw_input_filename, dim_vector dv,
              int max_x, int max_y, int max_z, int channel_nr, double histo_correct[])
{
    octave_value_list sort, result_acc, result_un, unique, conv, imwrite, accumarray, hist_uniq, conv_result;

    Matrix uniq, vec_i, vec_j, N, histogram, Nhist;

    FloatNDArray img(dv);

    double histo_treshold = histo_correct[1];
    int square_eq = histo_correct[2];

    double NPixel=max_x*max_y*max_z;

    int bit = 16;
    const int mp=pow(2, bit)-1;
    double hist[mp+1];
    double PDFwt[mp+1];
    double cdf[mp+1];

    for (int i = 0; i < NPixel; ++i)
    {
        img(i) =  0;
    }

    double hist_max;

    std::string str_ch_nr;

    unique(0) = index.matrix_value();
    result_un = feval ("unique", unique, 4);
    uniq = result_un(0).matrix_value();
    vec_i = result_un(1).matrix_value();
    vec_j = result_un(2).matrix_value();

    accumarray(0) = octave_value(vec_j);
    accumarray(1) = octave_value(1);
    result_acc = feval ("accumarray", accumarray, 2);
    N = result_acc(0).matrix_value();

    hist_max = N.column_max().max();
    logFile << "hist_max = " << hist_max<<  std::endl;
   // hist_max/=20;
    //raw hist
    logFile << "uniq.rows = " << uniq.rows() << std::endl;
    logFile << "round(N(j)*((mp)/hist_max)= " << N(uniq.rows()-1)*((mp)/hist_max) << std::endl;
    logFile << "uniq(uniq.rows()-1) = " << uniq(uniq.rows()-1)<< std::endl;
    logFile << "uniq(0) = " << uniq(0)<< std::endl;
    //logFile << uniq << std::endl;


    for (int j = 0; j < uniq.rows(); ++j)
    {
        //
        if(uniq(j)>=0 && uniq(j)<=NPixel)
        {
        img(uniq(j)) = floor(N(j)*((mp)/hist_max));           //linear assignment to 16bit
        }
        else
        {
           //std::cout << uniq(j) << "    " << N(j) << std::endl;
        }

    }


if (histo_correct!=0)
{
   //std::cout << "hist after conv" << std::endl;
    for (long g=0; g < (mp+1); ++g)                             //init histogram
        hist[g]=0;

    for (int i = 0; i < NPixel; ++i)                            //fill histogram
    {
        hist[(int)img(i)]+=1;
    }

    PDFwt[0]=(double)hist[0];
    PDFwt[mp]=(double)hist[mp];

    double Pmax = hist[0]/NPixel;
    double Pl = hist[0]/NPixel;

    for (int g=0; g < mp; ++g)              //calculate normal PDF [0,1]
    {
        PDFwt[g]=hist[g]/NPixel;
        if(PDFwt[g]<Pl)
            Pl = hist[g]/NPixel;
        if(PDFwt[g]>Pmax)
            Pmax = hist[g]/NPixel;
    }

    double Pu = histo_treshold * Pmax;                 //thresholding the PDF

    for (int g=0; g < mp; ++g)              //calculate weightet PDF
    {
        PDFwt[g]=pow( (PDFwt[g]-Pl)/(Pu-Pl), histo_correct[0]) * Pu;
        PDFwt[g]*=NPixel;                   //transform to histogram range [0,NPixel]
    }

    if (square_eq==1)
    {
        cdf[0]=sqrt(PDFwt[0]);
        for (long b=1; b<(mp+1); ++b)   //from PDF to CDF with square root of hist values
        {
            cdf[b]=cdf[b-1]+sqrt(PDFwt[b]);
        }
    }
    else{
        cdf[0]=PDFwt[0];
        for (long b=1; b<(mp+1); ++b)   //from PDF to CDF, classic linear approach
        {
            cdf[b]=cdf[b-1]+PDFwt[b];
        }

    }

    for (int i = 0; i < NPixel; ++i)         //for each pixel remap
    {

        img(i)=round(  mp * (cdf[(int)img(i)]-cdf[0]) / (cdf[mp]-cdf[0]));
        //img(i)=round(  (mp+1) * (cdf[(int) img(i)]-cdf[0]) );
    }

}
/*
std::cout << "1hist after hist eq. \n\n" << std::endl;
    for (long g=0; g < (mp+1); ++g)                             //init histogram
        hist[g]=0;
std::cout << "2hist after hist eq." << std::endl;
    for (int i = 0; i < NPixel; ++i)                            //fill histogram
    {
        hist[(int)img(i)]+=1;
    }
std::cout << "3hist after hist eq." << std::endl;
    for (long g=0; g < (mp+1); ++g)                             //std:cout histo
        logFile  << hist[g] << std::endl;

*/
    clock_t init, stat;
    init=clock();
    std::cout <<"load convolution..." << std::endl;
    //std::cout << octave_value(kernel).matrix_value() << std::endl;
    //logFile <<"WARNING: convolution is off!" << std::endl;


   /* conv(0) = octave_value(img);                               //convolution
    conv(1) = octave_value(kernel);
    conv(2) = "same";
    conv_result = feval("convn", conv, 3);
    img = conv_result(0).array_value();*/
    //std::cout << octave_value(img).matrix_value().column_max().max() << std::endl;
    //NDArray conv_img = convolve(img, kernel, convn_same);
    FloatNDArray conv_img = convolve(img, kernel, convn_same);

    //NDArray uint8convimg = octave_value(conv_img).uint8_array_value();
    //std::cout << octave_value(conv_img).matrix_value().column_max().max() << std::endl;
    stat=clock()-init;
    logFile << "convolution finished" << std::endl;
    std::cout << "convolution finished, rescaling.." << std::endl;
    //logFile << "convolution executed in " << (double)stat / ((double)CLOCKS_PER_SEC) << " seconds" << std::endl;

    for (int i = 0; i < NPixel; ++i)         //back to 8bit
    {
        //img(i) = floor(img(i)*(255.0/(mp)));
        conv_img(i) = (conv_img(i)*(1.0/(mp)));
    }
    std::cout << "rescale finished" << std::endl;
    logFile << "rescale finished" << std::endl;
    /*std::ofstream txtImg;
    txtImg.open("/home/gtadeus/txtImg.out", std::ios::out);
    txtImg << octave_value(img).matrix_value() << std::endl;
    txtImg.close();

    std::ofstream txtImgConv;
    txtImgConv.open("/home/gtadeus/txtImgConv.out", std::ios::out);
    txtImgConv << octave_value(conv_img).matrix_value() << std::endl;
    txtImgConv.close();*/
   // std::cout << kernel << std::endl;



    //dim_vector dvTiff(max_x, max_y, 1, max_z);              // multilayer Tiff: (x, y, RGB, z)    RGB = 1,2,3
/*
    uint8NDArray TiffImg(dvTiff);
    std::cout << "start reshape " << std::endl;
    TiffImg = uint8convimg.reshape(dvTiff);                          // reshape is easy for gray scale, z = 1 for 2D
*/
    str_ch_nr = static_cast<std::ostringstream*>( &(std::ostringstream() << (channel_nr)) )->str(); //get out an int from string (crazy..)

    std::string output_image= raw_input_filename.append("_ch").append(str_ch_nr);//.append(".tif");
    //std::string output_image2= raw_input_filename.append("2.tif");

    std::cout << "start writing output channel nr. " << channel_nr << " (this may take some time...) " << std::endl;
    /*
    imwrite(0) = octave_value (TiffImg);
    imwrite(1) = output_image.c_str();
    result = feval("imwrite", imwrite, 3);*/


    //double * pix_p=(double *)conv_img.fortran_vec();
 /*   Magick::Image im(Magick::Geometry (conv_img.cols(), conv_img.rows()), "black");
    im.compressType(Magick::LZWCompression);
    im.type (Magick::GrayscaleType);
    im.classType (Magick::PseudoClass);
    im.depth (16);
    Magick::ColorGray c;
    dim_vector dsizes = conv_img.dims();
    Array<octave_idx_type> idx (dim_vector (dsizes.length (), 1));
    //Magick::Image im( img.cols(), img.rows(), "R", Magick::DoublePixel, pix_p );
    for (int y = 0; y < conv_img.cols(); y++)
    {
        idx(1) = y;
        for (int x=0; x < conv_img.rows(); x++)
        {
            idx(0) = x;
            //c.shade(conv_img(x, y));
            c.shade (conv_img(idx));
            im.pixelColor (y, x, c);
        }
    }
*/
   // for ( int i = 0; i < max_z ; i++)
    {
        int i = 0;
        std::vector<Magick::Image> imvec;
        //logFile << "start encode_uint_image, frame " << i + 1 << std::endl;
       // if( (i % 18) == 0)
        //    std::cout << "writing frame " << i + 1 <<  " (" << ((float)i / (float)max_z)*100 << " %).." << std::endl;

        encode_uint_image<FloatNDArray> (imvec, conv_img/*.page(i)*/, true);
        std::string fmt = "tiff";
        //logFile << "try Magick::writeImages " << std::endl;
        std::stringstream ss;
        ss << output_image << /*"_frame_" << i + 1 << */".tif";
        try
        {
            Magick::writeImages(imvec.begin() , imvec.end(), fmt + ":" + ss.str());
        }
        catch (Magick::Warning& w)
        {
            warning ("Magick++ warning: %s", w.what ());
        }
        catch (Magick::ErrorCoder& e)
        {
            warning ("Magick++ coder error: %s", e.what ());
        }
        catch (Magick::Exception& e)
        {
            error ("Magick++ exception: %s", e.what ());
        }
    }

    //im.write(output_image);

}

FloatNDArray gaussFilter(double sigma_xy)
{
    int size_xy = ceil(6*sigma_xy);
    float sum = 0;
    float max = 0;
    float min = 10;

    octave_value_list linspace;
    octave_value result;
    Matrix x, y;

    linspace(0) = -size_xy/2;
    linspace(1) = size_xy/2;
    //linspace(2) = size_xy;
    linspace(2) = 5;

    result = ov_linspace(linspace);
    x = result.matrix_value();
    y = result.matrix_value();

    dim_vector dv (x.cols(), y.cols());
    FloatNDArray kernel(dv);

    for (int i = 0; i < x.cols(); ++i)
    {
        for (int j = 0; j < y.cols(); ++j)
        {

kernel(i, j) = exp(  (-1*pow(x(0,i), 2)/(2*pow(sigma_xy, 2))) +   (-1*pow(y(0,j), 2)/(2*pow(sigma_xy, 2)))    );
           sum += kernel(i,j);
           if (max < kernel(i,j))
            {
                max = kernel(i,j);
            }
        if (min > kernel(i,j))
            {
                min = kernel(i,j);
            }

        }
    }

    return (kernel/sum);
    //return (kernel/max);
    //return (kernel/min);

}

FloatNDArray gaussFilter3D(double sigma_xy, double sigma_z, double z_scale)
{
    int size_xy = ceil(6*sigma_xy);
    int size_z = ceil(6*sigma_z/z_scale);
    float sum = 0;

    octave_value_list linspace;
    octave_value result;
    Matrix x, y, z;

    linspace(0) = -size_xy/2;
    linspace(1) = size_xy/2;
    linspace(2) = size_xy;

    result = ov_linspace(linspace);
    x = result.matrix_value();
    y = result.matrix_value();

    //std::cout << x << std::endl;

    linspace(0) = -size_z/2;
    linspace(1) = size_z/2;
    linspace(2) = size_z;

    result = ov_linspace(linspace);
    z = result.matrix_value()*z_scale;
    //std::cout << z << std::endl;

    dim_vector dv (x.cols(), y.cols(), z.cols());
    FloatNDArray kernel(dv);

    for (int i = 0; i < x.cols(); ++i)
    {
        for (int j = 0; j < y.cols(); ++j)
        {
            for (int k = 0; k < z.cols(); ++k)
            {
kernel(i, j, k) = exp(  (-1*pow(x(0,i), 2)/(2*pow(sigma_xy, 2))) +   (-1*pow(y(0,j), 2)/(2*pow(sigma_xy, 2))) + (-1*pow(z(0,k), 2)/(2*pow(sigma_z, 2)))   );
            sum += kernel(i,j, k);
            }
        }
    }

    return kernel/sum;

    //return kernel;
}

