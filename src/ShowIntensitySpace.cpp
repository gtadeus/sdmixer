
#include "sdmixer.h"


int ShowIntensitiesMain(Matrix& xInt, Matrix& yInt, Matrix Filter[], int nr_of_filters, std::string inputfile, int max_x, int max_y)
{
    octave_value_list result, sub2ind, rot90, result_rot90, horzcat, imwrite, args;
    octave_value ov_result;

    std::string oldfilename = inputfile;

    octave_value_list sort, unique, conv, accumarray, hist_uniq;

    Matrix uniq, vec_i, vec_j, N, histogram, Nhist;

    Matrix size_of_img(1,2);

    //Matrix ind;
   /* Matrix img_rot(max_y, max_x);
    img_rot.fill(0.0);*/

    int NPixel = max_x*max_y;


    size_of_img(0,0) = max_x; //* precision;                       //size of output image
    size_of_img(1,0) = max_y; //* precision;
    dim_vector dvTiff(max_x, max_y);
    uint16NDArray img(dvTiff);

    img.fill(0);
    for (int i = 0; i < nr_of_filters; ++i )
    {
        rot90(0) = octave_value(Filter[i]);
        rot90(1) = -1;
        result_rot90 = feval("rot90", rot90, 2);
        NDArray curr_filter =  result_rot90(0).array_value();
        //std::cout << "curr_filter: rows = " << curr_filter.rows() << " cols = " <<curr_filter.cols() << " " << std::endl;
        try
        {
            //std::cout << curr_filter.max() << std::endl;
            //logFile << curr_filter << std::endl;
            //curr_filter=(100);
            //logFile << octave_value(curr_filter).uint16_array_value() << std::endl;
            //std::cout << img(0) << std::endl;
            //std::cout << curr_filter(0) << std::endl;
            img = octave_value(octave_value(curr_filter).matrix_value()+octave_value(img).matrix_value()).uint16_array_value();
            //std::cout << img(0) << std::endl;
            //img = img + (Filter[i]/(2+5*(i+1)));
        }
        catch ( octave_execution_exception )
        {
            std::cout << "\n\nERROR:Your Filter does not fit with Max xInt/yInt from data!\nResize Filter File!\n\n" << std::endl;
            return 1;
        }
        //img_rot = img_rot + (Filter[i]+((i+1)*200/nr_of_filters));
    }

    args(0) = octave_value(size_of_img);
    args(1) = octave_value(xInt);
    args(2) = octave_value(yInt);
    ov_result = ov_sub2ind(args);
/*
    sub2ind(0) = octave_value(size_of_img);
    sub2ind(1) = octave_value(xInt);
    sub2ind(2) = octave_value(yInt);
    result = feval("sub2ind", sub2ind, 3);
*/
    //std::cout << size_of_img << std::endl;

    Matrix ind = ov_result.matrix_value();

    int bit = 16;
    int mp=pow(2, bit)-1;
    double hist[mp+1];
    double PDFwt[mp+1];
    double cdf[mp+1];

    unique(0) = ov_result.matrix_value();
    result = feval ("unique", unique, 4);
    uniq = result(0).matrix_value();
    vec_i = result(1).matrix_value();
    vec_j = result(2).matrix_value();

    accumarray(0) = octave_value(vec_j);
    accumarray(1) = octave_value(1);
    result = feval ("accumarray", accumarray, 3);
    N = result(0).matrix_value();

    int hist_max = N.column_max().max();

    for (int j = 0; j < uniq.rows(); ++j)
    {
        if(uniq(j) >= 0)
            img(uniq(j)) = ceil(N(j)*((double)(mp)/(double)hist_max));           //linear assignment to 8bit

    }


    for (long g=0; g < (mp+1); ++g)                             //init histogram
        hist[g]=0;

    for (int i = 0; i < NPixel; ++i)                            //fill histogram
        hist[(int)img(i)]+=1;


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

    double Pu = 1 * Pmax;                 //thresholding the PDF
    for (int g=0; g < mp; ++g)              //calculate weightet PDF
    {
        PDFwt[g]=pow( (PDFwt[g]-Pl)/(Pu-Pl), 0.8) * Pu;
        PDFwt[g]*=NPixel;                   //transform to histogram range [0,NPixel]
    }

    cdf[0]=sqrt(PDFwt[0]);
    for (long b=1; b<(mp+1); ++b)   //from PDF to CDF with square root of hist values
        cdf[b]=cdf[b-1]+sqrt(PDFwt[b]);


    for (int i = 0; i < NPixel; ++i)         //for each pixel remap
        img(i)=ceil(  mp * (double) (cdf[(int)img(i)]-cdf[0]) / (cdf[mp]-cdf[0]));



   // std::cout << "img: rows = " << img.rows() << " cols = " <<img.cols() << " " << std::endl;
   // std::cout << "Filter: rows = " << Filter[0].rows() << " cols = " <<Filter[0].cols() << " " << std::endl;

  /*  imwrite(0) = octave_value(img);
    imwrite(1) = inputfile.append("_IntSpace.tif");
    result = feval("imwrite", imwrite, 2);
*/


    //std::cout << "rotiert" << std::endl;

    rot90(0) = octave_value(img);
    rot90(1) = 1;
    result_rot90 = feval("rot90", rot90, 2);
    NDArray img_out = result_rot90(0).array_value();

    for (int i = 0; i < NPixel; ++i)         //back to 8bit
    {
        img_out(i) = (img_out(i)*(1.0/(mp)));
    }

    //Magick::InitializeMagick(NULL);
    std::vector<Magick::Image> imvec;
    logFile << "start encode_uint_image " << std::endl;
    encode_uint_image<NDArray> (imvec, img_out, true);
    std::string fmt = "tif";
    logFile << "try Magick::writeImages " << std::endl;
    try
    {
        Magick::Image Mimg = imvec[0];
        Mimg.write(fmt + ":" + inputfile.append("_IntSpace.tif"));
        //Magick::writeImages(imvec.begin() , imvec.end(), fmt + ":" + inputfile.append("_IntSpace.tif"));
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

  /*  imwrite(0) = octave_value(img_rot);
    imwrite(1) = inputfile.append("_IntSpace.tif");
    result = feval("imwrite", imwrite, 2);
*/
    horzcat(0) = octave_value(xInt);
    horzcat(1) = octave_value(yInt);
    result = feval("horzcat", horzcat, 2);

    std::ofstream ofsFile;
    std::string xy_outputFile = oldfilename.append("_IntSpace.out");

    ofsFile.open(xy_outputFile.c_str());

    if(!std::ifstream(xy_outputFile.c_str()))
    {
        std::cout << "\nERROR: ShowIntensities could not write in " << xy_outputFile.c_str() << std::endl;
    }

    else
    {
        ofsFile << result(0).matrix_value();
        ofsFile.close();
    }

    return 0;

}

