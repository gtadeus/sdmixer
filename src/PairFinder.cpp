//#include "sdmixer.h"

#include "sdmixer.h"
//input must be sorted by frames!!
int FindPairs(Matrix& input, Matrix& output, int dimensions, double Offset[], double Epsilon[], int FrameColumnC, int& NrOfDifferentFrames, int& multiple_counter, bool create_output = false, int row_stop = -1)
{
        int pair=0;
        int increment=1;
        NrOfDifferentFrames=0;
        multiple_counter=0;
        int numpairs = 0;
        double EllipsoidSum = 0;

        int last_frame = -1;

        const int raw_data_cols = input.cols();
        const int raw_data_rows = input.rows();
        int endrow = raw_data_rows;

        if(row_stop != -1)
            endrow = row_stop;


        double (*pFullData) = new double[raw_data_cols*raw_data_rows];
        double *fv_pFullData = input.fortran_vec();

        memcpy(pFullData, fv_pFullData, raw_data_cols*raw_data_rows*sizeof(double));

        std::vector< std::vector<double> > vecPairsOut;

        std::vector<int> grouped_rows;

        while (pair < endrow)
        {
            double EllipsoidSumR=0;
            double EllipsoidSumL=0;
            if (last_frame != pFullData[pair+FrameColumnC*raw_data_rows])
            {
                last_frame = pFullData[pair+FrameColumnC*raw_data_rows];
                NrOfDifferentFrames++;
            }

            if( (pFullData[pair+FrameColumnC*raw_data_rows] == pFullData[pair+increment+FrameColumnC*raw_data_rows]) )
            {
                    for (int d = 0; d < dimensions; ++d)
                    {
                        double tempL = ((pFullData[pair+d*raw_data_rows] - Offset[d]) - pFullData[pair+increment+d*raw_data_rows]);
                        double tempR = ((pFullData[pair+increment+d*raw_data_rows] - Offset[d]) - pFullData[pair+d*raw_data_rows]);

                        tempL*=tempL;
                        tempR*=tempR;
                        tempL/=(Epsilon[d]*Epsilon[d]);
                        tempR/=(Epsilon[d]*Epsilon[d]);

                        EllipsoidSumL += tempL;
                        EllipsoidSumR += tempR;

                    }
                if (EllipsoidSumL <= 1)
                {
                    std::vector<double> temp;

                    for(int i = 0; i < raw_data_cols; ++i)
                        temp.push_back(pFullData[pair+i*raw_data_rows]);

                    for(int i = 0; i < raw_data_cols; ++i)
                        temp.push_back(pFullData[pair+increment+i*raw_data_rows]);

                    vecPairsOut.push_back(temp);

                    numpairs++;

                    grouped_rows.push_back(pair);
                    grouped_rows.push_back(pair+increment);
                }
                else if (EllipsoidSumR <= 1)
                {
                    std::vector<double> temp;

                    for(int i = 0; i < raw_data_cols; ++i)
                        temp.push_back(pFullData[pair+increment+i*raw_data_rows]);

                    for(int i = 0; i < raw_data_cols; ++i)
                        temp.push_back(pFullData[pair+i*raw_data_rows]);

                    vecPairsOut.push_back(temp);
                    numpairs++;

                    grouped_rows.push_back(pair);
                    grouped_rows.push_back(pair+increment);
                }

                if (increment != raw_data_rows)
                    increment++;

                 EllipsoidSumR=0;

            }
            else
            {
                pair++;
                increment = 1;
                EllipsoidSumR=0;
            }

        }

    std::sort ( grouped_rows.begin(), grouped_rows.end() );
    multiple_counter = std::abs(std::distance(std::unique ( grouped_rows.begin(), grouped_rows.end()), grouped_rows.end()));

    if ( create_output )
    {
    output.resize(numpairs, raw_data_cols*2);

    for(int i = 0; i < vecPairsOut.size(); ++i)
         for(int j=0; j< vecPairsOut[i].size(); ++j)
               output(i,j)= vecPairsOut[i][j] ;

    }


    delete[] pFullData;

    return numpairs;
}

int PairFinderMain(boost::property_tree::ptree *pt, fs::path *config_file_full_path, std::string inputfile)
{

    // we start with the definitions from the config file
    std::string str_xOffset;
    std::string str_yOffset;
    std::string str_zOffset;
    std::string str_NM_PER_PX;
    std::string str_eps_x;
    std::string str_eps_y;
    std::string str_eps_z;
    std::string str_fishing;
    std::string str_fishing_end;
    std::string str_fishing_inc;
    std::string str_fishing_subset;

    std::string str_grouping;
    std::string str_x_length;
    std::string str_y_length;
    std::string str_z_length;
    std::string str_histogram_out;
    std::string str_groupingFileStamp;
    std::string str_PairFinderFileStamp;

    //new since v1.9 CameraOrientation, LongShort reconstruction
    std::string str_CameraOrientation;


    try{
        str_xOffset = pt->get<std::string>("PairFinder.xOffset");
        str_yOffset = pt->get<std::string>("PairFinder.yOffset");
        str_zOffset = pt->get<std::string>("PairFinder.zOffset");
        str_NM_PER_PX = pt->get<std::string>("General.NM_PER_PX");
        str_eps_x = pt->get<std::string>("PairFinder.eps_x");
        str_eps_y = pt->get<std::string>("PairFinder.eps_y");
        str_eps_z = pt->get<std::string>("PairFinder.eps_z");
        str_PairFinderFileStamp = pt->get<std::string>("PairFinder.FileStamp");

        // new since 1.3 fishing...
        str_fishing = pt->get<std::string>("xy-fishing.fishing_on");
        str_fishing_inc = pt->get<std::string>("xy-fishing.fishing_inc");
        str_fishing_end = pt->get<std::string>("xy-fishing.fishing_end");
        str_fishing_subset = pt->get<std::string>("xy-fishing.fishing_subset");

        // ver1.5 grouping
        str_grouping = pt->get<std::string>("grouping.grouping_on");
        str_groupingFileStamp = pt->get<std::string>("grouping.OutputFileStamp");
        str_x_length = pt->get<std::string>("grouping.x_length");
        str_y_length = pt->get<std::string>("grouping.y_length");
        str_z_length = pt->get<std::string>("grouping.z_length");
        str_histogram_out = pt->get<std::string>("grouping.histogram_out");
    }
    catch( const std::exception & ex ){

        std::cout << "\nERROR: PairFinder found syntax error in config file: "  << std::endl << config_file_full_path->string() << std::endl << ex.what() << std::endl;
        return 1;

    }
    //new since v1.9 CameraOrientation, LongShort reconstruction
    /*try{
        str_CameraOrientation = pt->get<std::string>("Camera.orientation");

    }
    catch( const std::exception & ex ){

        logFile << "\nWARNING: no or wrong Camera Orientation Information found in config file, using left-right (LR)"  << std::endl;
        std::cout << "\nWARNING: no or wrong Camera Orientation Information found in config file, using left-right (LR)"  << std::endl;
        //return 1;

    }*/


//*****************************
    int PairFinderStatus = 0;
    int NM_PER_PX = atoi(str_NM_PER_PX.c_str());
    double xOffset = atof(str_xOffset.c_str()) * NM_PER_PX;
    double yOffset = atof(str_yOffset.c_str()) * NM_PER_PX;
    double zOffset = atof(str_zOffset.c_str());

    Matrix epsilon = Matrix(1,3);
    epsilon(0,0) = atof(str_eps_x.c_str()) * NM_PER_PX;
    epsilon(0,1) = atof(str_eps_y.c_str()) * NM_PER_PX;
    epsilon(0,2) = atof(str_eps_z.c_str()) * NM_PER_PX;
    double Epsilon[3];
    Epsilon[0] = atof(str_eps_x.c_str()) * NM_PER_PX;
    Epsilon[1] = atof(str_eps_y.c_str()) * NM_PER_PX;
    Epsilon[2] = atof(str_eps_z.c_str()) * NM_PER_PX;
    double Offset[3];
    double tempOffset[3];
    Offset[0]=xOffset;
    Offset[1]=yOffset;
    Offset[2]=zOffset;

    // new since v1.3 fishing...
    int fishing_on = atoi(str_fishing.c_str());
    bool fishing = false;
    if (fishing_on != 0)
    {
        fishing = true;
    }
    int fishing_end = atoi(str_fishing_end.c_str());
    double fishing_inc = atof(str_fishing_inc.c_str());
    int fishing_subset = atoi(str_fishing_subset.c_str());

    Matrix fishing_results, fishing_FoundPairs, fishing_OffsetX, fishing_OffsetY;

    int size_fishingMatrix = 2*fishing_end+1;

    if(fishing)
    {
        fishing_results.resize(size_fishingMatrix, size_fishingMatrix, 0);
        fishing_FoundPairs.resize(size_fishingMatrix, size_fishingMatrix, 0);
        fishing_OffsetX.resize(1, size_fishingMatrix, 0);
        fishing_OffsetY.resize(1, size_fishingMatrix, 0);

        for (int i = 0; i < size_fishingMatrix; ++i)
        {
            fishing_OffsetX(0,i) = Offset[0] + ( -1 * fishing_end*fishing_inc + fishing_inc*i) * NM_PER_PX;
            fishing_OffsetY(0,i) = Offset[1] + ( -1 * fishing_end*fishing_inc + fishing_inc*i) * NM_PER_PX;
        }

    }

    // since v1.5 grouping
    int grouping_on = atoi(str_grouping.c_str());
    bool grouping = false;

    if(grouping_on != 0)
    {
        grouping = true;
    }

    double grouping_matrix[3];
    grouping_matrix[0] = atof(str_x_length.c_str());
    grouping_matrix[1] = atof(str_y_length.c_str());
    grouping_matrix[2] = atof(str_z_length.c_str());

    int histogram_out = atoi(str_histogram_out.c_str());
    bool grouping_histogram_out = false;

    if(histogram_out != 0)
    {
        grouping_histogram_out = true;
    }


    int xColumnOctave = 1;                      // in octave array (matrices) start from 1
    int xColumnC = 0;                           // c++ style arrays start from 0. be aware, it's a trap



    std::string outputFile=inputfile;
    outputFile.append(str_PairFinderFileStamp.c_str());

//###################################################

    Matrix rsRawDataByFrame;

    int FrameColumnC, totalColNr, FrameColumnOctave;
    unsigned long totalRowNr ;


    logFile << "loading file " << inputfile << std::endl << std::endl;

    load_file(inputfile, rsRawDataByFrame);

    std::ifstream ifs(inputfile.c_str());
    std::string firstLine;
    getline (ifs, firstLine);                           // load first line (header) and

    const int dimensions = ReturnDimensions("Position", firstLine);       //count the occurence of "Position" --> this is the nr of dimensions!
    logFile << "Dimensions:\n " << dimensions << std::endl;
    logFile << "Scale (nm per px):\n " << NM_PER_PX << std::endl;
    logFile << "Ellipsoid Dimensions (nm):\n" << epsilon << std::endl;
    logFile << "Offset Matrix (nm):\n " << xOffset << " " << yOffset << " " << zOffset << std::endl;


    FrameColumnC = dimensions;       // C++ Frame Col !             here we have to insert the actual dimensions of the data set
    FrameColumnOctave = dimensions + 1;        // C++ Frame Col !

    sortrows(rsRawDataByFrame, FrameColumnC);

    totalRowNr = rsRawDataByFrame.rows();
    totalColNr = rsRawDataByFrame.cols();

    logFile << "\n" << totalRowNr << " rows; ";

    int NrOfDifferentFrames = 0;
    int multiple_counter = 0;
    std::vector<double> vecPairsOut;
    Matrix PairsOut;

    std::ofstream ofsFile;

//goto finishpairfinder;

    if(std::ifstream(outputFile.c_str()))
    {
        remove(outputFile.c_str());                 // if there was a previos output file, it will be deleted
                                                    // no error handling... has to be added.. :(
    }

    ofsFile.open(outputFile.c_str(), std::ios::app);         // we write the output file now

    int numpairs;

        //this is redundant, but still... we do the fishing in a seperate loop.
        // we need a function for the PF main loop.. :(
        // new since v1.3 fishing:
        if(fishing)
        {
            std::vector<int> frames;
            for(int i = 0; i < totalRowNr; ++i)
                frames.push_back(rsRawDataByFrame(i, FrameColumnC));

            int endrow = std::abs(std::distance(std::find(frames.begin(), frames.end(), fishing_subset+1), frames.begin()));

            if (endrow >= totalRowNr)
            {
                std::cout << "WARNING: chosen fishing subset > NrOfFrames, fishing on whole dataset" << std::endl;
                logFile << "WARNING: chosen fishing subset > NrOfFrames, fishing on whole dataset" << std::endl;
            }
            for (int ii = 0; ii < size_fishingMatrix; ++ii)
            {
                for (int jj = 0; jj < size_fishingMatrix; ++jj)
                {
                    tempOffset[0] = fishing_OffsetX(0,ii);
                    tempOffset[1] = fishing_OffsetY(0,jj);
                    tempOffset[2] = Offset[2];                          // no fishing in z

                    fishing_FoundPairs(ii,jj)   =  FindPairs(rsRawDataByFrame, PairsOut, dimensions, tempOffset, Epsilon, FrameColumnC, NrOfDifferentFrames, multiple_counter, false, endrow);

                }
            }

            // here we choose our best offset...
            int max = fishing_FoundPairs(0,0) ;       // start with max = first element
            int max_xOffset = Offset[0];
            int max_yOffset = Offset[1];
            for (int ii = 0; ii < size_fishingMatrix; ++ii)
            {
                for (int jj = 0; jj < size_fishingMatrix; ++jj)
                {
                    if(fishing_FoundPairs(ii,jj) > max)                 // Matrix(row, col)
                    {
                        max = fishing_FoundPairs(ii,jj);
                        max_xOffset = Offset[0] + ( -1 * fishing_end*fishing_inc + fishing_inc*ii) * NM_PER_PX;
                        max_yOffset = Offset[1] + ( -1 * fishing_end*fishing_inc + fishing_inc*jj) * NM_PER_PX;
                    }
                }
            }

            //set the Offset!
            Offset[0] = max_xOffset;
            Offset[1] = max_yOffset;
            logFile << "\n**fishing:" << std::endl;
            logFile << "inc = " << fishing_inc << ";  end = " << fishing_end << ";  subset = " << fishing_subset << std::endl;
            logFile << "results:" << std::endl;
            logFile << fishing_FoundPairs << std::endl;
            logFile << "max = " << max << std::endl;
            logFile << "New Offset Matrix (nm):\n " << Offset[0] << " " << Offset[1] << " " << Offset[2] << "\n" << std::endl;
            logFile << "New Offset Matrix (px):\n " << Offset[0]/NM_PER_PX << " " << Offset[1]/NM_PER_PX << " " << Offset[2]/NM_PER_PX << "\n\n" << std::endl;

        }


    //PF v2*********************************

    double startTime = omp_get_wtime();

    numpairs = FindPairs(rsRawDataByFrame, PairsOut, dimensions, Offset, Epsilon, FrameColumnC, NrOfDifferentFrames, multiple_counter, true);

    double endTime = omp_get_wtime();


    logFile  << "PairFinder finished in " << (endTime - startTime) << " seconds" << std::endl;
    logFile << NrOfDifferentFrames << " different frames" << std::endl;
    if (numpairs != 0)
    {
        ofsFile << PairsOut << std::endl;
        ofsFile.close();
        logFile << numpairs << " pairs" << std::endl;
        logFile << multiple_counter << " multiple hits (" << (double) multiple_counter / numpairs * 100 << " %)"  << std::endl;
    }
    else
        std::cout << "WARNING: found " << numpairs << " pairs. maybe wrong xOffset? (current: " << Offset[0]/NM_PER_PX << " px)" << std::endl;


//finishpairfinder:

    // gesamte pairs_out ist in PairsOut, nach frames sortiert geladen.



    if(grouping)
    {
            bool GroupingMatrixIncomplete = false;
    for (int i = 0; i < dimensions; i++)
    {
        if(grouping_matrix[i] < 0.1)
        {
            GroupingMatrixIncomplete = true;
        }
    }
    if(GroupingMatrixIncomplete)
    {
        PairFinderStatus = 2;
        std::cout << "\nWARNING: Grouping Matrix contains a very small value,\ngrouping may be futile with this " << dimensions << "D data set."<< std::endl;
        std::cout << "Grouping Matrix (nm):" <<std::endl;
        for (int i = 0; i < dimensions; ++i)
        {
            std::cout << " " << grouping_matrix[i];
        }
        std::cout <<"\n"<< std::endl;
    }

        std::cout << "grouping... " << std::endl;
        std::ofstream groupFile, histogramFile, wtFile;
        std::string groupHistogram = outputFile+"_grouphist.out";
        std::string weightedCoordinates = outputFile+ str_groupingFileStamp;
        int group_counter = 0;

        int xColIntC = (PairsOut.cols()/2) + dimensions + 1;
        int yColIntC = dimensions + 1;
        int PairsOutTotalCols = PairsOut.cols();

        const int cols = PairsOut.cols();
        const int rows = PairsOut.rows();

        double (*pData) = new double[rows*cols];
        double *pD1 = PairsOut.fortran_vec();

        memcpy(pData, pD1, rows*cols*sizeof(double));
        double ctime1 = omp_get_wtime();


        std::vector< std::vector<int> > local_nodes;
        std::vector< std::vector<double> > centers;
        int j = 0;
        std::vector<int> active_old;
        std::vector<int> active_new;
        bool found = false;

        for (int i = 0; i < rows; ++i)
        {
            if ( i == 0 || pData[i+FrameColumnC*rows] != pData[(i-1)+FrameColumnC*rows] )
            {
                active_old = active_new;
                active_new.clear();
            }
            found = false;
            for_each(active_old.begin(), active_old.end(), [&](int k)
            {
              if(found == false)
                {
                    double nodes_weight;
                    double dx[dimensions];
                    int EllipsoidSum=0;

                    for (int d = 0; d < dimensions; ++d)
                        dx[d] = centers[k][d]-pData[i+d*rows];

                    for (int d = 0; d < dimensions; ++d)
                        EllipsoidSum += pow(dx[d]/grouping_matrix[d] , 2);

                   if (EllipsoidSum < 1 && pData[local_nodes[k][local_nodes[k].size()-1]+FrameColumnC*rows] - pData[i+FrameColumnC*rows]<= 1)
                    {
                        double temp[dimensions];
                        nodes_weight = 0;
/*
                        for(int nodeIdx : local_nodes[k])
                            nodes_weight += PairsOut(nodeIdx, yColIntC);*/
                        for_each(local_nodes[k].begin(), local_nodes[k].end(), [&](int nodeIdx)
                        {
                            nodes_weight += pData[nodeIdx+yColIntC*rows];
                        });
                        for (int d = 0; d < dimensions; ++d)
                        {
                            temp[d]=nodes_weight*centers[k][d];
                            temp[d]+= (pData[i+yColIntC*rows] * pData[i+d*rows]);
                            temp[d]/=(nodes_weight + pData[i+yColIntC*rows]);

                            centers[k][d] = temp[d];//joined_centers[jc_counter+d];
                        }

                        local_nodes[k].push_back(i);
                        active_new.push_back(k);
                        found = true;

                    }

                }
            });


           if(found == false)
            {
                std::vector<int> tmp;
                tmp.push_back(i);
                local_nodes.push_back(tmp);
                active_new.push_back(j);

                std::vector<double> tmp2;
                for (int d = 0; d < cols; ++d)
                {
                    tmp2.push_back(pData[i+d*rows]);
                }
                centers.push_back(tmp2);
                j++;
           }
        }

        double ctime2 = omp_get_wtime();

        std::cout << "finished grouping in "<< ctime2-ctime1 << " seconds\n" << std::endl;
        logFile << "finished grouping in "<< ctime2-ctime1 << " seconds\n" << std::endl;

            std::vector<double> vecIntLeft;
            std::vector<double> vecIntRight;
            std::vector<double> groupSize;
            std::vector<double> realGroupSize;


            for ( std::vector<  std::vector<int> >::size_type u = 0; u < local_nodes.size(); u++)
            {
                double Navg=0;
                double Nsum_second = 0;
                double tmp_groupSize=0;
                for ( std::vector<int>::size_type v = 0; v < local_nodes[u].size(); v++)
                {
                    Navg+=PairsOut(local_nodes[u][v], yColIntC);
                    Nsum_second += PairsOut(local_nodes[u][v], xColIntC);
                    tmp_groupSize+=1;
                }
                vecIntLeft.push_back(Navg);
                vecIntRight.push_back(Nsum_second);
                groupSize.push_back(tmp_groupSize);
                if(tmp_groupSize>1)
                    realGroupSize.push_back(tmp_groupSize);

            }

            Matrix IntLeft, IntRight, frequency, frequency2;
            IntLeft = vector2matrix(vecIntLeft);
            IntRight = vector2matrix(vecIntRight);
            frequency = vector2matrix(groupSize);
            frequency2 = vector2matrix(realGroupSize);


            Matrix GroupFile(dim_vector (IntRight.rows(), 0));

            Matrix glob_centers(dim_vector (IntRight.rows(), cols ));
            for(int i = 0; i < centers.size(); ++i)
            {
                for(j=0; j< centers[i].size(); ++j)
                    glob_centers(i,j)= centers[i][j] ;

            }

            GroupFile=GroupFile.append(glob_centers.extract(0, 0, glob_centers.rows()-1, yColIntC-1));
            GroupFile=GroupFile.append(IntLeft);
            GroupFile=GroupFile.append(glob_centers.extract(0, yColIntC+1, glob_centers.rows()-1, xColIntC-1));
            GroupFile=GroupFile.append(IntRight);
            GroupFile=GroupFile.append(glob_centers.extract(0, xColIntC+1, glob_centers.rows()-1, glob_centers.cols()-1));

        std::string groupFileName = outputFile+"_group.out";
        groupFile.open(groupFileName.c_str(), std::ios::out);
        groupFile << GroupFile << std::endl;
        groupFile.close();
       // std::cout << "grouping ready" << std::endl;
        Matrix GroupStat;
/*
        for(int i = 0; i < local_nodes.size(); ++i)
        {


            for(j=0; j< local_nodes[i].size(); ++j)
            std::cout << local_nodes[i][j] << "  ";// << std::endl;

            std::cout << std::endl;
        }
        std::cout << std::endl;*/

        logFile << "\nGrouping stats:" <<std::endl;
        logFile << "Grouping Matrix (nm):" <<std::endl;
        for (int i = 0; i < dimensions; ++i)
        {
            logFile << " " << grouping_matrix[i];
        }
        logFile <<std::endl;

        int grouped_pairs=frequency2.rows();
        if(grouped_pairs == 0)
        {
                PairFinderStatus = 2;

                std::cout << "\nWARNING: found "<< grouped_pairs << " pairs which were \"on\" for more than one frame.\nUsing "
                << str_PairFinderFileStamp.c_str() << ", for further analysis\n" << std::endl;
                if (dimensions == 3 && grouping_matrix[2] < 0.1)
                {
                    std::cout << "it seems as if you forgot to provide a z_length value\nin the [grouping] section of the config file for this 3D data set! " << std::endl;
                }else{
                    std::cout << "maybe Grouping Matrix too small?" << std::endl;
                }
                std::cout << "Grouping Matrix (nm):" <<std::endl;
                for (int i = 0; i < dimensions; ++i)
                {
                    std::cout << " " << grouping_matrix[i];
                }
                std::cout <<"\n"<< std::endl;
                logFile << "\nWARNING: found "<< grouped_pairs <<" pairs which were \"on\" for more than one frame.\nUsing "
                << str_PairFinderFileStamp.c_str() << ", for further analysis\n" << std::endl;
        }
        else{
                wtFile.open(weightedCoordinates.c_str(), std::ios::out);
                wtFile << GroupFile;
                wtFile.close();

                std::string stat_out_string;
                if (grouped_pairs< 2)
                {
                    std::cout << "\nWARNING: found only " << grouped_pairs << " group." << std::endl;
                    logFile << "\nWARNING: found only " << grouped_pairs << " group." << std::endl;

                }else{

                    statistics1D(frequency2, GroupStat, stat_out_string);


                    logFile << grouped_pairs << " pairs (" << (double) grouped_pairs / numpairs * 100 << " %)"  << " were \"on\" for more than one frame." << std::endl;
                    logFile << "statistical analysis excluding pairs which were \"on\" for one frame:\n" << std::endl;
                    logFile << stat_out_string << std::endl;
                }


                statistics1D(frequency, GroupStat, stat_out_string);

                logFile << "statistical analysis including all pairs:\n" << std::endl;
                logFile << stat_out_string << std::endl;

                if(grouping_histogram_out)
                {
                    Matrix histogram;
                    octave_value_list accumarray;

                    accumarray(0) = octave_value(frequency);
                    accumarray(1) = octave_value(1);
                    octave_value_list result = feval ("accumarray", accumarray, 3);
                    histogram = result(0).matrix_value();
                    logFile << "writing histogram into " << groupHistogram.c_str() << "\n" <<std::endl;


                    Matrix groupsize=Matrix(histogram.rows(),1);
                    for(int i = 0; i < histogram.rows(); ++i)
                    {
                        groupsize(i, 0) = i+1;
                    }
                    groupsize = groupsize.append(histogram);
                    histogramFile.open(groupHistogram.c_str(), std::ios::out);
                    histogramFile << "#frames_on  nr_of_pairs" << std::endl;
                    histogramFile << groupsize << std::endl;
                    histogramFile.close();
                }
        }


    }

    return PairFinderStatus;
}

Matrix vector2matrix(std::vector<double>& vec)
{
    double* a = &vec[0];
    Matrix m(dim_vector (vec.size(), 1));
    double *pm = m.fortran_vec();
    memcpy(pm, a, vec.size()*sizeof(double));

    return m;
}

void deleteRow(Matrix& m, const int& row)
{
    idx_vector iv(row);
    m.delete_elements(0, iv);
}
void deleteAllRows(Matrix& m)
{
    idx_vector iv(0, m.rows(), 1);
    m.delete_elements(0, iv);
}

int CalculateFrequency(Matrix& raw_data, Matrix& Frequency)
{
    octave_value_list unique, result, accumarray;
    Matrix uniq, vec_i, vec_j;

    unique(0) = octave_value(raw_data);
    result = feval ("unique", unique, 4);
    uniq = result(0).matrix_value();
    vec_i = result(1).matrix_value();
    vec_j = result(2).matrix_value();

    accumarray(0) = octave_value(vec_j);
    accumarray(1) = octave_value(1);
    result = feval ("accumarray", accumarray, 3);
    Frequency = result(0).matrix_value();

    return 0;
}
int CentroidPoints(Matrix& coordinates, Matrix& intensities, Matrix& centroid)
{

    double sum = 0;
    double sum_int = 0;
    for (int i = 0; i < intensities.rows(); i++)
    {
        sum_int+=intensities(i, 0);
    }

    for(int i = 0; i < coordinates.cols(); ++i)
    {
        sum = 0;
        for (int j = 0; j < coordinates.rows(); ++j)
        {
            sum += coordinates(j, i)*(intensities(j, 0)/sum_int);
        }
        centroid(0, i) = sum;

    }

    //centroid /= coordinates.rows();

    return 0;                                               //return
}

int statistics1D(Matrix& row_vec, Matrix& statisticsOut, std::string& formatedOutput)
{
    octave_value_list statistics, result;

    statistics(0) = octave_value(row_vec);

    result = feval("statistics", statistics);

    statisticsOut = result(0).matrix_value();

    std::stringstream ssOutput;
    ssOutput
            << "min\t=\t" << statisticsOut(0, 0) << "\n"
            << "median\t=\t" << statisticsOut(2, 0) << "\n"
            << "max\t=\t" << statisticsOut(4, 0) << "\n"
            << "mean\t=\t" << statisticsOut(5, 0) << "\n"
            << "stdev\t=\t" << statisticsOut(6, 0) << "\n"
            << "skewness=\t" << statisticsOut(7, 0) <<  std::endl;

    formatedOutput = ssOutput.str();

    //min
    //Q.25
    //Q.50
    //Q.75
    //max
    //mean
    //stdv
    //skewness
    //kurtosis

    return 0;
}

int appendMatrix(Matrix& m, Matrix& subMatrix)
{
    octave_value_list result, vertcat;

    vertcat(0) = octave_value(m);
    vertcat(1) = octave_value(subMatrix);

    result = feval ("vertcat", vertcat);                // concatenate vertical

    m = result(0).matrix_value();

    return 0;
}
