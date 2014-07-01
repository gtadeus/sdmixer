
#include "sdmixer.h"


#define NAME "sdmixer"
#define SDVERSION "1.95"


//using namespace std;

std::string PairFinderFileStamp;        // bad place for definition here..
std::string PairFinderFileStampGroupingFailed;
std::string FilterFileStamp;            // one should not use global defs in multi file projects
                                        // we don't know when they are set, but it seems to work in this case..
std::ofstream logFile;
//LogFile* LogFile::m_pInstance = NULL;


int main(int argcc, char* argvv[])
{
    Magick::InitializeMagick(NULL);
                                        // welcome to the main  :)
    std::stringstream ss_welcome;

    ss_welcome << "*** " << NAME <<" v. " << SDVERSION  << " - Analysis of 2D/3D multicolor SD-dSTORM data ***" ;



    std::cout << "\n" << ss_welcome.str() << std::endl;

    int dir_match_count = 0;
    int file_match_count = 0;

    std::string default_config_filename = NAME".conf";      // if u haven't specified it in the command-line
    std::string logfilename = get_date().append("_" NAME ".log");

    //std::string strDimensions;

    unsigned long file_count = 0;
    unsigned long dir_count = 0;
    unsigned long err_count = 0;

    std::stringstream ssusage_string;
    ssusage_string  << "\nUsage: " << argvv[0] <<
    " PATH [OPTIONS] \nruns "<< NAME << " in working directory \n\n" <<
    "optional arguments:\n" <<
    "\t-d \t\t : if specified, 2D image output of a 3D data set\n" <<
    "\t-o \t\t : one-color flag, overrides filter information\n" <<
    "\t-c FILENAME \t : specify an alternative configuration file\n" <<
    "\t-s STAGE_VALUE \t : skip stage, value has to be 1, 2 or 3.\n\t\t\t   1: skip PairFinder" <<
    "\n\t\t\t   2: skip Filter" <<
    "\n\t\t\t   3: skip Reconstructor\n" <<
    "\t-e STAGE_VALUE \t : end after stage, value has to be 1 or 2.\n\t\t\t   1: end after PairFinder" <<
    "\n\t\t\t   2: end after Filter\n" <<
    "\t-h \t\t : print this usage information\n"
    << std::endl;

    int nonopt_count = 0;                          // counter for optional argumnets
    bool force_dim = false;                     // -d
    int hflag = 0;                              // -h
    bool one_color = false;
    char *opt_conf_file = NULL;

    int skip_stage = 0;                         // value of -s
    int end_stage = 0;                          // value of -e

    int index;
    int c;
    opterr = 0;

    while ((c = getopt (argcc, argvv, "dc:s:e:oh")) != -1)
    {
        switch (c)
        {
            case 'd':
                force_dim = true;                   // force 2d analysis of 3d data set (why should one do this :)?
                break;
            case 'c':
                opt_conf_file = optarg;             // c-style string to our optional config file
                break;
            case 's':
                skip_stage = atoi(optarg);
                if (skip_stage != 1 && skip_stage != 2 && skip_stage != 3)
                {
                    std::cout << "\nERROR: skip stage: value has to be 1, 2 or 3\n" << std::endl;
                    return 1;
                }
                break;
            case 'e':
                end_stage = atoi(optarg);
                if (end_stage != 1 && end_stage != 2)
                {
                    std::cout << "\nERROR: end after stage: value has to be 1 or 2\n" << std::endl;
                    return 1;
                }
                break;
            case 'o':
                one_color = true;                   // one-color flag, overrides filter information
                break;
            case 'h':
                hflag = 1;
                std::cout << ssusage_string.str() ;
                return 1;
                break;
            case '?':
                if (optopt == 'c' || optopt == 's')
                {
                    std::cout << "\nERROR: Option '-" << (char) optopt << "' requires an argument\nrun " << argvv[0] << " -h for help\n";
                }

                else if (isprint (optopt))
                {
                    std::cout << "\nERROR: Unknown option '-" << (char) optopt << "'\nrun " << argvv[0] << " -h for help\n";
                }
                else                         // who the f*** would enter "umlaute" und/oder ß as an argument???
                {
                    std::cout << "\nERROR: Unknown optional character '-" << (wchar_t)optopt <<
                    "', Non-ASCII characters are not valid arguments!\nrun " << argvv[0] << " -h for help\n";
                }
                return 1;
            default:
            abort ();
        }
    }

    if (optind == argcc && hflag == 0)                      // jeah, genius give me an working dir or an help flag!
    {
        std::cout << "\nERROR: no path to working dir specified!\nrun " << argvv[0] << " -h for help\n";
        return 1;
    }

    for (index = optind; index < argcc; index++)
    {
        ++nonopt_count;
    }

    if (nonopt_count > 1)                               // one working dir is enough!
    {
        std::cout << "\nERROR: only one path to working dir is allowed!\nrun " << argvv[0] << " -h for help\n";
        return 1;
    }


    std::vector<std::string> dir_list;
    std::vector<std::string> file_list;


    boost::property_tree::ptree ptree_config_file;
    fs::path full_path( fs::initial_path<fs::path>() );
    full_path = fs::system_complete( fs::path( argvv[optind] ) );           // full path is the specified working dir

                                                    // this is seriously strange work-around to run octave main
                                                                        // but it works.
                                                                        // it would cause a segmentation fault when octave_main()
                                                                       // would be run with not the same nr. of args that were passed to the
                                                                        //sdmixer main.. ???
    string_vector om(argcc);                                            //
    for (int a = 0; a<argcc; ++a)
    {
         om(a) = "--silent";
         //om(a) = "-x";
    }
    octave_main(argcc, om.c_str_vec(), 1);                              // now its ok
                                                                        // if you have a better solution write to: g.tadeus@fu-berlin.de



    if ( !fs::exists( full_path ) )
    {

        std::cout << "\nERROR: working directory not found: " << full_path.string() << std::endl << std::endl;
        return 1;
    }

    fs::path config_file_full_path = operator/(full_path, default_config_filename);     // add name of default config file to working dir
                                                                                    // really nice construct this "operator/" it adds
                                                                                    // the right directory_seperator for unix "/"
                                                                                    // let's see if it works with windows...


    if (opt_conf_file != NULL)                                            // user specified an other location for the config file
    {
        config_file_full_path = opt_conf_file;                             // redefinition of the path above

        if (!std::ifstream(config_file_full_path.string().c_str()))
        {
            std::cout << "\nERROR: failed opening specified configuration file " << config_file_full_path.string()  << std::endl;
            return 1;
        }

    }
    else
    {
        if (!std::ifstream(config_file_full_path.string().c_str()))
        {
            std::cout << "\nERROR: no configuration file found in the working directory." << std::endl;
            return 1;
        }

    }


    fs::path log_file_full_path = operator/(full_path,logfilename);
   // fs::path log_file_full_path2 = operator/(full_path,"log.out");

    logFile.open(log_file_full_path.string().c_str());

    //LogFile::GetInstance()->OpenLogFile(log_file_full_path2.string().c_str());
/*
    LogFile &lg = *LogFile::GetInstance();
    lg.OpenLogFile(log_file_full_path2.string().c_str());

    lg << "Test\n";

*/


    if (!std::ifstream(log_file_full_path.string().c_str()))
    {
            std::cout << "\nERROR: unable to write logfile in working dir." << std::endl;
            return 1;
    }

    //logFile <<  "\n*** " << NAME <<" v. " << SDVERSION  << " - Analysis of 2D/3D multicolor SD-dSTORM data ***\n" << std::endl;
    logFile << ss_welcome.str() << "\n" << std::endl;

    std::string search_pattern_dir;
    std::string search_pattern_file;


    try                                                                 // try to read the config file
    {
        boost::property_tree::ini_parser::read_ini(config_file_full_path.string(), ptree_config_file);
    }
    catch( const std::exception & ex )
    {

        std::cout << "\nERROR: syntax error in config file: " << config_file_full_path.string() << std::endl << ex.what() << std::endl;

        return 1;

    }
    try                                                               // try to load some config-file entries
    {
        search_pattern_dir = ptree_config_file.get<std::string>("General.search_pattern_dir");
        search_pattern_file = ptree_config_file.get<std::string>("General.search_pattern_file");
        PairFinderFileStamp = ptree_config_file.get<std::string>("PairFinder.FileStamp");
        FilterFileStamp = ptree_config_file.get<std::string>("Filter.FileStamp");

        //crap!
        std::string str_grouping = ptree_config_file.get<std::string>("grouping.grouping_on");
        std::string str_groupingFileStamp = ptree_config_file.get<std::string>("grouping.OutputFileStamp");
        //CHANGE ME!!
        int grouping_on = atoi(str_grouping.c_str());
        PairFinderFileStampGroupingFailed = PairFinderFileStamp;        //just in case there are no groups..
        if(grouping_on != 0)
            PairFinderFileStamp = PairFinderFileStamp+str_groupingFileStamp;
    }
    catch( const std::exception & ex )
    {

        std::cout << "\nERROR: syntax error in config file: " << config_file_full_path.string() << std::endl << ex.what() << std::endl;

        return 1;

    }

    boost::regex pattern_dir(search_pattern_dir.c_str(), boost::regex_constants::icase);
    boost::regex pattern_file(search_pattern_file.c_str(), boost::regex_constants::icase);


    if ( fs::is_directory( full_path ) )                            // the working dir should be an directory
    {
        try
        {
            fs::directory_iterator end_iter;

            //std::cout << "\nWorking directory: " << full_path.string() << "\n" << std::endl;


            for ( fs::directory_iterator dir_itr( full_path ); dir_itr != end_iter; ++dir_itr )
            {
                try
                {
                    if ( fs::is_directory( dir_itr->status() ) )        // found a directory
                    {
                        ++dir_count;

                        std::stringstream ssfile;
                        ssfile << dir_itr->path().filename();

                        if (boost::regex_match(ssfile.str(), pattern_dir))      //dir matches the search pattern from config file
                        {                                                        //
                            std::string match_dir = ssfile.str();               //

                            fs::path full_path_matched_dir = operator/(full_path.string(), dir_itr->path().filename());

                            dir_list.push_back(full_path_matched_dir.string());         // add to vector which contains all our matched dirs

                            ++dir_match_count;
                        }

                    }
                    else if ( fs::is_regular_file( dir_itr->status() ) )            // found a subdir in the working dir
                    {
                        ++file_count;

                        std::stringstream ssfile;
                        ssfile << dir_itr->path().filename();

                        if (boost::regex_match(ssfile.str(), pattern_file))             // pattern matching?
                        {

                            std::string match_file = ssfile.str();

                            fs::path full_path_matched_file = operator/(full_path.string(), dir_itr->path().filename());
                            file_list.push_back(full_path_matched_file.string());       // yes, add to vector containing all the files

                            ++file_match_count;
                        }
                    }

                }
                catch ( const std::exception & ex )
                {
                    ++err_count;
                    std::cout << dir_itr->path().filename() << " " << ex.what() << std::endl;
                    return 1;
                }
            }
        }

        catch ( const std::exception & ex )
        {
            std::cout << "ERROR: failed opening directory " << full_path.string() << std::endl << ex.what() << std::endl;
            return 1;
        }
    }
    else
    {
        std::cout << "\nERROR: " << full_path.string() << " is not a directory" << std::endl << std::endl;
        return 1;
    }

    if ( file_match_count == 0 && dir_match_count != 0)
    {

        std::sort(dir_list.begin(), dir_list.end());
        //std::cout << "\nstarting analysis of " << dir_match_count << " subdirectorie(s). lean back, this may take a while...\n\n";
        AnalyseDirs(ptree_config_file, config_file_full_path, force_dim, one_color, dir_list, skip_stage, full_path.string(), end_stage);

    }

    if ( file_match_count != 0 && dir_match_count == 0)
    {
        int nr_of_dir = 0;
        //std::cout << "\nstarting analysis of " << file_match_count << " file(s) in the working directory.\n";
        std::sort(file_list.begin(), file_list.end());
        AnalyseFiles(ptree_config_file, config_file_full_path, force_dim, one_color, file_list, skip_stage, full_path.string(), end_stage, dir_list, nr_of_dir);

    }

    if ( file_match_count != 0 && dir_match_count != 0)
    {
        std::cout << "\nERROR: found " << file_match_count << " file(s) and "
        << dir_match_count << " directorie(s) that matched the search pattern\n" <<
        "please make sure that your working directory contains either subdirectories\nwhich content has to be processed " <<
        "OR files that have to be processed. Not both!\n" <<
        "don't know what you want me to do. terminating\n\n";

        return 1;
    }

    if ( file_match_count == 0 && dir_match_count == 0)
    {
        std::cout << "\nfound no files or directories that match the regex search pattern.\nNothing to do, terminating..\n\n";
        return 0;
    }

}


int AnalyseFiles(boost::property_tree::ptree pt, fs::path config_file_full_path, bool force2D, bool one_color,std::vector<std::string> file_list,
                 int skip_stage, std::string working_dir, int end_stage, std::vector<std::string> input_dirs, int& nr_of_dir )
{

    std::string old_file_name;
    std::string pf_input;
    std::string filter_input;
    std::string recon_input;
    std::string only_filename;

    octave_value_list loadfile, result;

    Matrix FullData;

    bool FullData_loaded = false;
    int file_count = 0;
    int PairFinderStatus;

    if (!input_dirs.empty())
    {
        std::cout << "\ndirectory [" << nr_of_dir << "/" << input_dirs.size() <<"]: "  << input_dirs.at(nr_of_dir-1) << std::endl;
    }


    for (std::vector<std::string>::iterator n = file_list.begin(); n != file_list.end(); ++n)
    {

        //std::cout << *n << std::endl;
        old_file_name = n->c_str();
        if (!std::ifstream(n->c_str()))
        {
            std::cout << "\nWARNING: failed opening file:\n" << *n  << std::endl << "check permissions! skipping...\n" << std::endl;

        }
        else{

    //************************ Here are the demixing stages ********************************
            file_count++;

            only_filename = old_file_name;

            const size_t last_slash_idx = only_filename.find_last_of("\\/");
            if (std::string::npos != last_slash_idx)
            {
                only_filename.erase(0, last_slash_idx + 1);
            }


            std::cout << "\ndemixing file [" << file_count << "/" << file_list.size() << "]: " << only_filename  << std::endl;
            logFile << "file: " << old_file_name << "\n" << std::endl;
            if (skip_stage < 1)
            {
                logFile << "\nPairFinder:\n\n" ;
                PairFinderStatus = PairFinderMain(&pt, &config_file_full_path, old_file_name);

                if (PairFinderStatus == 1)
                {
                    std::cout << "\nWARNING: PairFinder failed" << std::endl;
                    logFile << "\nPairFinder failed\n" ;
                    return 1;               // PF failed, want to continue with the next file, so we exit here
                }
                if (PairFinderStatus == 2)
                {
                    PairFinderFileStamp=PairFinderFileStampGroupingFailed;
                }

                if (end_stage == 1)
                {
                    logFile << "exit stage : " << end_stage << " terminating after PairFinder" << std::endl;
                    //return 0;
                }

            }
            else{
            logFile << "skipping PairFinder...\n";
            }
            filter_input = n->append(PairFinderFileStamp);
            recon_input = n->append(FilterFileStamp);
            if (skip_stage < 2 && end_stage != 1)
            {
                logFile << "\nFilter:\n\n" ;


                std::filebuf fb;
                fb.open (filter_input.c_str(),std::ios::in);
                std::istream is(&fb);

                octave_value temp;
                read_mat_ascii_data (is, filter_input, temp);

                FullData = temp.matrix_value();

                FullData_loaded = true;

                if (FilterMain(&pt, &config_file_full_path, filter_input, working_dir, old_file_name, FullData, one_color) == 1)
                {
                    std::cout << "\nWARNING: Filter failed" << std::endl;
                    logFile << "\nFilter failed\n" ;
                    return 1;               // same as above
                }
                if (end_stage == 2)
                {
                    logFile << "exit stage : " << end_stage << " terminating after Filter" << std::endl;
                    //return 0;
                }

            }
            else{
            logFile << "skipping Filter...\n";
            }
            if (skip_stage < 3 && end_stage != 1 && end_stage != 2)
            {
                if(!FullData_loaded)
                {
                    std::filebuf fb;
                    fb.open (recon_input.c_str(),std::ios::in);
                    std::istream is(&fb);

                    octave_value temp;
                    read_mat_ascii_data (is, recon_input, temp);

                    FullData = temp.matrix_value();

                    fb.close();
                }

                logFile << "\nReconstructor:\n\n" ;
                if (ReconstructorMain(&pt, &config_file_full_path, recon_input.c_str(), working_dir, old_file_name.c_str(), force2D, FullData) == 1)
                {
                    std::cout << "\nWARNING: Reconstructor failed" << std::endl;
                    logFile << "\nReconstructor failed\n" ;
                    return 1;                   // is this futile?
                }

            }
            else{
            logFile << "skipping Reconstructor...\n";
            }
            if (n != (file_list.end() - 1))
            {
                //std::cout << "\ndemixing finished. "  << (file_list.end() - n)-1 << " files in directory remaining"  << std::endl;
                logFile << "\n*********************** demixing finished, next file\n" << std::endl;
            }else if (nr_of_dir != 0 )
            {
                try{
                logFile << "\n*****demixing finished. next directory: " << input_dirs.at(nr_of_dir) << "\n" << std::endl;
                }
                catch(std::exception e){}
            }



    //************************ Here is the end of the demixing stages ********************************
        }
    }
    return 0;

}


int AnalyseDirs(boost::property_tree::ptree pt, fs::path config_file_full_path, bool force2D, bool one_color, std::vector<std::string> inputdirs,
                int skip_stage, std::string working_dir, int end_stage)
{
    std::string search_pattern_file;
    int file_count;
    int file_match_count;
    int err_count;
    fs::path full_path( fs::initial_path<fs::path>() );
    int dir_count = 0;

    for (std::vector<std::string>::iterator n = inputdirs.begin(); n != inputdirs.end(); ++n)
    {
        dir_count++;
        file_count = 0;
        file_match_count = 0;
        err_count = 0;
        full_path = fs::system_complete( fs::path(n->c_str()) );
        std::vector<std::string> file_list;

        if (!std::ifstream(n->c_str()))
        {
            std::cout << "\nWARNING: failed entering directory:\n" << *n  << std::endl << "check permissions! skipping...\n" << std::endl;
            //return 1;
        }
        else
        {
            try
            {
                search_pattern_file = pt.get<std::string>("General.search_pattern_file");
            }
            catch( const std::exception & ex )
            {

                std::cout << "\nERROR: syntax error in config file: " << config_file_full_path.string() << std::endl << ex.what() << std::endl;

                return 1;

            }

            boost::regex pattern_file(search_pattern_file.c_str(), boost::regex_constants::icase);

            try
            {
                fs::directory_iterator end_iter;

                for ( fs::directory_iterator dir_itr( n->c_str() ); dir_itr != end_iter; ++dir_itr )
                {
                    try
                    {
                        if ( fs::is_directory( dir_itr->status() ) )        //ignoring subdirs
                            {


                            }
                            else if ( fs::is_regular_file( dir_itr->status() ) )
                            {
                                ++file_count;

                                std::stringstream ssfile;
                                ssfile << dir_itr->path().filename();

                                if (boost::regex_match(ssfile.str(), pattern_file))
                                {

                                    std::string match_file = ssfile.str();

                                    fs::path full_path_matched_file = operator/(full_path.string(), dir_itr->path().filename());
                                    file_list.push_back(full_path_matched_file.string());

                                    ++file_match_count;
                                }
                            }

                        }
                        catch ( const std::exception & ex )
                        {
                                    ++err_count;
                                    std::cout << dir_itr->path().filename() << " " << ex.what() << std::endl;
                                    return 1;
                        }
                    }
                }

                catch ( const std::exception & ex )
                {
                        std::cout << "ERROR: failed opening directory " << full_path.string() << std::endl << ex.what() << std::endl;
                        //return 1;
                }

        }

        if ( file_match_count != 0)
        {

            //std::cout << "\nstarting analysis of " << file_match_count << " file(s) in the subdirectory: " << n->c_str() << "\n\n";
            std::sort(file_list.begin(), file_list.end());
            AnalyseFiles(pt, config_file_full_path, force2D, one_color, file_list, skip_stage, working_dir, end_stage, inputdirs, dir_count);

        }
    }


    return 0;
}

std::string get_date(void)
{
   time_t now;
   const int MAX_DATE=16;
   char the_date[MAX_DATE];

   the_date[0] = '\0';

   now = time(NULL);

   if (now != -1)
   {
      strftime(the_date, MAX_DATE, "%Y%m%d_%H%M%S", localtime(&now));
   }

   return std::string(the_date);
}

int ReturnDimensions( const std::string & pattern, const std::string & firstLine )
{
    int n = 0;
    std::string ::size_type pos = 0;
    while( (pos = firstLine.find( pattern, pos ))
                 != std::string::npos ) {
    	n++;
    	pos += pattern.size();
    }
    return n;
}


std::vector<double> readMaxValues( std::string xmlConfigString)
{
    // populate tree structure pt
    using boost::property_tree::ptree;
    ptree pt;
    // Read the XML config string into the property tree. Catch any exception
    xmlConfigString.erase (xmlConfigString.begin());
    //std::cout << xmlConfigString << std::endl;
    std::stringstream ss; ss << xmlConfigString;
    read_xml(ss, pt);

    std::vector<double> maxValueList;

    // traverse pt
    loc ans;
    BOOST_FOREACH( ptree::value_type const& v, pt.get_child("localizations") )
    {
        if( v.first == "field" )
        {
            try
            {
                Localizations f;
                double value;

                f.identifier = v.second.get<std::string>("<xmlattr>.identifier");
                //f.syntax = v.second.get<std::string>("<xmlattr>.syntax");
                //f.semantic =v.second.get<std::string>("<xmlattr>.semantic");
                //f.unit = v.second.get<std::string>("<xmlattr>.unit");
                //f.min = v.second.get<std::string>("<xmlattr>.min");
                f.max = v.second.get<std::string>("<xmlattr>.max");
                std::istringstream iss(f.max);
                iss >> value;
                maxValueList.push_back(value);
                //std::cout << f.max << std::endl;
                ans.push_back(f);
            }
            catch(boost::exception const& ex){}
        }
    }

   /* for (int i = 0; i < maxValueList.size(); i++)
    {
        std::cout << maxValueList[i] << std::endl;
    }*/
    return maxValueList;
}

