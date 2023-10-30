//!#####################################################################
//! \file Utilities.cpp
//!#####################################################################
#include "Utilities.h"
#include <cassert>
#include <fstream>
#include <random>
//######################################################################
// Number_Of_Intervals
//######################################################################
unsigned long Utilities::
Number_Of_Intervals(int argc, char** argv)
{
    unsigned long nentries = 0;

    for(int i=1;i<argc;++i){
        FILE* fp = fopen(argv[i], "r");
        assert(fp != nullptr);

        fseek(fp, 0, SEEK_END);
        unsigned long current_nentries = ftell(fp);
        current_nentries /= (3*sizeof(double));
        printf("number of intervals in file %d = %lu\n", i, current_nentries);
        nentries += current_nentries;
        fclose(fp);
    }
    return nentries;
}
//######################################################################
// Read_Intervals_From_Files
//######################################################################
void Utilities::
Read_Intervals_From_Files(int argc,char** argv,interval_data* intervals)
{
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    unsigned long j=0;
    for(int i=1;i<argc;++i){
        FILE* fp = fopen(argv[i], "r");
        assert(fp != nullptr);

        fseek(fp, 0, SEEK_END);
        unsigned long nentries = ftell(fp);
        nentries /= (3*sizeof(double));
        fseek(fp, 0, SEEK_SET);

        for(unsigned long k = 0; k < nentries; ++k){
            double data_entry[3];
            size_t bytes = fread(data_entry, sizeof(double), 3, fp);
            intervals[j] = interval_data(data_entry[0],data_entry[1],data_entry[2],
                                         1.0,distribution(generator),i-1);
            ++j;
        }
        fclose(fp);
    }
}
//######################################################################
