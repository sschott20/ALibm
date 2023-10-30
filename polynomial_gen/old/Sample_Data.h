//!#####################################################################
//! \file Sample_Data.h
//!#####################################################################
// Class Sample_Data
//######################################################################
#ifndef __Sample_Data__
#define __Sample_Data__

class sample_data
{
  public:
    double x;
    double lb;
    double ub;
    double orig_lb;  /* remember lb before narrowing */
    double orig_ub;  /* remember ub before narrowing */
    double w;
    double u;
    double k;        /* key computed as 1/u^w */
    int i;

    sample_data() {}
    ~sample_data() {}
};
#endif
