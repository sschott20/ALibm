//!#####################################################################
//! \file Interval_Data.h
//!#####################################################################
// Class Interval_Data
//######################################################################
#ifndef __Interval_Data__
#define __Interval_Data__

class interval_data
{
  public:
    double x;         /* original input */
    double lb;        /* lower bound */
    double ub;        /* upper bound */
    double w;         /* weight */
    double u;         /* uniform random value */
    int i;            /* index in the powers table */

    interval_data() {}

    interval_data(const double& x_input,const double& lb_input,const double& ub_input,
                  const double& w_input,const double& u_input,const int i_input)
        :x(x_input),lb(lb_input),ub(ub_input),w(w_input),u(u_input),i(i_input)
    {}

    ~interval_data() {}

    interval_data& operator=(const interval_data& rhs)
    {
        if(this==&rhs) return *this;

        x = rhs.x;
        lb = rhs.lb;
        ub = rhs.ub;
        w = rhs.w;
        u = rhs.u;
        i = rhs.i;
        return *this;
    }
};
#endif
