// CIFlow is a very flexible configuration interaction program
// Copyright (C) Ghent University 2014-2015
//
// This file is part of CIFlow.
//
// CIFlow is developed by Mario Van Raemdonck <mario.vanraemdonck@ugent.be>
// a member of the Ghent Quantum Chemistry Group (Ghent University).
// See also : http://www.quantum.ugent.be
//
// At this moment CIFlow is not yet distributed.
// However this might change in the future in the hope that
// it will be useful to someone.
//
// For now you have to ask the main author for permission.
//
//--
#include <math.h>
#include <stdio.h>

#include "binomial.h"


/*Algorithm to calculate the binomial coefficient, found in numerical recipes for C. */

double gammln(double xx) /*    Returns the value ln[Γ( xx)] for xx > 0 .*/
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

double factln(int n)
{
    static double a[101];
    if (n < 0) printf("Negative factorial in routine factln");    /*A static array is automatically initialized to zero. */
    if (n <= 1) return 0.0;
    if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
    else return gammln(n+1.0);
}

double bico(int n, int k){  /*Returns the binomial coefficient n k as a floating-point number. */
	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}
