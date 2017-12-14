////////////////////////////////////////////////////////////////////////////////
// File: simpsons_rule_tab.c                                                  //
// Routines:                                                                  //
//    Simpsons_Rule_Tab_Sum_LR                                                //
//    Simpsons_Rule_Tab_Sum_RL                                                //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     A (closed) Newton-Cotes formula is an interpolatory quadrature formula //
//     using equally spaced nodes for approximating the integral of a         //
//     function over a closed and bounded interval.                           //
//                                                                            //
//     Simpson's rule is a (closed) Newton-Cotes formula for which the        //
//     interpolating polynomial is a quadratic, the parabola passing through  //
//     (x,f(x)) at the endpoints of the interval and through (xm,f(xm))       //
//     where xm is the midpoint of the interval.                              //
//     The estimation of the integral then is simply the area beneath the     //
//     parabola formed by the line segment joining the points (x, f(x)) at the//
//     endpoints, the line segments joining (x,0) to (x,f(x)) at the endpoints//
//     and the interval of integration.                                       //
//                                                                            //
//     In order to integrate a function using Simpson's rule over a closed    //
//     amd bounded interval [a,b], divide the interval [a,b] into N           //
//     subintervals of length h = (b-a) / N.  The integral of the function    //
//     over the interval [a,b] is the sum of the integrals of the function    //
//     over the subintervals [a,a+h], [a+h, a+2h],...,[b-h,b].  The integral  //
//     is approximated then by applying the Simpson's rule to each            //
//     subinterval.  By abuse of language, this composite Simpson's rule      //
//     is simply also called Simpson's rule.                                  //
//                                                                            //
//     Let I(f) be the approximation of the integral of f(x) over [a,b] then  //
//     the Simpson's rule is:                                                 //
//       I(f) = h/6 { f(a) + 4f(a+h/2) + 2f(a+h) + 4f(a+3h/2) + ... +         //
//                                                       4f(b-h/2) + f(b)) }  //
//     where h is the subinterval length, h = (b - a) / n, n is the number of //
//     subintervals.                                                          //
//     The truncation error estimate for Simpson's rule is                    //
//                     - n * h^5  f''''(xm) / (90 * 32)                       //
//     where f''''(xm) is the value of the fourth derivative evaluated at     //
//     some (unknown) point a <= xm <= b.                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//  double Simpsons_Rule_Tab_Sum_LR( double h, int n, double (*f)(double) );  //
//                                                                            //
//  Description:                                                              //
//     Assuming that a is the lower limit of integration and the ith element  //
//     of the array f[] is the value of the function evaluated at the left    //
//     endpoint of the ith subinterval for 0 <= i < n and f[n] is the right   //
//     endpoint of the nth interval, i.e. f[i] = f( a + i*h ), i = 0,...,n,   //
//     then  this routine integrates f(x) using Simpson's rule by summing     //
//     from the left end of the interval to the right end.                    //
//                                                                            //
//  Arguments:                                                                //
//     double h   The length of each subinterval, h > 0.                      //
//     int    n   The number of subintervals, the upper limit of integration  //
//                is a + n * h, where a is the lower limit of integration.    //
//     double f[] Array of function values the ith element of which is the    //
//                f[i] = f(a + i*h/2), i = 0,...,2n, i.e. the dimension of f  //
//                is 2n + 1.                                                  //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) from a to (a + n * h).                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Simpsons_Rule_Tab_Sum_LR( double h, int n, double f[] ) {
 
   double integral = f[0] + 4.0 * f[1];
   int i;

   for (i = 2; i < (n + n); i += 2) 
      integral += 2.0 * f[i] + 4.0 * f[i+1];
 
   return 0.1666666666666666666667 *  h * ( integral + f[i] );
}


////////////////////////////////////////////////////////////////////////////////
//  double Simpsons_Rule_Tab_Sum_RL( double h, int n, double f[] );           //
//                                                                            //
//  Description:                                                              //
//     Assuming that a is the lower limit of integration and the ith element  //
//     of the array f[] is the value of the function evaluated at the left    //
//     endpoint of the ith subinterval for 0 <= i < n and f[n] is the right   //
//     endpoint of the nth interval, i.e. f[i] = f( a + i*h ), i = 0,...,n,   //
//     then  this routine integrates f(x) using Simpson's rule by summing     //
//     from the left end of the interval to the right end.                    //
//                                                                            //
//  Arguments:                                                                //
//     double h   The length of each subinterval, h > 0.                      //
//     int    n   The number of subintervals, the upper limit of integration  //
//                is a + n * h, where a is the lower limit of integration.    //
//     double f[] Array of function values the ith element of which is the    //
//                f[i] = f(a + i*h/2), i = 0,...,2n, i.e. the dimension of f  //
//                is 2n + 1.                                                  //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) from a to (a + n * h).                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Simpsons_Rule_Tab_Sum_RL( double h, int n, double f[] ) {
 
   double integral  = f[n+n] + 4.0 * f[n+n-1];
   int i;

   for (i = n+n-2; i > 0; i -= 2) integral += 2.0 * f[i] + 4.0 * f[i-1];
 
   return 0.1666666666666666666667 *  h * ( integral + f[0] );
}