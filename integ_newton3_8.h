////////////////////////////////////////////////////////////////////////////////
// File: newtons_3_8_rule_tab.c                                               //
// Routines:                                                                  //
//    Newtons_3_8_Rule_Tab_Sum_LR                                             //
//    Newtons_3_8_Rule_Tab_Sum_RL                                             //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     A (closed) Newton-Cotes formula is an interpolatory quadrature formula //
//     using equally spaced nodes for approximating the integral of a         //
//     function over a closed and bounded interval.                           //
//                                                                            //
//     Newton's 3/8ths rule is a (closed) Newton-Cotes formula for which the  //
//     interpolating polynomial is a cubic, the cubic passing through         //
//     (x,f(x)) at the endpoints of the interval and through (x1,f(x1)) and   //
//     (x2,f(x2)) where x1 is 1/3 of the distance from to left endpoint to    //
//     the right endpoint and x2 is 1/3 of the distance from the right enpoint//
//     to the left endpoint.                                                  //
//     The estimation of the integral then is simply the area beneath the     //
//     cubic from the left endpoint to the right endpoint.                    //
//                                                                            //
//     In order to integrate a function using Newton's 3/8ths rule over a     //
//     closed and bounded interval [a,b], divide the interval [a,b] into N    //
//     subintervals of length h = (b-a) / N.  The integral of the function    //
//     over the interval [a,b] is the sum of the integrals of the function    //
//     over the subintervals [a,a+h], [a+h, a+2h],...,[b-h,b].  The integral  //
//     is approximated then by applying the Newton's 3/8ths rule to each      //
//     subinterval.  By abuse of language, this composite Newton's 3/8ths     //
//     rule is simply also called Newton's 3/8ths rule.                       //
//                                                                            //
//     Let I(f) be the approximation of the integral of f(x) over [a,b] then  //
//     Newton's 3/8ths rule is:                                               //
//       I(f) = h/8 { f(a) + 3f(a+h/3) + 3f(a+2h/3) + 2f(a+h) + ... +         //
//                                          3f(b-2h/3) + 3 f(b-h/3) + f(b) }  //
//     where h is the subinterval length, h = (b - a) / n, n is the number of //
//     subintervals.                                                          //
//                                                                            //
//     The truncation error estimate for Newton's 3/8ths rule is              //
//                    -  * n * h^5  f''''(xm) / 6480                          //
//     where f''''(xm) is the value of the fourth derivative evaluated at     //
//     some (unknown) point a <= xm <= b.                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//  double Newtons_3_8_RuleTab_Sum_LR( double h, int n, double f[] );         //
//                                                                            //
//  Description:                                                              //
//     Assuming that a is the lower limit of integration and the ith element  //
//     of the array f[] is the value of the function evaluated at the left    //
//     endpoint of the ith subinterval for 0 <= i < n and f[n] is the right   //
//     endpoint of the nth interval, i.e. f[i] = f( a + i*h ), i = 0,...,n,   //
//     then  this routine integrates f(x) using Newton's 3/8 rule by summing  //
//     from the left end of the interval to the right end.                    //
//                                                                            //
//  Arguments:                                                                //
//     double h   The length of each subinterval, h > 0.                      //
//     int    n   The number of subintervals, the upper limit of integration  //
//                is a + n * h, where a is the lower limit of integration.    //
//     double f[] Array of function values the ith element of which is the    //
//                f[i] = f(a + i*h/3), i = 0,...,3n, i.e. the dimension of f  //
//                is 3n + 1.                                                  //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) from a to (a + n * h).                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
double Newtons_3_8_Rule_Tab_Sum_LR( double h, int n, double f[] ) {
 
   double integral = f[0] + 3.0 * ( f[1] + f[2] );
   int i;
   
   for (i = 3; i < ( n + n + n ); i += 3)
      integral += 2.0 * f[i] + 3.0 * ( f[i+1] + f[i+2] );

   return 0.125 * h * ( integral + f[i] );
}


////////////////////////////////////////////////////////////////////////////////
//  double Newtons_3_8_RuleTab_Sum_RL( double h, int n, double f[] );         //
//                                                                            //
//  Description:                                                              //
//     Assuming that a is the lower limit of integration and the ith element  //
//     of the array f[] is the value of the function evaluated at the left    //
//     endpoint of the ith subinterval for 0 <= i < n and f[n] is the right   //
//     endpoint of the nth interval, i.e. f[i] = f( a + i*h ), i = 0,...,n,   //
//     then  this routine integrates f(x) using Newton's 3/8 rule by summing  //
//     from the right end of the interval to the left end.                    //
//                                                                            //
//  Arguments:                                                                //
//     double h   The length of each subinterval, h > 0.                      //
//     int    n   The number of subintervals, the upper limit of integration  //
//                is a + n * h, where a is the lower limit of integration.    //
//     double f[] Array of function values the ith element of which is the    //
//                f[i] = f(a + i*h/3), i = 0,...,3n, i.e. the dimension of f  //
//                is 3n + 1.                                                  //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) from a to (a + n * h).                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
double Newtons_3_8_Rule_Tab_Sum_RL( double h, int n, double f[] ) { 
  
   int i = n + n + n; 
   double integral = f[i] + 3.0 * ( f[i-1] + f[i-2] );

   for (i -= 3 ; i > 0; i -= 3) 
      integral += 2.0 * f[i] + 3.0 * ( f[i-1] + f[i-2] );

   return 0.125 * h * ( integral + f[0] );
}