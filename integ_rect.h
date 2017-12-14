////////////////////////////////////////////////////////////////////////////////
// File: rectangle_rule_tab.c                                                 //
// Routines:                                                                  //
//    Rectangle_Rule_Tab_Sum_LR                                               //
//    Rectangle_Rule_Tab_Sum_RL                                               //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     A Newton-Cotes formula is an interpolatory quadrature formula using    //
//     the values of the integrand (function being integrated) at equally     //
//     spaced nodes to approximate the integral of the function over a        //
//     closed and bounded interval.                                           //
//                                                                            //
//     The rectangle rule is an (open) Newton-Cotes formula with a single     //
//     node at the midpoint of the closed and bounded interval.  The          //
//     interpolating polynomial is a constant, the constant being the value   //
//     of the function evaluated at the node.                                 //
//     The estimation of the integral then is simply the function evaluated   //
//     at the mid-point times the length of the interval.                     //
//                                                                            //
//     In order to integrate a function using the rectangle rule over a       //
//     closed and bounded interval [a,b], divide the interval [a,b] into N    //
//     subintervals of length h = (b-a) / N.  The integral of the function    //
//     over the interval [a,b] is the sum of the integrals of the function    //
//     over the subintervals [a,a+h], [a+h, a+2h],...,[b-h,b].  The integral  //
//     is approximated then by applying the rectangle rule to each            //
//     subinterval.  By abuse of language, this composite rectangle rule      //
//     is simply also called the rectangle rule.                              //
//                                                                            //
//     Let I(f) be the approximation of the integral of f(x) over [a,b] then  //
//     the rectangle rule is:                                                 //
//       I(f) = h { f( [x[0] + x[1]] / 2) + ... + f( [x[n-1] + x[n]] / 2) }   //
//     where h is the step size, h = x[i] - x[i-1]; n is the number of        //
//     steps, n = abs(b - a) / h and x[0] = a, x[i] = x[i-1] + h for          //
//     i = 1,...,n-1 so that x[n] = b.                                        //
//                                                                            //
//     For a function f(x) which is at least be differentiable of class C2    //
//     on the interval of integration the truncation error estimate is        //
//                            n * h^3 * f''(xm) / 24                          //
//     where f''(xm) is the value of the second derivative evaluated at some  //
//     (unknown) point, a <= xm <= b.                                         //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//  double Rectangle_Rule_Tab_Sum_LR( double h, int n, double f[] )           //
//                                                                            //
//  Description:                                                              //
//     Assuming that a is the lower limit of integration and the ith element  //
//     of the array f[] is the value of the function evaluated at the midpoint//
//     of the ith subinterval, i.e. f[i] = f( a + (i+0.5)*h ), this routine   //
//     integrates f(x) using the rectangle rule by summing from the left      //
//     end of the interval to the right end.                                  //
//                                                                            //
//  Arguments:                                                                //
//     double h   The length of each subinterval, h > 0.                      //
//     int    n   The number of subintervals, the upper limit of integration  //
//                is a + n * h, where a is the lower limit of integration.    //
//     double f[] Array of function values the ith element of which is the    //
//                function evaluated at a + (i + 0.5)*h.                      //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) from a to (a + n * h).                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
double Rectangle_Rule_Tab_Sum_LR( double h, int n, double f[] ) { 

   double integral = f[0];
   int i;

   for (i = 1; i < n; i++) integral += f[i];
 
   return h * integral;
}


////////////////////////////////////////////////////////////////////////////////
//  double Rectangle_Rule_Tab_Sum_RL( double h, int n, double f[] )           //
//                                                                            //
//  Description:                                                              //
//     Assuming that a is the lower limit of integration and the ith element  //
//     of the array f[] is the value of the function evaluated at the midpoint//
//     of the ith subinterval, i.e. f[i] = f( a + (i+0.5)*h ), this routine   //
//     integrates f(x) using the rectangle rule by summing from the right     //
//     end of the interval to the left end.                                   //
//                                                                            //
//  Arguments:                                                                //
//     double h   The length of each subinterval, h > 0.                      //
//     int    n   The number of subintervals, the upper limit of integration  //
//                is a + n * h.                                               //
//     double f[] Array of function values the ith element of which is the    //
//                function evaluated at a + (i + 0.5)*h.                      //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) from a to (a + n * h).                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
double Rectangle_Rule_Tab_Sum_RL( double h, int n, double f[] ) { 

   double integral = f[n-1];
   int i;

   for (i = n - 2; i >= 0; i--) integral += f[i];
 
   return h * integral;
}