/********************************************************************************
	BLOCKSIM é um programa limitado que realiza algumas simulações de 
	sistemas dinâmicos que podem ser representados por diagramas de blocos.
    Copyright (C) 2004, Renato Mikio Nakagomi

    Este arquivo é parte do programa BLOCKSIM. 

	As funções aqui apresentadas neste arquivo foram encontradas em DOMÍNIO PÚBLICO
	no site http://www.mymathlib.webtrellis.net/index.html
	No momento estamos tentando entrar em contato com o autor pois, apesar do 
	código-fonte destas funções estarem em DOMÍNIO PÚBLICO, acreditamos que devemos 
	honrar nossos colaboradores, prestando a devida homenagem e respeitando a
	inteligência alheia.

    Este programa é software livre; você pode redistribuí-lo e/ou
    modificá-lo sob os termos da Licença Pública Geral GNU, conforme
    publicada pela Free Software Foundation; tanto a versão 2 da
    Licença como (a seu critério) qualquer versão mais nova.

    Este programa é distribuído na expectativa de ser útil, mas SEM
    QUALQUER GARANTIA; sem mesmo a garantia implícita de
    COMERCIALIZAÇÃO ou de ADEQUAÇÃO A QUALQUER PROPÓSITO EM
    PARTICULAR. Consulte a Licença Pública Geral GNU para obter mais
    detalhes.
 
    Você deve ter recebido uma cópia da Licença Pública Geral GNU
    junto com este programa; se não, escreva para a Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
    02111-1307, USA.

********************************************************************************
	BLOCKSIM is a limited software that runs some simulations of dinamic systems
	that can be represented by block diagrams.
    Copyright (C) 2004, Renato Mikio Nakagomi

    This file is part of BLOCKSIM software. 

	The functions presented in this file were found in Public Domain at site
	http://www.mymathlib.webtrellis.net/index.html
	Although the source-code of those functions are in Public Domain,
	we are trying to contact the author because we believe we must honour our 
	colaborators, respecting their intelligence.

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

********************************************************************************
	Renato Mikio Nakagomi
	renato.nakagomi@poli.usp.br
********************************************************************************/


////////////////////////////////////////////////////////////////////////////////
// File: trapezoidal_rule_tab.c                                               //
// Routines:                                                                  //
//    Trapezoidal_Rule_Tab_Sum_LR                                             //
//    Trapezoidal_Rule_Tab_Sum_RL                                             //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     A Newton-Cotes formula is an interpolatory quadrature formula using    //
//     the values of the integrand (function being integrated) at equally     //
//     spaced nodes to approximate the integral of the function over a        //
//     closed and bounded interval.                                           //
//                                                                            //
//     The trapezoidal rule is a (closed) Newton-Cotes formula for which the  //
//     interpolating polynomial is a line, the line joining the points        //
//     (x,f(x)) at the endpoints of the interval.                             //
//     The estimation of the integral then is simply the area of the trapezoid//
//     formed by the line segment joining the points (x, f(x)) at the         //
//     endpoints, the line segments joining (x,0) to (x,f(x)) at the endpoints//
//     and the interval of integration.                                       //
//                                                                            //
//     In order to integrate a function using the trapezoidal rule over a     //
//     closed and bounded interval [a,b], divide the interval [a,b] into N    //
//     subintervals of length h = (b-a) / N.  The integral of the function    //
//     over the interval [a,b] is the sum of the integrals of the function    //
//     over the subintervals [a,a+h], [a+h, a+2h],...,[b-h,b].  The integral  //
//     is approximated then by applying the trapezoidal rule to each          //
//     subinterval.  By abuse of language, this composite trapezoidal rule    //
//     is simply also called the trapezoidal rule.                            //
//                                                                            //
//     Let I(f) be the approximation of the integral of f(x) over [a,b] then  //
//     the trapezoidal rule is:                                               //
//       I(f) = h { f(x[0]) / 2 + f(x[1]) + ... + f(x[n-1]) + f(x[n]) / 2) }  //
//     where h is the step size, h = x[i] - x[i-1]; n is the number of        //
//     steps, n = abs(b - a) / h and x[0] = a, x[i] = x[i-1] + h for          //
//     i = 1,...,n-1 so that x[n] = b.                                        //
//                                                                            //
//     The truncation error estimate for the trapezoidal rule is              //
//                           n * h^3   f''(xm) / 12                           //
//     where f''(xm) is the value of the second derivative evaluated at some  //
//     (unknown) point a <= xm <= b.                                          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//  double Trapezoidal_Rule_Tab_Sum_LR( double h, int n, double f[] );        //
//                                                                            //
//  Description:                                                              //
//     Assuming that a is the lower limit of integration and the ith element  //
//     of the array f[] is the value of the function evaluated at the left    //
//     endpoint of the ith subinterval for 0 <= i < n and f[n] is the right   //
//     endpoint of the nth interval, i.e. f[i] = f( a + i*h ), i = 0,...,n,   //
//     then  this routine integrates f(x) using the trapezoidal rule by       //
//     summing from the left end of the interval to the right end.            //
//                                                                            //
//  Arguments:                                                                //
//     double h   The length of each subinterval, h > 0.                      //
//                a + n * h, if a is the lower limit of integration.          //
//     int    n   The number of subintervals, the upper limit of integration  //
//                is a + n * h, where a is the lower limit of integration.    //
//     double f[] Array of function values the ith element of which is the    //
//                f[i] = f(a + i*h), i = 0,...,n, i.e. the dimension of f is  //
//                n + 1.                                                      //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) from a to (a + n * h).                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
long double Trapezoidal_Rule_Tab_Sum_LR( double h, int n, double f[] ) { 
//float Trapezoidal_Rule_Tab_Sum_LR( float h, int n, float f[] ) { 

   double integral = 0.5 * f[0];
//   float integral = 0.5 * f[0];
   int i;
   
	//printf(" f[0] = %lf \n",f[0]);

   for (i = 1; i < n; i++) {
	   integral += f[i];
		//printf(" f[%d] = %lf \n",i,f[i]);
	   
   }

   return h * ( integral + 0.5 * f[n] );
}


////////////////////////////////////////////////////////////////////////////////
//  double Trapezoidal_Rule_Tab_Sum_RL( double h, int n, double f[] );        //
//                                                                            //
//  Description:                                                              //
//     Assuming that a is the lower limit of integration and the ith element  //
//     of the array f[] is the value of the function evaluated at the left    //
//     endpoint of the ith subinterval for 0 <= i < n and f[n] is the right   //
//     endpoint of the nth interval, i.e. f[i] = f( a + i*h ), i = 0,...,n,   //
//     then  this routine integrates f(x) using the trapezoidal rule by       //
//     summing from the right end of the interval to the left end.            //
//                                                                            //
//  Arguments:                                                                //
//     double h   The length of each subinterval, h > 0.                      //
//     int    n   The number of subintervals, the upper limit of integration  //
//                is a + n * h, where a is the lower limit of integration.    //
//     double f[] Array of function values the ith element of which is the    //
//                f[i] = f(a + i*h), i = 0,...,n, i.e. the dimension of f is  //
//                n + 1.                                                      //
//                                                                            //
//  Return Values:                                                            //
//     The integral of f(x) from a to (a + n * h).                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
double Trapezoidal_Rule_Tab_Sum_RL( double h, int n, double f[] ) { 
//float Trapezoidal_Rule_Tab_Sum_RL( float h, int n, float f[] ) { 

   double integral = 0.5 * f[n];
//   float integral = 0.5 * f[n];
   int i;

   for (i = n - 1; i > 0; i--) integral += f[i];

   return h * ( integral + 0.5 * f[0] );
}