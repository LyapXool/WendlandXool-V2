/* WendlandXool, version 2, is a program to construct and evaluate Wendland functions,
 ->
 ->
 -> This program is free software; you can redistribute it and/or
 -> modify it under the terms of the GNU General Public License
 -> as published by the Free Software Foundation; either version 3
 -> of the License, or (at your option) any later version.
 ->
 -> This program is distributed in the hope that it will be useful,
 -> but WITHOUT ANY WARRANTY; without even the implied warranty of
 -> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 -> GNU General Public License for more details.
 ->
 -> You should have received a copy of the GNU General Public License
 -> along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ->
 -> Authors: Carlos Argáez, Peter Giesl, Sigurdur Freyr Hafstein
 */

#ifndef wendland_hpp
#define wendland_hpp

#include <stdio.h>
#include <armadillo>
class WendRBF {
	arma::Mat<double> wdlfunction, wdlf1, wdlf2;
	arma::Mat<double> wdlfunctionfac, wdlf1fac, wdlf2fac;
	int powl0, powl1, powl2;
	double c;
public:
	WendRBF(int l, int k, double _c,bool _printreport=false);
	double operator()(double r);
	double aux1(double r);
	double aux2(double r);
};

namespace wendland{

    void wendlandfunction(arma::Mat<double> &wdlf,int l,int k);

    void wendlandderivative(arma::Mat<double> &wdlfinput, int l, int k, arma::Mat<double> &wdlf1, int wendlandorder);

    double evawdlfn(double r, double c, arma::Mat<double> const &wdlfn, int wendlandorder);

    double sdsds(double r, double c, arma::Mat<double> const &wdlfn, int wendlandorder);

    double sdsdssdsds(double r, double c, arma::Mat<double> const &wdlfn, int wendlandorder);

    void pascal(int currentdimension, int l, arma::Mat <double> &vector1);

    void tprod(arma::Mat <double> &cpower);

    void fixindex(arma::Mat <double> &vector);

    void printwendland(int j, arma::Mat <double> &vector1, int wendlandorder);

    void printwendlandr(bool finalization, arma::Mat <double> &vector1, int wendlandorder);

    void printwendlandf(arma::Mat <double> &vector1, int wendlandorder);

    void printwendlandff(arma::Mat <double> &vector1, int wendlandorder);

    char signaturing(int i, double unknownsign);

    void cleanzeros(arma::Mat<double> &vector);

    void writeequation(const int classifier);

    void mcmp(long long unsigned value1, long long unsigned value2, long long unsigned &valueout);

    void gcmp(long long unsigned value1, long long unsigned value2, long long unsigned &valueout);

    long long unsigned getmcm(arma::Mat<double> &vector);

    long long unsigned getmcd(arma::Mat<double> &vector);

	bool syndiv1(int l, int k, const arma::Mat<double> &wdlin, arma::Mat<double> &wdlout, int *pow1mr, int order);

    void startsection(std::string valor);

    void startreport(int wendlandorder);

    extern std::ofstream wendlandreport;
    
    extern bool factorise;
};

#endif /* wendland_hpp */
