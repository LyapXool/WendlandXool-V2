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

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <iterator>
#include <armadillo>
#include <algorithm>
#include <iomanip>      // std::setprecision
#include "wendland.hpp"



using namespace std;
using namespace arma;
span const All=span::all;
bool printreport=false;

std::ofstream wendland::wendlandreport;
bool wendland::factorise;




WendRBF::WendRBF(int l, int k, double _c, bool _printreport) : c(_c) {
    printreport = _printreport;
    powl0 = powl1 = powl2 = 0;
    if (printreport == true) {
        wendland::wendlandreport.open("wendlandreport.tex", fstream::out);
        wendland::startreport(0);
    }
    wendland::wendlandfunction(wdlfunction,l,k);
    wendland::wendlandderivative(wdlfunction,l,k,wdlf1,1);
    wendland::wendlandderivative(wdlf1,l,k,wdlf2,2);
    if(wendland::factorise)
    {
        wendland::syndiv1(l,k,wdlfunction, wdlfunctionfac, &powl0, 0);
        wendland::syndiv1(l,k,wdlf1, wdlf1fac, &powl1, 1);
        wendland::syndiv1(l,k,wdlf2, wdlf2fac, &powl2, 2);
    }else{
        cout << "               PROGRAM STOPPED!!!! " << endl;
        cout << "There has been overflow in the coefficients of the Wendland function." << endl;
        cout << "This only happens if the parameters l or k or both are too big." << endl;
        cout << "The authors of this program recommend using lower parameters to those you have chosen." << endl;
        exit(9);
    }
    
    if (printreport == true) {
        wendland::startreport(1);
        wendland::wendlandreport.close();
    }
}

double WendRBF::operator()(double r) {
    return pow(1.0 - c * r, powl0)*wendland::evawdlfn(r, c, wdlfunctionfac, 0);
}

double WendRBF::aux1(double r) {
    return pow(1.0 - c * r, powl1)*wendland::evawdlfn(r, c, wdlf1fac, 1);
}

double WendRBF::aux2(double r) {
    return pow(1.0 - c * r, powl2)*wendland::evawdlfn(r, c, wdlf2fac, 2);
}


void wendland::wendlandfunction(Mat<double> &wdlf,int l,int k)
{
    //NOTE: comodin is a vector that inside the loop works the divisors of the fractions that the integration produces: Then computes the LCD. After the loop is finished. Then it works the numerators of such fractions and computes the MCD.
    ///wdlf(0,All) IS THE EXPONENTS
    ///wdnf(1,All) IS THE COEFFICIENTS
    int wendlandorder=0;
    
    wall_clock timer;
    startsection("Construction and all steps for the construction of the Wendland function.");
    wendlandreport << "\\textbf{Wendland function }$\\Psi_{"<< l<<","<< k << "}^0$" << "\\\\[1cm]" << endl;
    timer.tic();
    
    int N=l+1;
    int finaldim=N+2*k;
    Mat <double> comodin(1,(int)(finaldim));
    wdlf.resize(2,finaldim);
    wdlf.zeros();
    comodin.zeros();
    comodin.fill(1.0);
    int intcounter=1;
    
    long mcm=0;
    long mcd=0;
    ///FIRST LET US COMPUTER (1-r)^l. NOTICE that it is being computing in the entire vector, so it will include several zeros.
    /// THE HIGHEST EXPONENT STARTS WITH THE VECTOR
    wdlf.resize(2,(int)(N));
    double factormcm=0.0;
    pascal(N,l, wdlf);
    if(printreport){
        wendlandreport << "\\textbf{First binomial} $(1-t_1)^"<< l << "$\\\\" << endl;
        writeequation(0);
        printwendland(1,wdlf,wendlandorder);
        writeequation(1);
    }
    int counter=N;
    for(int j=0; j<k; ++j)
    {
        factormcm=0.0;
        /// PRODUCT BY t
        tprod(wdlf);
        ++counter;
        wdlf.resize(2,(int)(counter));
        fixindex(wdlf);
        if(printreport){
            wendlandreport << "\\textbf{Multiplying by} $t_{"<<intcounter<<"}$:\\\\" << endl;
            writeequation(0);
            printwendland(intcounter,wdlf,wendlandorder);
            writeequation(1);
        }
        ///// Integrating
        tprod(wdlf);
        ++counter;
        wdlf.resize(2,(int)(counter));
        comodin.resize(1,(int)(counter));
        fixindex(wdlf);
        fixindex(comodin);
        wdlf.resize(2,counter);
        wdlf(0,0)=0;
        comodin(0,0)=1;
        
        for(int i=counter-1; i>=1; --i)
        {
            wdlf(1,i)=wdlf(1,i)/wdlf(0,i);
            comodin(0,i)=(i)*comodin(0,i-1);
            factormcm+=wdlf(1,i);
        }
        factormcm+=wdlf(0,0);
        wdlf(1,All)*=-1.0;
        wdlf(1,0)=factormcm;
        mcm=getmcm(comodin);
        comodin(0)=mcm;
        wendlandreport << "\\textbf{Computing integration}: "<< j+1 << "\\\\" << endl;
        ++intcounter;
        
        if(j<k-1)
        {
            if(printreport){
                writeequation(0);
                printwendland(intcounter,wdlf,wendlandorder);
                writeequation(1);
            }
        }else{
            if(printreport){
                wendlandreport << "\\textbf{Integrating from} $r$ \\textbf{to} $1$:\\\\" << endl;
                writeequation(0);
                printwendlandr(false,wdlf,wendlandorder);
                writeequation(1);
            }
        }
    }
    
    //NOW: comodin works the numerators
    comodin.resize((int)(counter));
    comodin.fill(1.0);
    comodin=round(mcm*wdlf(1,All));
    mcd=getmcd(comodin);
    // mcm or mcd are =0 if there was overflow (too large numbers)
    if((mcm>0)&&(mcd>0))
    {
        wendland::factorise=true;
        wdlf(1,All)=round((mcm/mcd)*wdlf(1,All));
        if(printreport){
            wendlandreport << "\\textbf{Using the factor} " << round(mcm/mcd) << " \\textbf{the Wendland function becomes for} $0\\leq cr \\leq 1$:" <<  "\\\\" << endl;
            writeequation(0);
            wendlandreport << "\\Psi_{"<< l<<","<< k << "}^{"<< wendlandorder << "}\\left(cr\\right)=" << endl;
            printwendlandf(wdlf,wendlandorder);
            writeequation(1);
        }
    }else{
        if(printreport){
            wendlandreport << "\\textbf{The Wendland function becomes for} $0\\leq cr \\leq 1$:" <<  "\\\\" << endl;
            writeequation(0);
            wendlandreport << "\\Psi_{"<< l<<","<< k << "}^{"<< wendlandorder << "}\\left(cr\\right)=" << endl;
            printwendlandf(wdlf,wendlandorder);
            writeequation(1);
        }
    }
}

void wendland::wendlandderivative(Mat<double> &wdlfinput, int l, int k,  Mat<double> &wdlf1, int wendlandorder)
{
    startsection("Construction and all steps for the construction of the auxiliar function.");
    
    int dim=(int)wdlfinput.n_cols;
    wdlf1.resize(2,dim);
    wdlf1.zeros();
    for(int i=0; i<dim; ++i)
    {
        wdlf1(0,i)=i-2;
        wdlf1(1,i)=wdlfinput(0,i)*wdlfinput(1,i);
    }
    cleanzeros(wdlf1);
    if(printreport){
        wendlandreport << "\\textbf{Order}: " << wendlandorder << endl;
        wendlandreport << "Wendland function derivative $\\Psi_{"<< l<<","<< k << "}^{"<< wendlandorder << "}$.\\\\ Derive $\\Psi_{"<< l<<","<< k << "}^{"<< wendlandorder-1 << "}\\left(cr\\right)$ by $r$ and divide the result by $r$:" << "\\\\" << endl;
        writeequation(0);
        wendlandreport << "\\Psi_{"<< l<<","<< k << "}^{"<< wendlandorder << "}\\left(cr\\right)=" << endl;
        printwendlandf(wdlf1,wendlandorder);
        writeequation(1);
        wendlandreport << "\\textbf{For} $0< cr \\leq 1$.\\\\[1cm]"<< endl;
        
    }
}


double wendland::evawdlfn(double r, double c, Mat<double> const &wdlfn, int wendlandorder)
{
    double checking=0.0;
    int dim = (int)wdlfn.n_cols;
    double wdlfvalue=0.0;
    double x=c*r;
    checking=1.0-x;
    if(checking>0.0)
    {
        for(int i=dim-1; i>=0; i--)
        {
            wdlfvalue=wdlfvalue*x+wdlfn(1,i);
        }
        wdlfvalue*=pow(x,wdlfn(0,0))*pow(c,2*wendlandorder);
    }else{
        wdlfvalue=0.0;
    }
    return wdlfvalue;
}



void wendland::pascal(int currentdimension, int l, Mat <double> &vector1)
{
    //    int negcounter=0;
    for(int i=0;i<=l;++i)
    {
        vector1(0,i)=i;
        double x=1.0;
        for(int h=0;h<=i;h++)
        {
            //negcounter=i-h;
            vector1(1,h)=pow(-1,h)*x;
            x = x * (i - h) / (h + 1.0);
        }
    }
}

void wendland::tprod(Mat <double> &vector)
{
    int locallength=(int)vector.n_cols;
    vector(0,All)+=ones(locallength).t();
}

void wendland::fixindex(Mat <double> &vector)
{
    int locallength=(int)vector.n_cols;
    int localwidth=(int)vector.n_rows;
    Mat<double> conmutevec;
    if(localwidth>1)
    {
        conmutevec.zeros(localwidth,locallength);
        for(int i=1; i<locallength; ++i)
        {
            conmutevec(0,i)=vector(0,i-1);
            conmutevec(1,i)=vector(1,i-1);
        }
        vector=conmutevec;
    }else{
        conmutevec.zeros(localwidth,locallength);
        for(int i=1; i<locallength; ++i)
        {
            conmutevec(0,i)=vector(0,i-1);
        }
        vector=conmutevec;
    }
}


void wendland::printwendland(int j, Mat<double> &vector1, int wendlandorder)
{
    //vector1=coefficients
    int localdimension=(int)vector1.n_cols;
    int i=0;
    double exponent;
    {
        for(i=0; i<=localdimension-1; ++i)
        {
            exponent=(int)std::floor(std::log10(std::fabs(vector1(1,i))));
            if(abs(vector1(1,i))==1.0)
            {
                if(exponent==0)
                {
                    wendlandreport <<signaturing(i,vector1(1,i))<<"t_{"<<j<<"}^{"<<vector1(0,i)<<"}";
                }else{
                    double mantissa=vector1(1,i)*pow(10,-exponent);
                    wendlandreport <<signaturing(i,vector1(1,i))<<abs(mantissa)<<"\\times10^{"<<exponent<<"}r^{"<<vector1(0,i)<<"}";
                }
            }else{
                if((abs(vector1(1,i))!=0.0)&&(abs(vector1(1,i))!=1.0))
                {
                    exponent=(int)std::floor(std::log10(std::fabs(vector1(1,i))));
                    if(exponent==0)
                    {
                        wendlandreport <<signaturing(i,vector1(1,i))<<abs(vector1(1,i))<<"t_{"<<j<<"}^{"<<vector1(0,i)<<"}";
                    }else{
                        double mantissa=vector1(1,i)*pow(10,-exponent);
                        wendlandreport <<signaturing(i,vector1(1,i))<<abs(mantissa)<<"\\times10^{"<<exponent<<"}t_{"<<j<<"}^{"<<vector1(0,i)<<"}";
                    }
                }else{
                    wendlandreport <<signaturing(i,vector1(1,i))<<abs(vector1(1,i))<<"t_{"<<j<<"}^{"<<vector1(0,i)<<"}";
                }
            }
        }
    }
}

void wendland::printwendlandr(bool finalization, Mat<double> &vector1, int wendlandorder)
{
    //vector1=coefficients
    int localdimension=(int)vector1.n_cols;
    int i=0;
    double exponent;
    if(finalization)
    {
        for(i=0; i<=localdimension-1; ++i)
        {
            exponent=(int)std::floor(std::log10(std::fabs(vector1(1,i))));
            if((abs(vector1(1,i))==1.0)&&(abs(vector1(1,i))!=0.0))
            {
                if(exponent==0)
                {
                    wendlandreport <<signaturing(i,vector1(1,i))<<"c^{"<<vector1(0,i)+2*wendlandorder<<"}r^{"<<vector1(0,i)<<"}";
                }else{
                    double mantissa=vector1(1,i)*pow(10,-exponent);
                    wendlandreport <<signaturing(i,vector1(1,i))<<abs(mantissa)<<"c^{"<<vector1(0,i)+2*wendlandorder<<"}\\times10^{"<<exponent<<"}r^{"<<vector1(0,i)<<"}";
                }
            }else{
                if(abs(vector1(1,i))!=0.0)
                {
                    exponent=(int)std::floor(std::log10(std::fabs(vector1(1,i))));
                    if(exponent==0)
                    {
                        wendlandreport <<signaturing(i,vector1(1,i))<<abs(vector1(1,i))<<"c^{"<<vector1(0,i)+2*wendlandorder<<"}r^{"<<vector1(0,i)<<"}";
                    }else{
                        double mantissa=vector1(1,i)*pow(10,-exponent);
                        wendlandreport <<signaturing(i,vector1(1,i))<<abs(mantissa)<<"\\times10^{"<<exponent<<"}c^{"<<vector1(0,i)+2*wendlandorder<<"}r^{"<<vector1(0,i)<<"}";
                    }
                }else{
                    wendlandreport <<signaturing(i,vector1(1,i))<<abs(vector1(1,i))<<"c^{"<<vector1(0,i)+2*wendlandorder<<"}r^{"<<vector1(0,i)<<"}";
                }
            }
        }
    }else{
        for(i=0; i<=localdimension-1; ++i)
        {
            exponent=(int)std::floor(std::log10(std::fabs(vector1(1,i))));
            if(abs(vector1(1,i))==1.0)
            {
                if(exponent==0)
                {
                    wendlandreport <<signaturing(i,vector1(1,i))<<"r^{"<<vector1(0,i)<<"}";
                }else{
                    double mantissa=vector1(1,i)*pow(10,-exponent);
                    wendlandreport <<signaturing(i,vector1(1,i))<<abs(mantissa)<<"\\times10^{"<<exponent<<"}r^{"<<vector1(0,i)<<"}";
                }
            }else{
                if((abs(vector1(1,i))!=0.0)&&(abs(vector1(1,i))!=1.0))
                {
                    exponent=(int)std::floor(std::log10(std::fabs(vector1(1,i))));
                    if(exponent==0)
                    {
                        wendlandreport <<signaturing(i,vector1(1,i))<<abs(vector1(1,i))<<"r^{"<<vector1(0,i)<<"}";
                    }else{
                        double mantissa=vector1(1,i)*pow(10,-exponent);
                        wendlandreport <<signaturing(i,vector1(1,i))<<abs(mantissa)<<"\\times10^{"<<exponent<<"}r^{"<<vector1(0,i)<<"}";
                    }
                }else{
                    wendlandreport <<signaturing(i,vector1(1,i))<<abs(vector1(1,i))<<"r^{"<<vector1(0,i)<<"}";
                }
            }
        }
    }
}

void wendland::printwendlandf(Mat<double> &vector1, int wendlandorder)
{
    int localdimension=(int)vector1.n_cols;
    int i=0;
    {
        if((abs(vector1(1,0))>=0)&&(abs(vector1(0,0))<=1.0e-10))
        {
            if((vector1(0,0)+2.0*wendlandorder)==0)
            {
                wendlandreport <<(int long)vector1(1,0);
            }else{
                wendlandreport <<(int long)vector1(1,0)<<"c^{"<<vector1(0,0)+2.0*wendlandorder<<"}";
            }
        }else{
            wendlandreport <<(int long)vector1(1,0)<<"c^{"<<vector1(0,0)+2.0*wendlandorder<<"}r^{"<<vector1(0,0)<<"}";
        }
        if(localdimension>1)
        {
            for(i=1; i<localdimension-1; ++i)
            {
                if(abs(vector1(1,i))==1.0)
                {
                    wendlandreport <<signaturing(i,vector1(1,i))<<"c^{"<<vector1(0,i)+2.0*wendlandorder<<"}r^{"<<vector1(0,i)<<"}";
                }else{
                    wendlandreport <<signaturing(i,vector1(1,i))<<(int long)abs(vector1(1,i))<<"c^{"<<vector1(0,i)+2.0*wendlandorder<<"}r^{"<<vector1(0,i)<<"}";
                }
            }
            if(abs(vector1(1,localdimension-1))==1.0)
            {
                wendlandreport <<signaturing(localdimension-1,vector1(1,localdimension-1))<<"c^{"<<vector1(0,localdimension-1)+2.0*wendlandorder<<"}r^{"<<vector1(0,localdimension-1)<<"}"<<endl;
            }else{
                wendlandreport <<signaturing(localdimension-1,vector1(1,localdimension-1))<<(int long)abs(vector1(1,localdimension-1))<<"c^{"<<vector1(0,localdimension-1)+2.0*wendlandorder<<"}r^{"<<vector1(0,localdimension-1)<<"}"<<endl;
            }
        }
    }
}


void wendland::printwendlandff(Mat<double> &vector1, int wendlandorder)
{
    int localdimension=(int)vector1.n_cols;
    int i=0;
    if(wendlandorder==0)
    {
        if((abs(vector1(1,0))>=0)&&(abs(vector1(0,0))<=1.0e-10))
        {
            if((vector1(0,0)+2.0*wendlandorder)==0)
            {
                wendlandreport <<(int long)vector1(1,0);
            }else{
                wendlandreport <<(int long)vector1(1,0)<<"c^{"<<vector1(0,0)+2.0*wendlandorder<<"}";
            }
        }else{
            wendlandreport <<(int long)vector1(1,0)<<"c^{"<<vector1(0,0)+2.0*wendlandorder<<"}r^{"<<vector1(0,0)<<"}";
        }
        if(localdimension>1)
        {
            for(i=1; i<localdimension-1; ++i)
            {
                if(abs(vector1(1,i))==1.0)
                {
                    wendlandreport <<signaturing(i,vector1(1,i))<<"c^{"<<vector1(0,i)+2.0*wendlandorder<<"}r^{"<<vector1(0,i)<<"}";
                }else{
                    wendlandreport <<signaturing(i,vector1(1,i))<<(int long)abs(vector1(1,i))<<"c^{"<<vector1(0,i)+2.0*wendlandorder<<"}r^{"<<vector1(0,i)<<"}";
                }
            }
            if(abs(vector1(1,localdimension-1))==1.0)
            {
                wendlandreport <<signaturing(localdimension-1,vector1(1,localdimension-1))<<"c^{"<<vector1(0,localdimension-1)+2.0*wendlandorder<<"}r^{"<<vector1(0,localdimension-1)<<"}"<<endl;
            }else{
                wendlandreport <<signaturing(localdimension-1,vector1(1,localdimension-1))<<(int long)abs(vector1(1,localdimension-1))<<"c^{"<<vector1(0,localdimension-1)+2.0*wendlandorder<<"}r^{"<<vector1(0,localdimension-1)<<"}"<<endl;
            }
        }
    }else{
        if(localdimension>1)
        {
            if((abs(vector1(1,0))>=0)&&(abs(vector1(0,0))<=1.0e-10))
            {
                if(abs(vector1(0,0))<=1.0e-10)
                {
                    wendlandreport <<(int long)vector1(1,0);
                }else{
                    wendlandreport <<(int long)vector1(1,0)<<"c^{"<<vector1(0,0)<<"}";
                }
            }else{
                wendlandreport <<(int long)vector1(1,0)<<"c^{"<<vector1(0,0)<<"}r^{"<<vector1(0,0)<<"}";
            }
            for(i=1; i<localdimension-1; ++i)
            {
                if(abs(vector1(1,i))==1.0)
                {
                    wendlandreport <<signaturing(i,vector1(1,i))<<"c^{"<<vector1(0,i)<<"}r^{"<<vector1(0,i)<<"}";
                }else{
                    wendlandreport <<signaturing(i,vector1(1,i))<<(int long)abs(vector1(1,i))<<"c^{"<<vector1(0,i)<<"}r^{"<<vector1(0,i)<<"}";
                }
            }
            if(abs(vector1(1,localdimension-1))==1.0)
            {
                wendlandreport <<signaturing(localdimension-1,vector1(1,localdimension-1))<<"c^{"<<vector1(0,localdimension-1)<<"}r^{"<<vector1(0,localdimension-1)<<"}"<<endl;
            }else{
                wendlandreport <<signaturing(localdimension-1,vector1(1,localdimension-1))<<(int long)abs(vector1(1,localdimension-1))<<"c^{"<<vector1(0,localdimension-1)<<"}r^{"<<vector1(0,localdimension-1)<<"}"<<endl;
            }
        }else{
            if((abs(vector1(1,0))>=0)&&(abs(vector1(0,0))<=1.0e-10))
            {
                    wendlandreport <<(int long)vector1(1,0);
                
            }else{
                wendlandreport <<(int long)vector1(1,0)<<"r^{"<<vector1(0,0)<<"}";
            }
        }
    }
}

char wendland::signaturing(int i, double unknownsign)
{
    char algebraicsymbol;
    if((i==0)&&(arma::sign(unknownsign)>=0))
    {
        algebraicsymbol='\~';
    }else{
        if(arma::sign(unknownsign)>=0)
        {
            algebraicsymbol='+';
        }else{
            algebraicsymbol='-';
        }
    }
    return algebraicsymbol;
}

void wendland::cleanzeros(Mat<double> &vector)
{
    
    int localdimension=(int)vector.n_cols;
    Mat<double> temp(2,localdimension);
    temp.zeros();
    int count=0;
    for(int i=0; i<localdimension; ++i)
    {
        if(vector(1,i)==0)
        {
            ++count;
        }else{
            break;
        }
    }
    if(count<localdimension)
    {
        temp.resize(2,localdimension-count);
        temp(All,All)=vector(All,span(count,localdimension-1));
        vector.resize(2,localdimension-count);
        vector=temp;
    }
}

void wendland::startreport(int wendlandorder)
{
    if(printreport)
    {
        if(wendlandorder==0)
        {
            wendlandreport <<"\\documentclass[a4paper,twoside]{report}" << endl;
            wendlandreport <<"\\usepackage{amssymb}" << endl;
            wendlandreport <<"\\usepackage{amstext}" << endl;
            wendlandreport <<"\\usepackage{amsmath}" << endl;
            wendlandreport <<"\\usepackage{amsthm}" << endl;
            wendlandreport <<"\\usepackage[colorlinks=true, linkcolor=black, urlcolor=blue, citecolor=blue, anchorcolor=blue]{hyperref}" << endl;
            wendlandreport <<"\\renewcommand{\\thesection}{\\arabic{section}}" << endl;
            wendlandreport <<"\\setcounter{chapter}{1}" << endl;
            wendlandreport <<"\\begin{document}" << endl;
            wendlandreport <<"\\begin{center}" << endl;
            wendlandreport <<"{\\bf{\\Large Construction of the Wendland function with \\\\ WendlandXool}} \\\\" << endl;
            wendlandreport <<"\\vspace*{.2in}" << endl;
            wendlandreport <<"\\begin{tabular}{cc}" << endl;
            wendlandreport <<"  Carlos Arg\\'aez\\textsuperscript{1,*}, Peter Giesl\\textsuperscript{2}, Sigurdur Freyr Hafstein\\textsuperscript{3}" << endl;
            wendlandreport <<"\\\\[0.25ex]" << endl;
            wendlandreport <<"{\\small \\textsuperscript{1} Science Institute, University of Iceland, carlos@hi.is} \\\\" << endl;
            wendlandreport <<"{\\small \\textsuperscript{2} Department of Mathematics, University of Sussex, P.A.Giesl@sussex.ac.uk} \\\\" << endl;
            wendlandreport <<"{\\small \\textsuperscript{3} Science Institute, University of Iceland, shafstein@hi.is} \\\\" << endl;
            wendlandreport <<"\\end{tabular}" << endl;
            wendlandreport <<"\\end{center}" << endl;
            wendlandreport <<"\\let\\thefootnote\\relax\\footnotetext{{\\small * corresponding author}}" << endl;
        }else{
            wendlandreport <<"\\begin{thebibliography}{9}" << endl;
            wendlandreport <<"\\bibitem{latexcompanion}" << endl;
            wendlandreport <<"Wendland, Holger (2005). \\textit{Scattered Data Approximation}. Cambridge: Cambridge University Press. pp. 11, 18-23, 64-66. ISBN 0521843359." << endl;
            wendlandreport <<"\\bibitem{wendlandxool1}" << endl;
            wendlandreport <<"Arg\\'aez, C., and Giesl, P., Hafstein, S.F. (2017). \\textit{Wendland Functions - A C++ Code to Compute Them.}. In Proceedings of the 7th International Conference on Simulation and Modeling Methodologies, Technologies and Applications (SIMULTECH 2017), pp. 323-330 ISBN: 978-989-758-265-3." << endl;
            wendlandreport <<"\\bibitem{wendlandxool2}" << endl;
            wendlandreport <<"Arg\\'aez, C., and Giesl, P., Hafstein, S.F. (2017). \\textit{WendlandXool: Simplified C++ code to compute Wendland functions}. Submitted." << endl;
            wendlandreport <<"\\end{thebibliography}" << endl;
            wendlandreport <<"\\end{document}" << endl;
        }
    }
}

void wendland::writeequation(const int classifier)
{
    if(classifier==0)
    {
        wendlandreport << "\\begin{flushleft}" << endl;
        wendlandreport << "$" << endl;
    }
    if(classifier==1)
    {
        wendlandreport << "$" << endl;
        wendlandreport << "\\end{flushleft}" << endl;
    }
}

void wendland::mcmp(long long unsigned value1, long long unsigned value2, long long unsigned &valueout)
{
    while (value2 > 0) {
        long long unsigned r = value1 % value2;
        value1 = value2;
        value2 = r;
    }
    valueout=value1;
}

void wendland::gcmp(long long unsigned value1, long long unsigned value2, long long unsigned &valueout)
{
    while (value2 > 0) {
        long long unsigned r = value1 % value2;
        value1 = value2;
        value2 = r;
    }
    valueout=value1;
}

long long unsigned wendland::getmcm(Mat<double> &vector)
{
    long long unsigned valueout=1;
    long long unsigned mcm=vector(0,0);
    vector=abs(vector);
    for(int i=0; i<=(int)vector.n_cols-1; ++i)
    {
        if((mcm == 0) || (vector(0,i) == 0)){
            break;
        }else{
            mcmp(mcm, vector(0,i),valueout);
            mcm=(mcm * vector(0,i)) / valueout; //
        }
    }
    return mcm;
}

long long unsigned wendland::getmcd(Mat<double> &vector)
{
    long long unsigned valueout=1;
    long long unsigned mcd=vector(0,0);
    vector=abs(vector);
    for(int i=1; i<(int)vector.n_cols;++i)
    {
        if(vector(0,i)!=0.0)
        {
            gcmp(mcd, (long long unsigned)abs(vector(0,i)),valueout);
            mcd=valueout;
        }
    }
    return mcd;
}

bool wendland::syndiv1(int l, int k, const arma::Mat<double> &wdlin, arma::Mat<double> &wdlout, int *pow1mr,int order) {
    *pow1mr = 0;
    if(printreport){
        if(order==0)
        {
            wendlandreport << "\\chapter\*\{Appendix\}" << endl;
            startsection("Functions presented in a factorised form.");
            wendlandreport << "Next, we present the factorised version of the Wendland function.\\\\" << endl;
        }
    }
    double tol = 1e-8;
    int m = wdlin.n_cols;
    vector<long long> intcoeff(m);
    for (int i = 0; i < m; i++) {
        intcoeff[m-i-1] = round(wdlin(1, i));
        if (abs(intcoeff[m - i - 1] - wdlin(1, i)) > tol) {
            wdlout = wdlin;
            return false;  // not integer coefficients
        }
    }
    
    while (true) {
        if (m - *pow1mr <= 1) {
            break;
        }
        for (int i = 1; i < m - *pow1mr; i++) {
            intcoeff[i] += intcoeff[i - 1];
        }
        if (intcoeff[m - *pow1mr - 1] == 0) {
            (*pow1mr)++;
        }
        else {
            break;
        }
    }
    // code intcoeff back
    for (int i = m - (*pow1mr)-1; i >= 1; i--) {
        intcoeff[i] -= intcoeff[i - 1];
    }
    
    int sign = ((*pow1mr) % 2 == 1 ? -1 : 1);
    wdlout.set_size(2, m - (*pow1mr));
    for (int i = 0; i < m - (*pow1mr); i++) {
        wdlout(0, i) = wdlin(0, i);
        wdlout(1, i) = double(sign*intcoeff[m-(*pow1mr)-1-i]);
    }
    if(printreport){
        writeequation(0);
        wendlandreport << "\\Psi_{"<< l <<","<< k << "}^{"<< order << "}\\left(cr\\right)=" << endl;
        if(order==0)
        {
            wendlandreport << "\\left(1-cr\\right)_+^{"<<*pow1mr<<"}\\left(" << endl;
        }else{
            wendlandreport << "\\left(1-cr\\right)_+^{"<<*pow1mr<<"}c^{"<<wdlout(0,0)+2.0*order<<"}\\left(" << endl;
        }
        printwendlandff(wdlout,order);
        wendlandreport << "\\right)" << endl;
        writeequation(1);
    }
    return true;
}

void wendland::startsection(string valor)
{
    wendlandreport  << "\\section{"<< valor << "}" << endl;
}

