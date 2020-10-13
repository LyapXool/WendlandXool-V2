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
#include <fstream>
#include <armadillo>
#include <list>
#include <filesystem>
#include <cstdio>
#include <cmath>
#include "wendland.hpp"


#include<armadillo>

using namespace std;




int main()
{
    system("pwd");
	double r = 0.3, c = 0.5;
	WendRBF wend(3,1,c,true);
	cout << wend(r) << " " << wend.aux1(r) << " " << wend.aux2(r) << endl;
    return 0;
}




