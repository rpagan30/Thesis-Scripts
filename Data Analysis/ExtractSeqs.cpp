/* ExtractSeqs.cpp

Description:
This program extracts DNA sequences from 'out' files.

Created by Rafael Pagan, rpagan30@gmail.com.

MIT License

Copyright(c) 2017 Rafael Pagan, Steven Massey, & Julian Velev.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

  USAGE:
  1) Place script in the same folder as the seed sequence and the out file produced by
     nuneutralrobustness.cpp.
  2) Type the following into the terminal: 
     $ g++ ExtractSeqs.cpp -o Extract
     $ ./Extract
  3) The extracted sequences will be found in "dnaseqs".
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

using namespace std ;


int main()
{
    remove("dnaseqs") ;

    string line = "" ;
    string originaldna ;

    ifstream file ;
    ifstream dna ;
    ofstream dnaseqs;
    
    //Gets original DNA size.
    dna.open("dna") ;
    
    dna >> originaldna ;
    
    int seq_Length = originaldna.length() ;

    file.open("out") ;

    dnaseqs.open("dnaseqs", ios::app ) ;

    while( file >> line )
    {
        if ( line.length() > seq_Length - 1 )
        {
	    dnaseqs << line << endl ;
        
            line = "" ;
        }
    }
    
    file.close() ;

    dnaseqs.close() ;

    return 0 ;
    
}
