/* DNA2AA.cpp

Description:
This script converts DNA sequences to AA sequences.

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

  1) Place this file in the same folder as "dnaseqs" created with ExtractSeqs.cpp
  2) Compile it and run: $ g++ DNA2AA.cpp -o DNA2AA
                         $ ./DNA2AA
  3) Amino acid sequences will be stored in "AASEQS".
*/ 

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <set>

using namespace std ;

int main()
{
    ifstream file ;
    ifstream dna ;
    ofstream matrix ;
    vector <string> indna ;
    vector <string> dnalist ;
    vector <double> SeqVec ;
    vector <double> distances ;
    string v1, v2 ;

    string dnastring = "" ;
    vector <string>  POPN ; //Population in vector form

    double NewValue =  0 ;
    double D2 = 0 ; //delta^2
    string line ;

    string originaldna ;  
    dna.open("dna") ;
    dna >> originaldna ;
    int Dimension = originaldna.length() / 3 ;
    
    file.open("dnaseqs") ;
    
    matrix.open("AASEQS") ;
    
    //Gets sequences from file
    while( file >> line )
    {

        if ( line.length() > (originaldna.length()-1) )
        {
            //cout << line << endl ;
            
            indna.push_back(line) ;
            
            line = "" ;

	}
    }

	set<string> s;
	unsigned size = indna.size();
	for( unsigned i = 0; i < size; ++i ) 
	{ s.insert( indna[i] ) ; }
	dnalist.assign( s.begin(), s.end() );

cout << "DNA LIST SIZE:" << dnalist.size() << endl;

             //cout << "hello 1" << endl ;
for( int N = 0 ; N < dnalist.size() ; N++ )
{
	for( int L = 0 ; L < originaldna.length() ; L += 3 )
	{   
                line = "" ;
                //cout << "hello 2" << endl ;
                line = dnalist[N] ; 
		//cout << line << endl ;

        if      (line.substr( L , 3 ) == "TCA")    { dnastring += "s";}
    	else if (line.substr( L , 3 ) == "TCC")    { dnastring += "s";}   
        else if (line.substr( L , 3 ) == "TCG")    { dnastring += "s";}
    	else if (line.substr( L , 3 ) == "TCT")    { dnastring += "s";}    
    	else if (line.substr( L , 3 ) == "TTC")    { dnastring += "f";}
    	else if (line.substr( L , 3 ) == "TTT")    { dnastring += "f";}    
    	else if (line.substr( L , 3 ) == "TTA")    { dnastring += "l";} 
    	else if (line.substr( L , 3 ) == "TTG")    { dnastring += "l";}   
    	else if (line.substr( L , 3 ) == "TAC")    { dnastring += "y";} 
    	else if (line.substr( L , 3 ) == "TAT")    { dnastring += "y";}  
    	else if (line.substr( L , 3 ) == "TAA")    { dnastring += "_";} 
    	else if (line.substr( L , 3 ) == "TAG")    { dnastring += "_";} 
    	else if (line.substr( L , 3 ) == "TGC")    { dnastring += "c";}  
    	else if (line.substr( L , 3 ) == "TGT")    { dnastring += "c";}  
    	else if (line.substr( L , 3 ) == "TGA")    { dnastring += "_";}  
    	else if (line.substr( L , 3 ) == "TGG")    { dnastring += "w";} 
    	else if (line.substr( L , 3 ) == "CTA")    { dnastring += "l";}  
    	else if (line.substr( L , 3 ) == "CTC")    { dnastring += "l";}  
    	else if (line.substr( L , 3 ) == "CTG")    { dnastring += "l";}  
    	else if (line.substr( L , 3 ) == "CTT")    { dnastring += "l";}  
    	else if (line.substr( L , 3 ) == "CCA")    { dnastring += "p";}   
    	else if (line.substr( L , 3 ) == "CCC")    { dnastring += "p";}   
    	else if (line.substr( L , 3 ) == "CCG")    { dnastring += "p";}   
    	else if (line.substr( L , 3 ) == "CCT")    { dnastring += "p";}  
    	else if (line.substr( L , 3 ) == "CAC")    { dnastring += "h";}  
    	else if (line.substr( L , 3 ) == "CAT")    { dnastring += "h";}  
    	else if (line.substr( L , 3 ) == "CAA")    { dnastring += "q";}   
    	else if (line.substr( L , 3 ) == "CAG")    { dnastring += "q";}   
    	else if (line.substr( L , 3 ) == "CGA")    { dnastring += "r";}   
    	else if (line.substr( L , 3 ) == "CGC")    { dnastring += "r";}  
    	else if (line.substr( L , 3 ) == "CGG")    { dnastring += "r";}   
    	else if (line.substr( L , 3 ) == "CGT")    { dnastring += "r";}   
    	else if (line.substr( L , 3 ) == "ATA")    { dnastring += "i";}  
    	else if (line.substr( L , 3 ) == "ATC")    { dnastring += "i";}   
    	else if (line.substr( L , 3 ) == "ATT")    { dnastring += "i";}   
    	else if (line.substr( L , 3 ) == "ATG")    { dnastring += "m";}   
    	else if (line.substr( L , 3 ) == "ACA")    { dnastring += "t";}    
    	else if (line.substr( L , 3 ) == "ACC")    { dnastring += "t";}    
    	else if (line.substr( L , 3 ) == "ACG")    { dnastring += "t";}   
    	else if (line.substr( L , 3 ) == "ACT")    { dnastring += "t";}   
    	else if (line.substr( L , 3 ) == "AAC")    { dnastring += "n";}    
    	else if (line.substr( L , 3 ) == "AAT")    { dnastring += "n";}    
    	else if (line.substr( L , 3 ) == "AAA")    { dnastring += "k";}   
    	else if (line.substr( L , 3 ) == "AAG")    { dnastring += "k";}   
    	else if (line.substr( L , 3 ) == "AGC")    { dnastring += "s";}   
    	else if (line.substr( L , 3 ) == "AGT")    { dnastring += "s";}   
    	else if (line.substr( L , 3 ) == "AGA")    { dnastring += "r";}   
    	else if (line.substr( L , 3 ) == "AGG")    { dnastring += "r";}  
    	else if (line.substr( L , 3 ) == "GTA")    { dnastring += "v";}   
    	else if (line.substr( L , 3 ) == "GTC")    { dnastring += "v";}   
    	else if (line.substr( L , 3 ) == "GTG")    { dnastring += "v";}    
    	else if (line.substr( L , 3 ) == "GTT")    { dnastring += "v";}   
    	else if (line.substr( L , 3 ) == "GCA")    { dnastring += "a";}    
    	else if (line.substr( L , 3 ) == "GCC")    { dnastring += "a";}   
    	else if (line.substr( L , 3 ) == "GCG")    { dnastring += "a";}    
    	else if (line.substr( L , 3 ) == "GCT")    { dnastring += "a";}    
    	else if (line.substr( L , 3 ) == "GAC")    { dnastring += "d";} 
    	else if (line.substr( L , 3 ) == "GAT")    { dnastring += "d";}  
    	else if (line.substr( L , 3 ) == "GAA")    { dnastring += "e";}  
    	else if (line.substr( L , 3 ) == "GAG")    { dnastring += "e";}    
    	else if (line.substr( L , 3 ) == "GGA")    { dnastring += "g";}   
    	else if (line.substr( L , 3 ) == "GGC")    { dnastring += "g";}    
    	else if (line.substr( L , 3 ) == "GGG")    { dnastring += "g";}    
    	else if (line.substr( L , 3 ) == "GGT")    { dnastring += "g";}

	}

	//cout << dnastring << endl ;
	matrix<< dnastring<<endl ;

	dnastring = "" ;

}

	file.close() ;
	dna.close() ;
	matrix.close() ;

	return 0 ;
    
}
