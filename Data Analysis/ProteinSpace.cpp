/* ProteinSpace.cpp

Description:
This script calculates the movement of populations over protein space. It requires
an input file with all of the DNA sequences of the simulation in order of their 
generation and produces a vector of amino acid distances. Each value represents the 
distance of each population from the original seed sequence.

It assigns a value from 0 to 20 to each codon according to its amino acid's
hydrophobicity using the Kyte-Doolittle scale. For N sequences in each population, the 
average hydrophobicity is calculated for index L. This is it's center of mass in that 
dimension. 

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
1) Place script in the same folder as the seed sequence and the dna file created with 
   ExtractSeqs.cpp
2) Type the following into the terminal: 
   $ g++ ExtractSeqs.cpp -o Extract
   $ ./Extract
3) The vector of calculated protein space distances will be output to "ProteinSpace.txt".
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

using namespace std ;

int main()
{
    remove("ProteinSpace.txt") ;
    
    ifstream file ;
    ifstream dna ;
    ofstream ProteinSpace ;
    vector<string> dnalist ;
    vector<double> coordinatevalue ;
    vector<double> CoM ; //Center of Mass vector for L/3 dimensions where L is the number of bases and L/3 is the number of amino acids.
    double NewValue =  0 ;
    double CoMSum = 0 ;
    
    double D2 = 0 ; //delta^2
    string line ;
    string originaldna ;
    
    //Converts original sequence to a vector in L/3-D space.
    dna.open("dna") ;
    
    dna >> originaldna ;
    
    vector<double> R0 ;// Vector that stores seed sequence as a vector in protein space.
    
    int Dimension = originaldna.length() / 3 ;
    
    //cout << "Dimension = " << Dimension << endl ;
    
    for( int L = 0 ; L < originaldna.length() ; L += 3 )
    {
        line = originaldna ; 
        
        if      (line.substr( L , 3 ) == "TCA")    { R0.push_back(10);   }
        else if (line.substr( L , 3 ) == "TCC")    { R0.push_back(10);   }
        else if (line.substr( L , 3 ) == "TCG")    { R0.push_back(10);   }
        else if (line.substr( L , 3 ) == "TCT")    { R0.push_back(10);   }
        else if (line.substr( L , 3 ) == "TTC")    { R0.push_back(16);   }
        else if (line.substr( L , 3 ) == "TTT")    { R0.push_back(16);   }
        else if (line.substr( L , 3 ) == "TTA")    { R0.push_back(17);   }
        else if (line.substr( L , 3 ) == "TTG")    { R0.push_back(17);   }
        else if (line.substr( L , 3 ) == "TAC")    { R0.push_back(8);    }
        else if (line.substr( L , 3 ) == "TAT")    { R0.push_back(8);    }
        else if (line.substr( L , 3 ) == "TGC")    { R0.push_back(15);   }
        else if (line.substr( L , 3 ) == "TGT")    { R0.push_back(15);   }
        else if (line.substr( L , 3 ) == "TGG")    { R0.push_back(9);    }
        else if (line.substr( L , 3 ) == "CTA")    { R0.push_back(17);   }
        else if (line.substr( L , 3 ) == "CTC")    { R0.push_back(17);   }
        else if (line.substr( L , 3 ) == "CTG")    { R0.push_back(17);   }
        else if (line.substr( L , 3 ) == "CTT")    { R0.push_back(17);   }
        else if (line.substr( L , 3 ) == "CCA")    { R0.push_back(7);    }
        else if (line.substr( L , 3 ) == "CCC")    { R0.push_back(7);    }
        else if (line.substr( L , 3 ) == "CCG")    { R0.push_back(7);    }
        else if (line.substr( L , 3 ) == "CCT")    { R0.push_back(7);    }
        else if (line.substr( L , 3 ) == "CAC")    { R0.push_back(6);    }
        else if (line.substr( L , 3 ) == "CAT")    { R0.push_back(6);    }
        else if (line.substr( L , 3 ) == "CAA")    { R0.push_back(4);    }
        else if (line.substr( L , 3 ) == "CAG")    { R0.push_back(4);    }
        else if (line.substr( L , 3 ) == "CGA")    { R0.push_back(0);    }
        else if (line.substr( L , 3 ) == "CGC")    { R0.push_back(0);    }
        else if (line.substr( L , 3 ) == "CGG")    { R0.push_back(0);    }
        else if (line.substr( L , 3 ) == "CGT")    { R0.push_back(0);    }
        else if (line.substr( L , 3 ) == "ATA")    { R0.push_back(19);   }
        else if (line.substr( L , 3 ) == "ATC")    { R0.push_back(19);   }
        else if (line.substr( L , 3 ) == "ATT")    { R0.push_back(19);   }
        else if (line.substr( L , 3 ) == "ATG")    { R0.push_back(14);   }
        else if (line.substr( L , 3 ) == "ACA")    { R0.push_back(11);   }
        else if (line.substr( L , 3 ) == "ACC")    { R0.push_back(11);   }
        else if (line.substr( L , 3 ) == "ACG")    { R0.push_back(11);   }
        else if (line.substr( L , 3 ) == "ACT")    { R0.push_back(11);   }
        else if (line.substr( L , 3 ) == "AAC")    { R0.push_back(2);    }
        else if (line.substr( L , 3 ) == "AAT")    { R0.push_back(2);    }
        else if (line.substr( L , 3 ) == "AAA")    { R0.push_back(1);    }
        else if (line.substr( L , 3 ) == "AAG")    { R0.push_back(1);    }
        else if (line.substr( L , 3 ) == "AGC")    { R0.push_back(10);   }
        else if (line.substr( L , 3 ) == "AGT")    { R0.push_back(10);   }
        else if (line.substr( L , 3 ) == "AGA")    { R0.push_back(0);    }
        else if (line.substr( L , 3 ) == "AGG")    { R0.push_back(0);    }
        else if (line.substr( L , 3 ) == "GTA")    { R0.push_back(18);   }
        else if (line.substr( L , 3 ) == "GTC")    { R0.push_back(18);   }
        else if (line.substr( L , 3 ) == "GTG")    { R0.push_back(18);   }
        else if (line.substr( L , 3 ) == "GTT")    { R0.push_back(18);   }
        else if (line.substr( L , 3 ) == "GCA")    { R0.push_back(13);   }
        else if (line.substr( L , 3 ) == "GCC")    { R0.push_back(13);   }
        else if (line.substr( L , 3 ) == "GCG")    { R0.push_back(13);   }
        else if (line.substr( L , 3 ) == "GCT")    { R0.push_back(13);   }
        else if (line.substr( L , 3 ) == "GAC")    { R0.push_back(3);    }
        else if (line.substr( L , 3 ) == "GAT")    { R0.push_back(3);    }
        else if (line.substr( L , 3 ) == "GAA")    { R0.push_back(5);    }
        else if (line.substr( L , 3 ) == "GAG")    { R0.push_back(5);    }
        else if (line.substr( L , 3 ) == "GGA")    { R0.push_back(12);   }
        else if (line.substr( L , 3 ) == "GGC")    { R0.push_back(12);   }
        else if (line.substr( L , 3 ) == "GGG")    { R0.push_back(12);   }
        else if (line.substr( L , 3 ) == "GGT")    { R0.push_back(12);   }
    }
    
    cout << "Original Sequence Vector" << endl ;
    
    cout << "(" ;
    for ( int i = 0 ; i < Dimension ; i++ )
    {
        cout << R0[i]  ;
        
        if ( i < Dimension - 1 )
        {
            cout << ", " ;
        }
    }
    cout << ")" << endl << endl ;
    
    file.open("out") ;
    
    ProteinSpace.open("ProteinSpace.txt") ;
    
    //Gets sequences from file
    while( file >> line )
    {
        if ( line.length() > originaldna.length()-2 )
        {
            //cout << line << endl ;
            
            dnalist.push_back(line) ;
            
            line = "" ;
        }
    }
    cout << " DNA List Size = " << dnalist.size() << endl << endl ;

    // Loops over all codons at position L per generation. Assigns a value from 0 to 20 to each codon according to its amino acid's hydrophobicity using the Kyte-Doolittle scale. The average of all codon values over N in position L yields the center of mass (or centroid) in that dimension/codon position. For each generation this yields a center of mass vector over all codon positions (in L/3-Dimension space).
    
    for( int A = 0 ; A < 2000000 ; A += 1000 )
    {
        
        int NF = A + 1000 - 1 ;
        
        //cout << "NF = " << NF << endl ;
        
        //cout << "dnalist[A] = " << dnalist [A] << endl ;
        
        for( int L = 0 ; L < originaldna.length() ; L += 3 )
        {
            
            CoMSum = 0 ;
            
	        coordinatevalue.clear() ;
            
            for( int N = A ; N < NF ; N++ )
            {
                line = "" ;
                
                line = dnalist[N] ;
                //cout << "N = " << N << ", dnalist[N] = " << line << endl ;
                //cout << "Substring: " << line.substr( L , 3 ) << endl ;
                
                if      (line.substr( L , 3 ) == "TCA")    { NewValue =  10 ; }
                else if (line.substr( L , 3 ) == "TCC")    { NewValue =  10 ; }
                else if (line.substr( L , 3 ) == "TCG")    { NewValue = 10 ; }
                else if (line.substr( L , 3 ) == "TCT")    { NewValue = 10 ; }
                else if (line.substr( L , 3 ) == "TTC")    { NewValue = 16 ; }
                else if (line.substr( L , 3 ) == "TTT")    { NewValue = 16 ; }
                else if (line.substr( L , 3 ) == "TTA")    { NewValue = 17 ; }
                else if (line.substr( L , 3 ) == "TTG")    { NewValue = 17 ; }
                else if (line.substr( L , 3 ) == "TAC")    { NewValue = 8 ;  }
                else if (line.substr( L , 3 ) == "TAT")    { NewValue = 8 ;  }
                else if (line.substr( L , 3 ) == "TGC")    { NewValue = 15 ; }
                else if (line.substr( L , 3 ) == "TGT")    { NewValue = 15 ; }
                else if (line.substr( L , 3 ) == "TGG")    { NewValue = 9 ;  }
                else if (line.substr( L , 3 ) == "CTA")    { NewValue = 17 ; }
                else if (line.substr( L , 3 ) == "CTC")    { NewValue = 17; }
                else if (line.substr( L , 3 ) == "CTG")    { NewValue = 17; }
                else if (line.substr( L , 3 ) == "CTT")    { NewValue = 17; }
                else if (line.substr( L , 3 ) == "CCA")    { NewValue = 7;  }
                else if (line.substr( L , 3 ) == "CCC")    { NewValue = 7;  }
                else if (line.substr( L , 3 ) == "CCG")    { NewValue = 7;  }
                else if (line.substr( L , 3 ) == "CCT")    { NewValue = 7;  }
                else if (line.substr( L , 3 ) == "CAC")    { NewValue = 6;  }
                else if (line.substr( L , 3 ) == "CAT")    { NewValue = 6;  }
                else if (line.substr( L , 3 ) == "CAA")    { NewValue = 4;  }
                else if (line.substr( L , 3 ) == "CAG")    { NewValue = 4;  }
                else if (line.substr( L , 3 ) == "CGA")    { NewValue = 0;  }
                else if (line.substr( L , 3 ) == "CGC")    { NewValue = 0;  }
                else if (line.substr( L , 3 ) == "CGG")    { NewValue = 0;  }
                else if (line.substr( L , 3 ) == "CGT")    { NewValue = 0;  }
                else if (line.substr( L , 3 ) == "ATA")    { NewValue = 19; }
                else if (line.substr( L , 3 ) == "ATC")    { NewValue = 19; }
                else if (line.substr( L , 3 ) == "ATT")    { NewValue = 19; }
                else if (line.substr( L , 3 ) == "ATG")    { NewValue = 14; }
                else if (line.substr( L , 3 ) == "ACA")    { NewValue = 11; }
                else if (line.substr( L , 3 ) == "ACC")    { NewValue = 11; }
                else if (line.substr( L , 3 ) == "ACG")    { NewValue = 11; }
                else if (line.substr( L , 3 ) == "ACT")    { NewValue = 11; }
                else if (line.substr( L , 3 ) == "AAC")    { NewValue = 2;  }
                else if (line.substr( L , 3 ) == "AAT")    { NewValue = 2;  }
                else if (line.substr( L , 3 ) == "AAA")    { NewValue = 1;  }
                else if (line.substr( L , 3 ) == "AAG")    { NewValue = 1;  }
                else if (line.substr( L , 3 ) == "AGC")    { NewValue = 10; }
                else if (line.substr( L , 3 ) == "AGT")    { NewValue = 10; }
                else if (line.substr( L , 3 ) == "AGA")    { NewValue = 0;  }
                else if (line.substr( L , 3 ) == "AGG")    { NewValue = 0;  }
                else if (line.substr( L , 3 ) == "GTA")    { NewValue = 18; }
                else if (line.substr( L , 3 ) == "GTC")    { NewValue = 18; }
                else if (line.substr( L , 3 ) == "GTG")    { NewValue = 18; }
                else if (line.substr( L , 3 ) == "GTT")    { NewValue = 18; }
                else if (line.substr( L , 3 ) == "GCA")    { NewValue = 13; }
                else if (line.substr( L , 3 ) == "GCC")    { NewValue = 13; }
                else if (line.substr( L , 3 ) == "GCG")    { NewValue = 13; }
                else if (line.substr( L , 3 ) == "GCT")    { NewValue = 13; }
                else if (line.substr( L , 3 ) == "GAC")    { NewValue = 3;  }
                else if (line.substr( L , 3 ) == "GAT")    { NewValue = 3;  }
                else if (line.substr( L , 3 ) == "GAA")    { NewValue = 5;  }
                else if (line.substr( L , 3 ) == "GAG")    { NewValue = 5;  }
                else if (line.substr( L , 3 ) == "GGA")    { NewValue = 12; }
                else if (line.substr( L , 3 ) == "GGC")    { NewValue = 12; }
                else if (line.substr( L , 3 ) == "GGG")    { NewValue = 12; }
                else if (line.substr( L , 3 ) == "GGT")    { NewValue = 12; }
                //else                                       { cout<<"Unrecognized codon: "<<codon<<endl; Unrecognized++; }
                //cout <<  "L: " << L << " " << "N: " << N << endl ;
                //cout << "New Value = " << NewValue << endl ;
                CoMSum += NewValue ;
                
                NewValue = 0 ;
                //cout << " L = " << L << " N = " << N << " CoMSum =" << CoMSum <<  endl ;
            }
            /**
             cout << "sum of x = " << CoMSum << endl << endl ;
             
             cout << "Number of elements = " << 1000 << endl << endl ;
             
             cout << "Sum of x / num of elements = " << ( CoMSum / 1000 ) << endl << endl ;
             **/
            
            CoM.push_back( CoMSum / 1000 ) ;
            
            
            if ( L == ( originaldna.length() - 3 ) )
            {
                
                
                 cout << "Center of Mass = ";
                 
                 cout << "(" ;
                 for ( int i = 0 ; i < Dimension ; i++ )
                 {
                 cout << CoM[i]  ;
                 
                 if ( i < Dimension - 1 )
                 {
                 cout << ", " ;
                 }
                 
                 }
                 cout << ")" << endl << endl ;
                 
                 cout << "Center of Mass Distance from Original Sequence = Sqrt[" ;
               
                for ( int i = 0 ; i < Dimension ; i++ )
                {
                    D2 += ( CoM[i] - R0[i] ) * ( CoM[i] - R0[i] ) ;
                    
                     cout << " (" << CoM[i]<< " " << "- " <<  R0[i] << ")^2" ;
                     
                     if ( i < (Dimension - 1) )
                     {
                     cout << " +" ;
                     
                     }
                     
                }
                cout << " ]" ;
                cout << "= " << sqrt( D2 ) << endl ;
                
                //cout <<"sqrt( D2 ) = " << sqrt( D2 ) << endl ;
                
                ProteinSpace << sqrt( D2 ) << endl ;
                
                CoM.clear() ;
                
                D2 = 0 	;
            }
        }
        
        
        cout << "A: " << A << endl ;
    }
    file.close() ;
    
    //coordinatevaluephobicity.close() ;
    
    //AASE.close() ;
    
    dna.close() ;
    ProteinSpace.close() ;

    return 0 ;
    
}
