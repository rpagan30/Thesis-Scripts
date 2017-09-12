/* AAEntropy&Hydrophobicity.cpp

Description:
This program calculates the Shannon Entropy and Kyte-Dolittle scale hydrophobicity for the 
out file of nuneutralrobustness.cpp.

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
  1) Place script in the same folder as the desired "out" file produced by nuneutralrobustness.cpp.
  2) Compile it and run: $ g++ AAEntropy&Hydrophobicity.cpp -o EntHyd
                         $ ./EntHyd
  3) Two text files get generated "hydrophobicity.txt" and "AASEntropy". These contain the hydrophobicity 
     and AA Entropy by sequence position throughout all generations.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

using namespace std ;



double Calculator( double value )
{
    double Entropy = 0 ;
    
    if( value == 0 )
    {
        Entropy = 0 ;
    }
    else
    {
        Entropy = value * log10( value ) ;
    }			


    return Entropy ;	
}



int main()
{
    remove("entropyout.txt") ;
    remove("hydrophobicity.txt") ;
    remove("AASEntropy.txt") ;

    string codon = "" ;

    //amino acid counters
    double  Ile, Val, Leu, Phe, Cys, Met, Ala, Gly, Thr, Ser, Trp, Tyr, Pro, His, Glu, Gln, Asp, Asn, Lys, Arg, Stop, AAs, AASEntropy;
    //amino acid probabilities
    double IleProb, ValProb, LeuProb, PheProb, CysProb, MetProb, AlaProb, GlyProb, ThrProb, SerProb, TrpProb,TyrProb, ProProb, HisProb, GluProb, GlnProb, AspProb, AsnProb, LysProb, ArgProb ;

    //nucleotide counters
    double A = 0 , T = 0, G = 0, C = 0 , nuc = 0 , Unrecognized = 0 ;
      
    ifstream file ;
    ofstream entropyout ;
    ofstream hydrophobicity ;
    ofstream AASE ;
    vector<string> dnalist ;
    vector<string> aalist ;
    vector<double> hydro ;
   // double MSEntropy = 0 ;
    double kdHydrophobicity = 0, hydrosum = 0 ;
    string line ;

    int S = 0 ;
    
    
    file.open("out") ;
    //entropyout.open("entropyout.txt", ios::app) ;
    hydrophobicity.open("hydrophobicity.txt", ios::app ) ;
    AASE.open("AASEntropy", ios::app ) ;
    
    while( file >> line )
    {
        if ( line.length() > 209 )
        {
            dnalist.push_back(line) ;
        
            //cout << line << endl;
        
            line = "" ;
        }
    }
    
    for( int L = 0 ; L < 210 ; L++ )
    {
        kdHydrophobicity = hydrosum = 0 ;

        Ile = Val = Leu = Phe = Cys = Met = Ala = Gly = Thr = Ser = Trp = Tyr = Pro = His = Glu = Gln = Asp = Asn = Lys = Arg = AAs = 0 ;

        hydro.clear() ;
        
        for( int N = 0 ; N < 500000 ; N++ )
        {
            line = "" ;
            
            line = dnalist[N] ;
            
            S = L + 1 ;

            if ( L == 0 || S % 3 == 0 )
            {


                    //cout << line.substr(L,3) << endl ;
                    if      (line.substr( L , 3 ) == "TCA")    { Ser++ ; hydro.push_back(-0.8); }
                    else if (line.substr( L , 3 ) == "TCC")    { Ser++ ; hydro.push_back(-0.8); }
                    else if (line.substr( L , 3 ) == "TCG")    { Ser++ ; hydro.push_back(-0.8); }
                    else if (line.substr( L , 3 ) == "TCT")    { Ser++ ; hydro.push_back(-0.8); }
                    else if (line.substr( L , 3 ) == "TTC")    { Phe++ ; hydro.push_back(2.8);  }
                    else if (line.substr( L , 3 ) == "TTT")    { Phe++ ; hydro.push_back(2.8);  }
                    else if (line.substr( L , 3 ) == "TTA")    { Leu++ ; hydro.push_back(3.8);  }
                    else if (line.substr( L , 3 ) == "TTG")    { Leu++ ; hydro.push_back(3.8);  }
                    else if (line.substr( L , 3 ) == "TAC")    { Tyr++ ; hydro.push_back(-1.3); }
                    else if (line.substr( L , 3 ) == "TAT")    { Tyr ++; hydro.push_back(-1.3); }
                    else if (line.substr( L , 3 ) == "TAA")    { Stop++; hydro.push_back(0.0) ; }
                    else if (line.substr( L , 3 ) == "TAG")    { Stop++; hydro.push_back(0.0) ; }
                    else if (line.substr( L , 3 ) == "TGC")    { Cys++ ; hydro.push_back(2.5);  }
                    else if (line.substr( L , 3 ) == "TGT")    { Cys++ ; hydro.push_back(2.5);  }
                    else if (line.substr( L , 3 ) == "TGA")    { Stop++; hydro.push_back(0.0) ; }
                    else if (line.substr( L , 3 ) == "TGG")    { Trp++ ; hydro.push_back(-0.9); }
                    else if (line.substr( L , 3 ) == "CTA")    { Leu++ ; hydro.push_back(3.8);  }
                    else if (line.substr( L , 3 ) == "CTC")    { Leu++ ; hydro.push_back(3.8);  }
                    else if (line.substr( L , 3 ) == "CTG")    { Leu++ ; hydro.push_back(3.8);  }
                    else if (line.substr( L , 3 ) == "CTT")    { Leu++ ; hydro.push_back(3.8);  }
                    else if (line.substr( L , 3 ) == "CCA")    { Pro++ ; hydro.push_back(-1.6); }
                    else if (line.substr( L , 3 ) == "CCC")    { Pro++ ; hydro.push_back(-1.6); }
                    else if (line.substr( L , 3 ) == "CCG")    { Pro++ ; hydro.push_back(-1.6); }
                    else if (line.substr( L , 3 ) == "CCT")    { Pro++ ; hydro.push_back(-1.6); }
                    else if (line.substr( L , 3 ) == "CAC")    { His++ ; hydro.push_back(-3.2); }
                    else if (line.substr( L , 3 ) == "CAT")    { His++ ; hydro.push_back(-3.2); }
                    else if (line.substr( L , 3 ) == "CAA")    { Gln++ ; hydro.push_back(-3.5); }
                    else if (line.substr( L , 3 ) == "CAG")    { Gln++ ; hydro.push_back(-3.5); }
                    else if (line.substr( L , 3 ) == "CGA")    { Arg++ ; hydro.push_back(-4.5); }
                    else if (line.substr( L , 3 ) == "CGC")    { Arg++ ; hydro.push_back(-4.5); }
                    else if (line.substr( L , 3 ) == "CGG")    { Arg++ ; hydro.push_back(-4.5); }
                    else if (line.substr( L , 3 ) == "CGT")    { Arg++ ; hydro.push_back(-4.5); }
                    else if (line.substr( L , 3 ) == "ATA")    { Ile++ ; hydro.push_back(4.5);  }
                    else if (line.substr( L , 3 ) == "ATC")    { Ile++ ; hydro.push_back(4.5);  }
                    else if (line.substr( L , 3 ) == "ATT")    { Ile++ ; hydro.push_back(4.5);  }
                    else if (line.substr( L , 3 ) == "ATG")    { Met++ ; hydro.push_back(1.9);  }
                    else if (line.substr( L , 3 ) == "ACA")    { Thr++ ; hydro.push_back(-0.7); }
                    else if (line.substr( L , 3 ) == "ACC")    { Thr++ ; hydro.push_back(-0.7); }
                    else if (line.substr( L , 3 ) == "ACG")    { Thr++ ; hydro.push_back(-0.7); }
                    else if (line.substr( L , 3 ) == "ACT")    { Thr++ ; hydro.push_back(-0.7); }
                    else if (line.substr( L , 3 ) == "AAC")    { Asn++ ; hydro.push_back(-3.5); }
                    else if (line.substr( L , 3 ) == "AAT")    { Asn++ ; hydro.push_back(-3.5); }
                    else if (line.substr( L , 3 ) == "AAA")    { Lys++ ; hydro.push_back(-3.9); }
                    else if (line.substr( L , 3 ) == "AAG")    { Lys++ ; hydro.push_back(-3.9); }
                    else if (line.substr( L , 3 ) == "AGC")    { Ser++ ; hydro.push_back(-0.8); }
                    else if (line.substr( L , 3 ) == "AGT")    { Ser++ ; hydro.push_back(-0.8); }
                    else if (line.substr( L , 3 ) == "AGA")    { Arg++ ; hydro.push_back(-4.5); }
                    else if (line.substr( L , 3 ) == "AGG")    { Arg++ ; hydro.push_back(-4.5); }
                    else if (line.substr( L , 3 ) == "GTA")    { Val++ ; hydro.push_back(4.2);  }
                    else if (line.substr( L , 3 ) == "GTC")    { Val++ ; hydro.push_back(4.2);  }
                    else if (line.substr( L , 3 ) == "GTG")    { Val++ ; hydro.push_back(4.2);  }
                    else if (line.substr( L , 3 ) == "GTT")    { Val++ ; hydro.push_back(4.2);  }
                    else if (line.substr( L , 3 ) == "GCA")    { Ala++ ; hydro.push_back(1.8);  }
                    else if (line.substr( L , 3 ) == "GCC")    { Ala++ ; hydro.push_back(1.8);  }
                    else if (line.substr( L , 3 ) == "GCG")    { Ala++ ; hydro.push_back(1.8);  }
                    else if (line.substr( L , 3 ) == "GCT")    { Ala++ ; hydro.push_back(1.8);  }
                    else if (line.substr( L , 3 ) == "GAC")    { Asp++ ; hydro.push_back(-3.5); }
                    else if (line.substr( L , 3 ) == "GAT")    { Asp++ ; hydro.push_back(-3.5); }
                    else if (line.substr( L , 3 ) == "GAA")    { Glu++ ; hydro.push_back(-3.5); }
                    else if (line.substr( L , 3 ) == "GAG")    { Glu++ ; hydro.push_back(-3.5); }
                    else if (line.substr( L , 3 ) == "GGA")    { Gly++ ; hydro.push_back(-0.4); }
                    else if (line.substr( L , 3 ) == "GGC")    { Gly++ ; hydro.push_back(-0.4); }
                    else if (line.substr( L , 3 ) == "GGG")    { Gly++ ; hydro.push_back(-0.4); }
                    else if (line.substr( L , 3 ) == "GGT")    { Gly++ ; hydro.push_back(-0.4); }
                    else                                       { cout<<"Unrecognized codon: "<<codon<<endl; Unrecognized++; }
            }

            hydrosum += hydro[N] ;

        }
        
        kdHydrophobicity = hydrosum / 500000 ;

        //cout << kdHydrophobicity << endl ;
             
        //Amino Acid probabilities are calculated.
        AAs =  Ile + Val + Leu + Phe + Cys + Met + Ala + Gly + Thr + Ser + Trp + Tyr + Pro + His + Glu + Gln + Asp + Asn + Lys + Arg ;

        IleProb = Ile/AAs; ValProb= Val/AAs; LeuProb=Leu/AAs; PheProb=Phe/AAs; CysProb=Cys/AAs; MetProb=Met/AAs; AlaProb=Ala/AAs; GlyProb=Gly/AAs; ThrProb=Thr/AAs; SerProb=Ser/AAs;

        TrpProb=Trp/AAs;TyrProb=Tyr/AAs; ProProb=Pro/AAs; HisProb=His/AAs; GluProb=Glu/AAs; GlnProb=Gln/AAs; AspProb=Asp/AAs; AsnProb=Asn/AAs; LysProb=Lys/AAs; ArgProb=Arg/AAs ;

        //Amino Acid entropy is calculated.

        AASEntropy = -1 *  ( Calculator(IleProb) +  Calculator(ValProb) +  Calculator(LeuProb) +  Calculator(PheProb) +  Calculator(CysProb) +
                             Calculator(MetProb) +  Calculator(AlaProb) +  Calculator(GlyProb) +  Calculator(ThrProb) +  Calculator(SerProb) +
                             Calculator(TrpProb) +  Calculator(TyrProb) +  Calculator(ProProb) +  Calculator(HisProb) +  Calculator(GluProb) +
                             Calculator(GlnProb) +  Calculator(AspProb) +  Calculator(AsnProb) +  Calculator(LysProb) +  Calculator(ArgProb) )/log10(20) ;

        //Amino Acid entropy is sent to a file.
        if ( L == 0 || S % 3 == 0 )
        {
            AASE << AASEntropy << endl ;
            hydrophobicity << kdHydrophobicity << endl ;
        }

        
    }
    
    file.close() ;
    
    //entropyout.close() ;

    hydrophobicity.close() ;

    AASE.close() ;


    cout << Unrecognized << " unrecognized codons were found." << endl ;
    
    return 0 ;
    
}
