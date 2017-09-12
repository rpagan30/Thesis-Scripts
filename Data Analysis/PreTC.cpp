/* PreTC.cpp

Description:
This script calculates the pre-termination codon bias for all generations.

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
$ g++ PreTC.cpp -o PreTC
$ ./PreTC
3) The files "PreTC.txt" and "PreTCSLAG.txt" will be produced in the same folder.
*/


#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <fstream>
using namespace std ;



int main()
{

    remove ("PreTC.txt") ;
    remove ("PreTCSLAG.txt") ;
    remove ("mix.txt") ;
    vector<string> dnalist ;
    vector<string> unique ;
    ifstream file ;
    ofstream ofile ;
    ofstream ofile1 ;
    ofstream mix ;
    ofstream fpretcmix ;
    ofstream fnonpretcmix ;
    ofstream funique ;

    // Pre termination codon counters
    double cAAA = 0, cAAG = 0 , cAGA = 0 , cCAA = 0 ,cCAG = 0 ,cCGA = 0 , cCGG = 0 , cGAA = 0 , cGAG = 0 , cGGA = 0 , cTAT = 0 , cTAC = 0 , cTCA = 0 , cTCG = 0 , cTTA = 0 , cTTG = 0 ,
            cTGG = 0 , cTGC = 0 , cTGT = 0 ;

    // Non Pre TC counters
    double cTTC = 0 , cTTT = 0 , cATT = 0 , cATC = 0 , cATG = 0 , cATA = 0 , cCTT = 0 , cCTG = 0 , cCTA = 0 , cCTC = 0 , cCCT = 0 , cCCC = 0 , cCCA = 0 , cCCG = 0 , cACT = 0 , cACC = 0 ,
            cACA = 0 , cACG = 0 , cGCT = 0 , cGCC = 0 , cGCA = 0 , cGCG = 0 , cCAT = 0 , cCAC = 0 , cAAT = 0 , cAAC = 0 , cGAT = 0 , cGAC = 0 ;

    // SlAG Non Pre TC counters
    double cTCT = 0 , cTCC = 0 , cAGT = 0 , cAGC = 0 ,cAGG = 0 , cCGT = 0 , cCGC = 0 , cGGT = 0 ,
            cGGC = 0 , cGGG = 0 , cGTT = 0 , cGTC = 0 , cGTA = 0 , cGTG = 0 , cNonPreTC = 0 , cNonPreTCSLAG = 0 , cNonPreTCMix = 0 ;

    //Stop Codons and Unrecognized Codons
    double cTAA = 0 , cTAG = 0 , cTGA = 0 , Unrecognized = 0 ;

    // Calculation variables
    double preTCSum = 0 , preTCSumSLAG = 0 , preTCBias = 0 , preTCSLAGBias = 0 , preTCMix = 0 ,PreTCMixBias = 0 ;

    int uvals = 0 ;// Stores unique value count.

    string codon1 , DNAseq , line ;

    line = "" ;

    file.open("out") ;

    ofile.open( "PreTC.txt" , ios::app ) ;

    ofile1.open( "PreTCSLAG.txt" , ios::app ) ;

    //mix.open( "mix.txt" , ios::app ) ;

    //fpretcmix.open( "fpretcmix.txt" , ios::app ) ;

    //fnonpretcmix.open( "fnonpretcmix.txt" , ios::app ) ;

    //funique.open( "uniqueseqs.txt" , ios::app ) ;

    while( file >> line )
    {
        if ( line.length() > 209 )// Edit this for different DNA sequence lengths. Given by length - 1.
        {
            dnalist.push_back(line) ;

            line = "" ;
        }
    }

    //cout << "This is dnalist size " << dnalist.size() << endl ;

    for ( int i = 0 ; i < 500000 ; i ++ )//dnalist.size()
    {
        DNAseq = dnalist[i] ;

        //cout << dnalist[i] << endl ;

        for( int j = 0 ; j < dnalist[i].length() ; j += 3 )
        {
            codon1 = DNAseq.substr( j , 3 ) ;
            //cout << codon1 << endl ;

            // 11 Pre TCs
            if     ( codon1 == "AAA" ){cAAA++ ;}//PreTC // Lys
            else if( codon1 == "AAG" ){cAAG++ ;}//PreTC // Lys
            else if( codon1 == "CAA" ){cCAA++ ;}//PreTC // Gln
            else if( codon1 == "CAG" ){cCAG++ ;}//PreTC // Gln
            else if( codon1 == "GAA" ){cGAA++ ;}//PreTC // Glu
            else if( codon1 == "GAG" ){cGAG++ ;}//PreTC // Glu
            else if( codon1 == "TAT" ){cTAT++ ;}//PreTC // Tyr
            else if( codon1 == "TAC" ){cTAC++ ;}//PreTC // Tyr
            else if( codon1 == "TGG" ){cTGG++ ;}//PreTC // Trp
            else if( codon1 == "TGC" ){cTGC++ ;}//PreTC // Cys
            else if( codon1 == "TGT" ){cTGT++ ;}//PreTC // Cys

            // 8 SLAG Pre TCs
            else if( codon1 == "AGA" ){cAGA++ ;}//SLAG PreTC // Arg
            else if( codon1 == "CGA" ){cCGA++ ;}//SLAG PreTC // Arg
            else if( codon1 == "CGG" ){cCGG++ ;}//SLAG PreTC // Arg
            else if( codon1 == "GGA" ){cGGA++ ;}//SLAG PreTC // Gly
            else if( codon1 == "TCA" ){cTCA++ ;}//SLAG PreTC // Ser
            else if( codon1 == "TCG" ){cTCG++ ;}//SLAG PreTC // Ser
            else if( codon1 == "TTA" ){cTTA++ ;}//SLAG PreTC // Leu
            else if( codon1 == "TTG" ){cTTG++ ;}//SLAG PreTC // Leu

            // 24 Non Pre TCs
            else if( codon1 == "TTT" ){cTTT++ ;}// Non PreTC // Phe
            else if( codon1 == "TTC" ){cTTC++ ;}// Non PreTC // Phe
            else if( codon1 == "ATT" ){cATT++ ;}// Non PreTC // Ile
            else if( codon1 == "ATC" ){cATC++ ;}// Non PreTC // Ile
            else if( codon1 == "ATG" ){cATG++ ;}// Non PreTC // Met
            else if( codon1 == "ATA" ){cATA++ ;}// Non PreTC // Ile
            else if( codon1 == "CCT" ){cCCT++ ;}// Non PreTC // Pro
            else if( codon1 == "CCC" ){cCCC++ ;}// Non PreTC // Pro
            else if( codon1 == "CCA" ){cCCA++ ;}// Non PreTC // Pro
            else if( codon1 == "CCG" ){cCCG++ ;}// Non PreTC // Pro
            else if( codon1 == "ACT" ){cACT++ ;}// Non PreTC // Thr
            else if( codon1 == "ACC" ){cACC++ ;}// Non PreTC // Thr
            else if( codon1 == "ACA" ){cACA++ ;}// Non PreTC // Thr
            else if( codon1 == "ACG" ){cACG++ ;}// Non PreTC // Thr
            else if( codon1 == "GCT" ){cGCT++ ;}// Non PreTC // Ala
            else if( codon1 == "GCC" ){cGCC++ ;}// Non PreTC // Ala
            else if( codon1 == "GCA" ){cGCA++ ;}// Non PreTC // Ala
            else if( codon1 == "GCG" ){cGCG++ ;}// Non PreTC // Ala
            else if( codon1 == "CAT" ){cCAT++ ;}// Non PreTC // His
            else if( codon1 == "CAC" ){cCAC++ ;}// Non PreTC // His
            else if( codon1 == "AAT" ){cAAT++ ;}// Non PreTC // Asn
            else if( codon1 == "AAC" ){cAAC++ ;}// Non PreTC // Asn
            else if( codon1 == "GAT" ){cGAT++ ;}// Non PreTC // Asp
            else if( codon1 == "GAC" ){cGAC++ ;}// Non PreTC // Asp

            // 13 SLAG Non Pre TCs
            else if( codon1 == "TCT" ){cTCT++ ;}// SLAG Non PreTC // Ser
            else if( codon1 == "TCC" ){cTCC++ ;}// SLAG Non PreTC // Ser
            else if( codon1 == "AGT" ){cAGT++ ;}// SLAG Non PreTC // Ser
            else if( codon1 == "AGC" ){cAGC++ ;}// SLAG Non PreTC // Ser
            else if( codon1 == "CTT" ){cCTT++ ;}// SLAG Non PreTC // Leu
            else if( codon1 == "CTG" ){cCTG++ ;}// SLAG Non PreTC // Leu
            else if( codon1 == "CTA" ){cCTA++ ;}// SLAG Non PreTC // Leu
            else if( codon1 == "CTC" ){cCTC++ ;}// SLAG Non PreTC // Leu
            else if( codon1 == "CGT" ){cCGT++ ;}// SLAG Non PreTC // Arg
            else if( codon1 == "CGC" ){cCGC++ ;}// SLAG Non PreTC // Arg
            else if( codon1 == "GGT" ){cGGT++ ;}// SLAG Non PreTC // Gly
            else if( codon1 == "GGC" ){cGGC++ ;}// SLAG Non PreTC // Gly
            else if( codon1 == "GGG" ){cGGG++ ;}// SLAG Non PreTC // Gly

            //14th SLAG Non Pre TC
            else if( codon1 == "AGG" ){cAGG++ ;}// SLAG Non PreTC  // Arg               //Not Mentioned Anywhere

            // 4 Unaccounted Codons //Val Amino Acid
            else if( codon1 == "GTT" ){cGTT++ ;}// Non PreTC // Val                     //Not Mentioned Anywhere // Added to Non PreTCs
            else if( codon1 == "GTC" ){cGTC++ ;}// Non PreTC // Val                     //Not Mentioned Anywhere // Added to Non PreTCs
            else if( codon1 == "GTA" ){cGTA++ ;}// Non PreTC // Val                     //Not Mentioned Anywhere // Added to Non PreTCs
            else if( codon1 == "GTG" ){cGTG++ ;}// Non PreTC // Val                     //Not Mentioned Anywhere // Added to Non PreTCs

            // Stop Codons and Unrecognized Codons
            else if( codon1 == "TAA" ){cTAA++ ;}// Stop Codon
            else if( codon1 == "TAG" ){cTAG++ ;}// Stop Codon
            else if( codon1 == "TGA" ){cTGA++ ;}// Stop Codon
            else                      {Unrecognized++ ; cout << codon1 << "Unrecognized codon found." <<  endl ; }// Unrecognized Codons






        }

        if ( (i+1) % 1000 == 0 )
        {


            //Calculations

            preTCSum = cAAA + cAAG + cCAA + cCAG + cGAA + cGAG + cTAT + cTAC + cTGG + cTGC + cTGT ;//11

            //cout << preTCSum << "preTC Sum" << endl ;

            preTCSumSLAG = cAGA + cCGG + cCGA + cGGA + cTCA + cTCG + cTTA + cTTG ;//8

            // Used to be : //cNonPreTC =  cNonPreTC + cTTC + cTTT + cATT + cATC + cATG + cATA + cCTT + cCTG + cCTA + cCTC + cCCT
            // + cCCC + cCCA + cCCG + cACT + cACC + cACA + cACG + cGCT + cGCC + cGCA + cGCG + cCAT + cCAC + cAAT + cAAC + cGAT + cGAC ;//24

            cNonPreTC = cTTC + cTTT + cATT + cATC + cATG + cATA + cCTT + cCTG + cCTA + cCTC + cCCT + cCCC + cCCA + cCCG +
                    cACT + cACC + cACA + cACG + cGCT + cGCC + cGCA + cGCG + cCAT + cCAC + cAAT + cAAC + cGAT + cGAC + cGTT + cGTC + cGTA + cGTG ;//28

            // Used to be : cNonPreTCSLAG =   CNonPreTCSLAG + cTCT + cTCC + cAGT + cAGC + cCGT + cCGC + cGGT + cGGC + cGGG + cGTT + cGTC + cGTA + cGTG ;//13

            cNonPreTCSLAG = cTCT + cTCC + cAGT + cAGC + cCGT + cCGC + cGGT + cGGC + cGGG + cGTT + cGTC + cGTA + cGTG + cAGG ;//14


            preTCBias = ( 28.0 / 11.0 ) * ( preTCSum / cNonPreTC ) ;

            preTCSLAGBias = ( 14.0 / 8.0 ) * ( preTCSumSLAG / cNonPreTCSLAG ) ;



            //preTCMix = preTCSum + preTCSumSLAG ;

            //cNonPreTCMix = cNonPreTC + cNonPreTCSLAG ;

            //PreTCMixBias = ( 42.0 / 19.0 ) * ( preTCMix / cNonPreTCMix ) ;


            ofile << preTCBias << endl ;

            ofile1 << preTCSLAGBias << endl ;

            //mix << PreTCMixBias << endl ;

            //fpretcmix << preTCMix << endl ;

            //fnonpretcmix << cNonPreTCMix << endl ;


            preTCSum = 0 ;

            preTCSumSLAG = 0 ;

            cNonPreTC = 0 ;

            cNonPreTCSLAG = 0 ;

            //preTCMix = 0 ;

            //cNonPreTCMix = 0 ;

            cAAA = cAAG = cAGA = cCAA = cCAG = cCGA = cCGG = cGAA = cGAG = cGGA = cTAT = cTAC = cTCA = cTCG = cTTA = cTTG = cTGG = cTGC = cTGT = 0 ;

            cTTC = cTTT = cATT = cATC = cATG = cATA = cCTT = cCTG = cCTA = cCTC = cCCT = cCCC = cCCA = cCCG = cACT = cACC = cACA = cACG = cGCT = cGCC = cGCA = cGCG = cCAT = cCAC = cAAT = cAAC = cGAT = cGAC = 0 ;

            cTCT = cTCC = cAGT = cAGC = cAGG =  cCGT = cCGC = cGGT = cGGC = cGGG = cGTT = cGTC = cGTA = cGTG = 0 ;

        }



    }

    cout << "PreTCBias: " << preTCBias << endl ;
    cout << "preTCSum: " << preTCSum << endl ;
    cout << "Non PreTC Bias: " << cNonPreTC << endl ;

    cout << Unrecognized << " unrecognized codons were found." << endl ;

    cout << "Stop Codons found: " << endl << "TAA: " << cTAA << endl
         << "TAG: " << cTAG << endl << "TGA: " << cTGA ;


    funique << "Unique Sequence Count: " << uvals << endl ;


    ofile.close() ;

    ofile1.close() ;

    file.close() ;

    //mix.close() ;

    //fpretcmix.close() ;

    //fnonpretcmix.close() ;

    //funique.close() ;

    return 0 ;

}
