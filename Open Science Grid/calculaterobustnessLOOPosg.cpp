#include<iostream>
#include<cmath>
#include<stdio.h>
#include<stdlib.h>
#include<fstream>
#include<string>
#include<iomanip>
#include<vector>
#include <sstream>

// this program will calculate the average delta delta G of a protein mutated at every position 
// need to give an argument for the output file
// will automatically generate the mutsizes file, which has a list of the ddG values
// the program will loop according to the size of the input file, this file is composed of 
// dna sequences

// robustness can be calculated for syn/nonsyn mutations or just nonsyn mutations

using namespace std;

char AA[20] = {'A','E','Q','D','N','L','G','K','S','V','R','T','P','I','M','F','Y','C','W','H'};


/////////////////the following generates the energy for the protein/////////////////////////

float energy(vector<float>& x, vector<float>& y, vector<float>& z,vector<int>& num, int n, vector<char>& seq, float contactenergy[20][20])
{
int aacontacts[1000][1000] = {0}; // this initializes the array as 0
float dist = 0.0;
float eng = 0.0;
int f=0,m=0;

//n = the number of lines that have x, y, z coordinates
//20.25 refers to 4.5 Angstroms squared (dist originally was a square root)

	for(int i=0;i<n;i++)
	{ 
		for(int j=0;j<n;j++)
		{
		dist = (x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]);
		if(dist<=20.25 && abs(num[i]-num[j])>=2) aacontacts[num[i]][num[j]] = 1;
         	
		}
	}

/////extract the amino acid sequence from num and seq////////////

char aminoacidseq[56];
int p=0;

for(int q=0;q<num.size();q++)
   {
    if(num[q]!=num[q+1] || num.size() < q + 2 )
      {aminoacidseq[p] = seq[q];
       p++;}
   }

//cout<<aminoacidseq[69]<<endl;
//////////calculate energy//////////////////////////////////////////////////////////

	for(int i=1;i<57;i++)
	{
		for(int j=1;j<57;j++)
		{
			if(aacontacts[i][j] == 1)

			{
                         //cout<<i<<" "<<j<<endl;
				for(int l=0;l<20;l++)
				{
					if(AA[l] == aminoacidseq[i-1])
						m = l;
					if(AA[l] == aminoacidseq[j-1])
						f = l;					
				}

				eng=eng+contactenergy[m][f];
				m=0;
				f=0;
			}
		}
	}
        eng = eng/2;
	return eng;
}

/////////////the following translates the DNA sequence/////////////////////////////////////////

char *translate(string codon)
{

	if      (codon == "TCA")    { return "s";}   
    	else if (codon == "TCC")    { return "s";}   
	else if (codon == "TCG")    { return "s";}  
    	else if (codon == "TCT")    { return "s";}    
    	else if (codon == "TTC")    { return "f";}
    	else if (codon == "TTT")    { return "f";}    
    	else if (codon == "TTA")    { return "l";} 
    	else if (codon == "TTG")    { return "l";}   
    	else if (codon == "TAC")    { return "y";} 
    	else if (codon == "TAT")    { return "y";}  
    	else if (codon == "TAA")    { return "_";} 
    	else if (codon == "TAG")    { return "_";} 
    	else if (codon == "TGC")    { return "c";}  
    	else if (codon == "TGT")    { return "c";}  
    	else if (codon == "TGA")    { return "_";}  
    	else if (codon == "TGG")    { return "w";} 
    	else if (codon == "CTA")    { return "l";}  
    	else if (codon == "CTC")    { return "l";}  
    	else if (codon == "CTG")    { return "l";}  
    	else if (codon == "CTT")    { return "l";}  
    	else if (codon == "CCA")    { return "p";}   
    	else if (codon == "CCC")    { return "p";}   
    	else if (codon == "CCG")    { return "p";}   
    	else if (codon == "CCT")    { return "p";}  
    	else if (codon == "CAC")    { return "h";}  
    	else if (codon == "CAT")    { return "h";}  
    	else if (codon == "CAA")    { return "q";}   
    	else if (codon == "CAG")    { return "q";}   
    	else if (codon == "CGA")    { return "r";}   
    	else if (codon == "CGC")    { return "r";}  
    	else if (codon == "CGG")    { return "r";}   
    	else if (codon == "CGT")    { return "r";}   
    	else if (codon == "ATA")    { return "i";}  
    	else if (codon == "ATC")    { return "i";}   
    	else if (codon == "ATT")    { return "i";}   
    	else if (codon == "ATG")    { return "m";}   
    	else if (codon == "ACA")    { return "t";}    
    	else if (codon == "ACC")    { return "t";}    
    	else if (codon == "ACG")    { return "t";}   
    	else if (codon == "ACT")    { return "t";}   
    	else if (codon == "AAC")    { return "n";}    
    	else if (codon == "AAT")    { return "n";}    
    	else if (codon == "AAA")    { return "k";}   
    	else if (codon == "AAG")    { return "k";}   
    	else if (codon == "AGC")    { return "s";}   
    	else if (codon == "AGT")    { return "s";}   
    	else if (codon == "AGA")    { return "r";}   
    	else if (codon == "AGG")    { return "r";}  
    	else if (codon == "GTA")    { return "v";}   
    	else if (codon == "GTC")    { return "v";}   
    	else if (codon == "GTG")    { return "v";}    
    	else if (codon == "GTT")    { return "v";}   
    	else if (codon == "GCA")    { return "a";}    
    	else if (codon == "GCC")    { return "a";}   
    	else if (codon == "GCG")    { return "a";}    
    	else if (codon == "GCT")    { return "a";}    
    	else if (codon == "GAC")    { return "d";} 
    	else if (codon == "GAT")    { return "d";}  
    	else if (codon == "GAA")    { return "e";}  
    	else if (codon == "GAG")    { return "e";}    
    	else if (codon == "GGA")    { return "g";}   
    	else if (codon == "GGC")    { return "g";}    
    	else if (codon == "GGG")    { return "g";}    
    	else if (codon == "GGT")    { return "g";}   
    	else
	{
	        cout<<"Unrecognized codon: "<<codon<<endl;
    	}

}

////////////////////////this is the main procedure///////////////////////////////////////////

int main(int argc, char *argv[])
{

// PDB file is the first argument
   string pdb_file;
if (argc < 2) {
   cout << "Not given pdb file to run.  Please give pdb as argument" << endl;
   exit(1);
} else {
   pdb_file = argv[1];
}

ofstream robustout("robustout");
ofstream mutsizes("mutsizes");
ofstream outfile("output");
ifstream energyfile;
float origenergy;
float contactenergy[20][20]; 
energyfile.open("vendruscolomatrix"); 

for (int i=0;i<20;i++)
    {
	for(int j=0;j<20;j++)
    	{
      		energyfile>>contactenergy[i][j];
    	}
     }
energyfile.close();



//// the input file  will be composed of DNA sequences ////////////////
ifstream stablednaseqs;
stablednaseqs.open("dnaseqs");
string ssline;

if(stablednaseqs.is_open()){cout << "helloe" << endl; }
{if(stablednaseqs.eof()){return 0;}
 while(!stablednaseqs.eof())
	{
	
	ofstream aastable("aastable");
	getline (stablednaseqs, ssline);
	string sscodon, ssprotein;

	for(int d=0; d<(ssline.length()-2);d += 3)
       		{
		sscodon = ssline.substr(d,3);
       		if(sscodon == "TAA"|| sscodon == "TAG" || sscodon =="TGA")
	  		{
	  		break;
	  		}
        	else ssprotein += translate(sscodon); 
   		}


	aastable<<ssprotein<<endl;
        
	/////////the following refines the initial pdb file////////////////////////////////////////////////
    stringstream comm;
    comm << "./scwrl4_lin/Scwrl4 -i " << pdb_file << " -o refined.pdb -s aastable -v -h > logfile";
	system(comm.str().c_str());

	//////the following extracts information from the pdb file ie. x,y,z, amino acid seq and numbers///

	ifstream pdbfile;
	pdbfile.open("refined.pdb");

	std::vector<int> num;
	std::vector<float> X;
	std::vector<float> Y;
	std::vector<float> Z;
	std::vector<char> aaseq;
	string line;

  	if (pdbfile.is_open())
  	{
    		while (! pdbfile.eof())
    		{
      		getline (pdbfile,line);
      		int Epos = line.find("ATOM",0);
      		int Epos2 = line.find("A",21);

      
       		if(Epos2 != string::npos && Epos != string::npos)
       		{
      		// cout<<line<<endl;
       		string str = line.substr(17,3);
      
       		if(str == "ALA"){/*cout<<"A"<<endl;*/ aaseq.push_back('A');}
       		if(str == "ARG"){/*cout<<"R"<<endl;*/ aaseq.push_back('R');}
       		if(str == "ASN"){/*cout<<"N"<<endl;*/ aaseq.push_back('N');}
       		if(str == "ASP"){/*cout<<"D"<<endl;*/ aaseq.push_back('D');}
       		if(str == "CYS"){/*cout<<"C"<<endl;*/ aaseq.push_back('C');}
       		if(str == "GLN"){/*cout<<"Q"<<endl;*/ aaseq.push_back('Q');}
       		if(str == "GLU"){/*cout<<"E"<<endl;*/ aaseq.push_back('E');}
       		if(str == "GLY"){/*cout<<"G"<<endl;*/ aaseq.push_back('G');}
       		if(str == "HIS"){/*cout<<"H"<<endl;*/ aaseq.push_back('H');}
       		if(str == "ILE"){/*cout<<"I"<<endl;*/ aaseq.push_back('I');}
       		if(str == "LEU"){/*cout<<"L"<<endl;*/ aaseq.push_back('L');}
       		if(str == "LYS"){/*cout<<"K"<<endl;*/ aaseq.push_back('K');}
       		if(str == "MET"){/*cout<<"M"<<endl;*/ aaseq.push_back('M');}
       		if(str == "PHE"){/*cout<<"F"<<endl;*/ aaseq.push_back('F');}
       		if(str == "PRO"){/*cout<<"P"<<endl;*/ aaseq.push_back('P');}
       		if(str == "SER"){/*cout<<"S"<<endl;*/ aaseq.push_back('S');}
       		if(str == "THR"){/*cout<<"T"<<endl;*/ aaseq.push_back('T');}
       		if(str == "TRP"){/*cout<<"W"<<endl;*/ aaseq.push_back('W');}
       		if(str == "TYR"){/*cout<<"Y"<<endl;*/ aaseq.push_back('Y');}
       		if(str == "VAL"){/*cout<<"V"<<endl;*/ aaseq.push_back('V');}

 	       string seqnum = line.substr(23, 3);
   	       for (int i = 0; i < seqnum.length(); i++)
     		     {
    	          if (seqnum[i] == ' ')
        	      {
         	     seqnum.replace(i,1,"");
         	     i = i - 1;
       	    		   } 
       		     }
      		int seqnumint = atoi(seqnum.c_str());
     		num.push_back(seqnumint);


      	 	string x = line.substr(31, 7);
         	 for (int i = 0; i < x.length(); i++)
        	  {
          	    if (x[i] == ' ')
          	    {
          	    x.replace(i,1,"");
          	    i = i - 1;
          	    } 
       	          }
      		float xint = atof(x.c_str());
      		X.push_back(xint);

 
  	       string y = line.substr(39,7);
   	       for (int i = 0; i < y.length(); i++)
   	       {
      	        if (y[i] == ' ')
      	        {
       	         y.replace(i,1,"");
        	 i = i - 1;
         	 } 
               }
    	       float yint = atof(y.c_str());
               Y.push_back(yint);


      		string z = line.substr(47,7);
         	 for (int i = 0; i < z.length(); i++)
        	  {
        	     if (z[i] == ' ')
          	     {
         	     z.replace(i,1,"");
          	     i = i - 1;
          	     } 
                  }
                float zint = atof(z.c_str());
                Z.push_back(zint);

               // cout<<seqnum<<"x"<<endl;
       	   	}
      	}
    	pdbfile.close();
    	//system("rm refined.pdb");
   	system("rm aastable");
        cout<<"initial aaseq size = "<<aaseq.size()<<endl;
   	}

	/////////////////////the following calls the energy calculation/////////////////////////////////

	origenergy = energy(X,Y,Z,num,aaseq.size(),aaseq, contactenergy);
	cout<<"initial refined energy = "<<origenergy<<endl;
	outfile<<"initial refined energy = "<<origenergy<<endl;
	///////////////the next part calculates mutational robustness//////////////////////////////////

	string dnainput, dnainput2, cod, dnanew, mutatedcodon;

	dnainput = ssline;

	cout<<"DNA Seq : "<<dnainput<<endl;
	outfile<<"DNA Seq : "<<dnainput<<endl;
	cout<<"AA Seq : "<<ssprotein<<endl;
	outfile<<"AA Seq : "<<ssprotein<<endl;

	//////////put the point mutated sequences into an array////////////////////////////////////////

	std::vector<string> dnalist;
	mutatedcodon ="";

	for(int i=0;i<(dnainput.length()-2);i += 3)
   	{
   	cod = dnainput.substr(i,3);
   	string one = cod.substr(0,1);
   	string two = cod.substr(1,1); 
   	string three = cod.substr(2,1);


	//////////first codon position:


	if(cod[0] == 'A')
 	  {
          mutatedcodon = "T" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = "C" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = "G" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

	if(cod[0] == 'T')
 	  {
          mutatedcodon = "A" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = "C" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = "G" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

    	if(cod[0] == 'G')
 	  {
          mutatedcodon = "A" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = "C" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = "T" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

   	if(cod[0] == 'C')
 	  {
          mutatedcodon = "A" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = "G" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = "T" + two + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }


	///second codon position:

	if(cod[1] == 'A')
 	  {
          mutatedcodon = one + "T" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one +"C" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + "G" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

	if(cod[1] == 'T')
 	  {
          mutatedcodon = one + "A" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + "C" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + "G" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

    	if(cod[1] == 'G')
 	  {
          mutatedcodon = one + "A" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + "C" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + "T" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

   	if(cod[1] == 'C')
 	  {
          mutatedcodon = one + "A" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + "G" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + "T" + three; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

	/////third codon position

	if(cod[2] == 'A')
 	  {
          mutatedcodon = one + two + "T"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + two + "C"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + two + "G"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

	if(cod[2] == 'T')
 	  {
          mutatedcodon = one + two + "A"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + two + "C"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + two + "G"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

    	if(cod[2] == 'G')
 	  {
          mutatedcodon = one + two + "A"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + two + "C"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + two + "T"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }

   	if(cod[2] == 'C')
 	  {
          mutatedcodon = one + two + "A"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
          
          mutatedcodon = one + two + "G"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
  	  
          mutatedcodon = one + two + "T"; 
          dnainput2 = dnainput;
          dnanew = dnainput2.replace(i,3,mutatedcodon);
  	  dnalist.push_back(dnanew);
          dnanew = "";
          mutatedcodon = "";
    	  }
    }

	cout<<"dnalist size = "<<dnalist.size()<<endl;
	//outfile<<"dnalist size = "<<dnalist.size()<<endl;

	///////////////////translate each DNA sequence//////////////////////////////////////////////

	std::vector<string> mutprotein;
	string codon, protein, dnaseq;

	for(int y=0; y<dnalist.size();y++)
    	{	
    	dnaseq = "";protein = "";
    	dnaseq = dnalist[y];
    	for(int z=0; z<(dnaseq.length()-2);z += 3)
       		{
		codon = dnaseq.substr(z,3);
       		if(codon == "TAA"|| codon == "TAG" ||codon =="TGA")
	  		{
	  		break;
	  		}
        	else protein += translate(codon); 
   		}

  	 mutprotein.push_back(protein);

	// the following code removes synonymous mutations from the calculation of robustness ///

	//  if(protein != ssprotein){mutprotein.push_back(protein);}

    }

	///////////////////calculate the energy of each mutant protein////////////////////////////////////

	std::vector<float> energyarray;
	string proteinseq;
	float energyoutput;
	string str;

	cout<<"Number of mutated proteins = "<<mutprotein.size()<<endl;

	for(int h=0; h<mutprotein.size();h++)
   	{
                proteinseq = "";
                proteinseq = mutprotein[h];
                //cout<<proteinseq<<endl;

                if(proteinseq.length() < (dnainput.length()/3))
                {
                energyoutput = 0.0; cout<<"Truncated protein "<<energyoutput<<endl;
                //outfile<<"Truncated protein "<<energyoutput<<endl;
                //energyarray.push_back(energyoutput);
                }  

                else if(proteinseq == ssprotein)
                {
                energyoutput = origenergy; //cout<<"Synonymous mutation"<<energyoutput<<endl;
                //outfile<<"Synonymous mutation: "<<endl;
                //outfile<<proteinseq<<endl;
                //outfile<<energyoutput<<endl;
                energyarray.push_back(energyoutput);
                }        

                else
                {
                ofstream outseq("mutatedprotein");
                outseq<<proteinseq<<endl;
stringstream comm;
comm << "./scwrl4_lin/Scwrl4 -i " << pdb_file << " -o nurefined.pdb -s mutatedprotein -v -h > logfile";
system(comm.str().c_str());
                //system("./scwrl4_lin/Scwrl4 -i 1ROP.pdb -o nurefined.pdb -s mutatedprotein -v -h > logfile");

                outseq.clear();
                outseq.close();
                
                //system("rm mutatedprotein");
		//cout<<"Size of X = "<<X.size()<<endl;                
		X.clear();
                Y.clear();
                Z.clear();
                num.clear();
		aaseq.clear();

                line = "";
                energyoutput = 0.0;

		ifstream pdbfile;
		pdbfile.open("nurefined.pdb");                
                //ofstream outx("x");
                //ofstream outy("y");
                //ofstream outz("z");
                //ofstream outnum("num");

                if(!pdbfile)
                  {
                   cerr << "Unable to open pdb file";
                   exit(1);
                   }


 		 if(pdbfile.is_open())
  		   { //outfile<<"File is open"<<endl;
    		  while(! pdbfile.eof())
    		     {
      		    getline (pdbfile,line);
      		    int Epos = line.find("ATOM",0);
       		    int Epos2 = line.find("A",21);

       			if(Epos2 != string::npos && Epos != string::npos)
       			  {
           	          // cout<<line<<endl;
                          str = "";
       			  str = line.substr(17,3);
      
     			  if(str == "ALA"){/*cout<<"A"<<endl;*/ aaseq.push_back('A');}
    			  if(str == "ARG"){/*cout<<"R"<<endl;*/ aaseq.push_back('R');}
       			  if(str == "ASN"){/*cout<<"N"<<endl;*/ aaseq.push_back('N');}
       		 	  if(str == "ASP"){/*cout<<"D"<<endl;*/ aaseq.push_back('D');}
       		 	  if(str == "CYS"){/*cout<<"C"<<endl;*/ aaseq.push_back('C');}
       		 	  if(str == "GLN"){/*cout<<"Q"<<endl;*/ aaseq.push_back('Q');}
      		 	  if(str == "GLU"){/*cout<<"E"<<endl;*/ aaseq.push_back('E');}
                 	  if(str == "GLY"){/*cout<<"G"<<endl;*/ aaseq.push_back('G');}
       			  if(str == "HIS"){/*cout<<"H"<<endl;*/ aaseq.push_back('H');}
       			  if(str == "ILE"){/*cout<<"I"<<endl;*/ aaseq.push_back('I');}
       			  if(str == "LEU"){/*cout<<"L"<<endl;*/ aaseq.push_back('L');}
      			  if(str == "LYS"){/*cout<<"K"<<endl;*/ aaseq.push_back('K');}
       		 	  if(str == "MET"){/*cout<<"M"<<endl;*/ aaseq.push_back('M');}
      		 	  if(str == "PHE"){/*cout<<"F"<<endl;*/ aaseq.push_back('F');}
      		  	  if(str == "PRO"){/*cout<<"P"<<endl;*/ aaseq.push_back('P');}
      			  if(str == "SER"){/*cout<<"S"<<endl;*/ aaseq.push_back('S');}
       			  if(str == "THR"){/*cout<<"T"<<endl;*/ aaseq.push_back('T');}
       			  if(str == "TRP"){/*cout<<"W"<<endl;*/ aaseq.push_back('W');}
       			  if(str == "TYR"){/*cout<<"Y"<<endl;*/ aaseq.push_back('Y');}
       			  if(str == "VAL"){/*cout<<"V"<<endl;*/ aaseq.push_back('V');}

 	       		string seqnum = line.substr(23, 3);
   	      		 for(int i = 0; i < seqnum.length(); i++)
     			    {
    	      		    if(seqnum[i] == ' ')
        		      {
         		      seqnum.replace(i,1,"");
         		      i = i - 1;
       	    		      } 
       			     }
      			int seqnumint = atoi(seqnum.c_str());
     			num.push_back(seqnumint);
                        
                        //outnum<<seqnumint<<endl;

      	 		string x = line.substr(31, 7);
         		for(int i = 0; i < x.length(); i++)
        	           {
          	          if(x[i] == ' ')
          	            {
          	            x.replace(i,1,"");
          	            i = i - 1;
          	            } 
       	                   }
      			float xint = atof(x.c_str());
      			X.push_back(xint);

                       //outx<<xint<<endl;

 
  	    	        string y = line.substr(39,7);
   	    	        for(int i = 0; i < y.length(); i++)
   	                   {
      	                   if(y[i] == ' ')
      	                     {
       	                     y.replace(i,1,"");
        	             i = i - 1;
         	             } 
                           }
    	       		float yint = atof(y.c_str());
              	        Y.push_back(yint);

                       // outy<<yint<<endl;

      			string z = line.substr(47,7);
                  	 for(int i = 0; i < z.length(); i++)
        		    {
        	            if(z[i] == ' ')
          	              {
         	              z.replace(i,1,"");
          	              i = i - 1;
          	              } 
                            }
                        float zint = atof(z.c_str());
                        Z.push_back(zint);

                       // outz<<zint<<endl;
                       // cout<<seqnum<<"x"<<endl;
          	        }
      	             }
                pdbfile.clear();
    	        pdbfile.close();
                system("rm nurefined.pdb");

              //  cout<<"aaseqsize = " << aaseq.size() << endl;
             //   cout<<"num size = "<<num.size()<<endl;
   	        }
             energyoutput = energy(X,Y,Z,num,aaseq.size(),aaseq, contactenergy);  
             energyarray.push_back(energyoutput);
 
             //outx.close();
             //outy.close();
             //outz.close();
   
             //outfile<<proteinseq<<endl;
             //outfile<<energyoutput<<endl;

             }
    	}

	cout<<"Size of energy array = "<<energyarray.size()<<endl;

	//////////////////////calculate the average delta delta G of a point mutation for a protein///////

	float ddG = 0.0;
	float avddG;

	for(int b=0; b<energyarray.size(); b++)
    	{
    	ddG = ddG + abs(origenergy - energyarray[b]); 
    	mutsizes<<abs(origenergy - energyarray[b])<<endl;
    	}

	avddG = ddG / energyarray.size();

	cout<<"average delta delta G = "<<avddG<<endl;
	robustout<<"average delta delta G = "<<avddG<<endl;
	outfile<<"average delta delta G = "<<avddG<<"\n"<<endl;
	aastable.close();


     }

}
robustout.close();
outfile.close();

return 0;
}

