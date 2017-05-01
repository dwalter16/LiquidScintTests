/************************************************************************
 *
 *  Filename: NewScan.cpp
 *
 *  Description: 
 *
 *	Author(s):
 *     Michael T. Febbraro
 *     David Walter     
 *
 *  Creation Date: 9/25/2016
 *  Last modified: 4/20/2017
 * 
 *  To compile: g++ -O3 -pedantic -o Scan.exe `root-config --cflags --libs` -lSpectrum NewScan.cpp
 *      - if errors copiling in Mac OSX
 *        - remove -03 option
 *        - clang++ instead of g++
 *        - $(root-config --cflags --libs) instead of 'root-config --cflags --libs'
 *
 *
 *  If "error while loading shared libraries: libcore.so: ..." occurs, type
 *  "source `root-config --prefix`/bin/thisroot.sh"
 * 
 * -----------------------------------------------------
 * 	Nuclear Astrophysics, Physics Division
 *  Oak Ridge National Laboratory
 * 
 */


#include <iostream>
#include <fstream>
#include <string>
#include <signal.h>
#include "PulseAnalysis.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TSpectrum.h"

using namespace std;

typedef struct 
{
  Float_t l;               // Long integral sqrt(ll*lr)
  Float_t s;               // Short integral sqrt(sl*sr)
  Float_t amp;             // Amplitude
  Float_t time;            // Trigger time 
  Float_t psd;             // PSD parameter s/l
  Float_t trg;             // Det trigger
  Float_t ap;              // Bar after pulsing detected
  Float_t td;              // Time since last trigger  (us)
  Float_t pp;              // Position of max amplitude 
  Float_t np;              // Number of peaks 

} Det;

int NewScan (){
	
	
  /** ---------------------------------------------------- 
   *	Variable declairation
   *   ---------------------------------------------------- 
   */
  static Det det1, det2;
		
  bool	beamON,
    trg,
    data;
	
  float	X;
  int	npeaks, multi, ap;
	
  ifstream fp[16];
			 
  string 	line;
	
  int		i,j, pposition,
    Tracelength = 100;
			
  float 	pulse[500],
    baseline[500];
			
  Float_t	amplitude,
    risetime,
    falltime,
    width,
    CFD,
    paraL,
    paraS,
    runtime;
			
  char 	filename[250],		
    prompt[10],
    openfile[250],
    prefix[100],
    interrputPrompt;
			
  Float_t trgtime, prevtime, difftime;
  Float_t prevtrgtime[10];
  long	TEvt = 0;

  TH1F *phist1 = new TH1F("h1","h1",Tracelength,0,Tracelength-1);
  TH1F *phist2 = new TH1F("h2","h2",Tracelength,0,Tracelength-1);

  TSpectrum *s = new TSpectrum();
	
  /** ---------------------------------------------------- 
   *	Calibrations and threshold
   *   ---------------------------------------------------- 
   */
        
  float cal[5] = 
    {   1.0,
	1.0,
	1.0,
	1.0,
	1.0
    }; // calibration (keVee / (bit * sample)) from manual check on calibration
    
    
  float threshold[10] = 
    {   15900,
	15900,
	15860,
	15920,
	15900,
	15900,
	15930,
	15930,
	15850,
	15850
    };

  // Upper electron band gate
  // Extracted from 'psdgate.cpp' using 5-sigma
  float egate[5][3] =
    {
      {2.207, 8.76E-6, 0.108},
      {1.964, 1.27E-6, 0.132},
      {2.094, 4.74E-6, 0.122},
      {2.13E-3, -3.24E-6, 0.397},
      {1.302, -2.05E-5, 0.252}
    };
        

  /** ---------------------------------------------------- 
   *	Get functions
   *   ---------------------------------------------------- 
   */
	
  PulseAnalysis *Analysis = new PulseAnalysis();
		
    
  /** ---------------------------------------------------- 
   *	Program start...
   *   ---------------------------------------------------- 
   */

  cout << " ------------------------------------------------ " << endl;
  cout << " | NewScan.cpp                                   |" << endl; 
  cout << " |   Experiment: Liquid Scintillator Tests       |" << endl; 
  cout << " |   Date: April 20th 2017                       |" << endl; 
  cout << " |   Calibration used: none                      |" << endl; 
  cout << " |   ORNL Nuclear Astrophysics                   |" << endl;
  cout << " ------------------------------------------------ " << endl;
	
  cout << "Root file name to be created: ";
  cin >> filename;
    
  cout << "Run file name: ";
  cin >> prefix;

  TFile *ff = new TFile(filename,"RECREATE");

  TTree *tt = new TTree("T","Liquid Scintillator Tests");
	
  tt->Branch("det1",&det1,"l:s:amp:time:psd:trg:ap:td:pp:np");
  tt->Branch("det2",&det2,"l:s:amp:time:psd:trg:ap:td:pp:np");
  
  tt->Branch("runtime",&runtime,"runtime/F");     // Runtime in ms
    
  	
  // Open files
  for (i = 0; i < 1; i++)
  {
    //sprintf(openfile, "%s_wave%d.txt", prefix, i);
    sprintf(openfile, "./%s", prefix);
    
    cout << " Opening file: " << openfile << endl;
    fp[i].open(openfile, std::ifstream::in);
    
  }
    
  data = 1;
  runtime = 0;
  prevtime = 0;
    
  /** ---------------------------------------------------- 
   *	Process liquid scintillator det events
   *   ----------------------------------------------------  
   */
	
  while (data)
    {
      multi = 0;
      
      X = -1;
      ap = 0;
      beamON = 0;
      for (j = 0; j < 1; j++)
	{
	  if(!fp[j].is_open()){data = 0; cout << "Could not open file!!" << endl;}
	  if(fp[j].is_open())
	    {
	      if (!getline(fp[j], line)) {data=0; break; } // Record length
	      getline(fp[j], line); // Channel number
	      getline(fp[j], line); // Event number
	      getline(fp[j], line); // Trigger time stamp
     	      
	      // Trigger time in 2 ADC clock cycles ( 8 ns )
	      if (j == 0)
                trgtime = atof(line.substr(line.find_last_of(" "),line.length()).c_str());
            
	      // Reset variables
	      CFD = -1;
	      amplitude = -1;
	      paraL = -1;
	      paraS = -1;
	      trg = 0;
	      pposition = -1;
	      npeaks = 0;
            
	      // Get traces
	      for (i = 0; i < Tracelength; i++)
		{
		  if (!getline(fp[j], line)) {data = 0; break;}
		  pulse[i] = 16383 - atof(line.c_str());

		  switch(j) {
		  case 0 : phist1->SetBinContent(i, pulse[i]); break;
		  case 1 : phist2->SetBinContent(i, pulse[i]); break;
		  }
                
		  if (j < 10)
                    if((16838 - pulse[i]) <= threshold[j]) {trg = 1;}
		}

	      /** Liquid det processing **/
	      if(Tracelength > 1 && j < 10)
		{
		  // Process trace
		  Analysis->Baseline_restore(pulse, baseline, Tracelength, 10, 3); 
		  Analysis->Parameters(pulse, Tracelength, 3, &CFD, &amplitude, &risetime, &falltime, &width);
                
		  // Note this method has been updated to including after pulsing detection
		  Analysis->PSD_Integration_Afterpulsing(pulse, Tracelength, 12, 50, 5, 15, 50, &pposition, &paraL, &paraS, &ap);

		}
	      
	    }
	  

	  switch(j) {
	  case 0 : 
	    det1.amp = amplitude;
	    det1.l = paraL;
	    det1.s = paraS;
	    det1.time = CFD;
	    det1.trg = trg;
			                    
	    det1.ap = ap;
   
	    det1.l = sqrt(det1.l * det1.l)*cal[0];
	    det1.s = sqrt(det1.s * det1.s)*cal[0];
	    det1.psd = det1.s / det1.l;
	    if (det1.trg){ det1.trg = 1;}
	    else {det1.trg = 0;}

	    // Runtime clock
	    difftime = trgtime - prevtime;
	    if(difftime < 0) { 
	      difftime = difftime + 2147483647;
	      runtime += ((8*difftime)/1.0E6);
	      prevtime = trgtime;
	    }
	    else {                
	      runtime += ((8*difftime)/1.0E6);
	      prevtime = trgtime;
	    }
                    
	    //Time trigger difference
	    det1.td = trgtime - prevtrgtime[0];
	    if(det1.td < 0) { det1.td = (8*(det1.td + 2147483647))/1.0E3; prevtrgtime[0] = trgtime;}
	    else { det1.td = (8*(det1.td))/1.0E3; prevtrgtime[0] = trgtime;}
                    
	    break;

	  case 1 : 
	    det2.amp = amplitude;
	    det2.l = paraL;
	    det2.s = paraS;
	    det2.time = CFD;
	    det2.trg = trg;
			    
	    det2.ap = ap;
                    
	    det2.l = sqrt(det2.l * det2.l)*cal[1];
	    det2.s = sqrt(det2.s * det2.s)*cal[1];
	    det2.psd = det2.s / det2.l;
	    if (det2.trg){ det2.trg = 1;}
	    else {det2.trg = 0;}
	                     
	    //Time trigger difference
	    det2.td = trgtime - prevtrgtime[1];
	    if(det2.td < 0) { det2.td = (8*(det2.td + 2147483647))/1.0E3; prevtrgtime[1] = trgtime;}
	    else { det2.td = (8*(det2.td))/1.0E3; prevtrgtime[1] = trgtime;}

	    // Noise and after pulsing detection
	    if (det2.trg && det2.amp<15500 && det2.psd > (egate[0][0]/sqrt(det2.l) + egate[0][1]*det2.l + egate[0][2]))
	      {
		det2.np = s->Search(phist1, 2.0, "", 0.05);
	      }
	    else { det2.np = -1; }
       	    break;

	  
	  
	  
	  }
               
            
            
        }
      
      //if (TEvt > 100) {data = 0;}
	
      tt->Fill();
      TEvt++;
      if (TEvt%1000==0) {cout << "\rEvent counter: " << TEvt << flush;}
	
    }
	
	
  for (i = 0; i < 15; i++)
    {
      if(fp[j].is_open())
	fp[j].close();
    }
  ff->cd(); tt->Write(); 
  ff->Close();
	
  cout << "\nFinsihed! " << TEvt << " events" << endl; 
    
  return 0;
}

int main() {
  return NewScan();
}
