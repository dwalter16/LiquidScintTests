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
  Float_t cfd;             // Trigger time 
  Float_t psd;             // PSD parameter s/l
  Float_t trg;             // Det trigger
  Float_t tac;             // Zero-crossing timing used for PSD
  Float_t pp;              // Position of max amplitude 

} Can;

int NewScan (){
	
	
  /** ---------------------------------------------------- 
   *	Variable declairation
   *   ---------------------------------------------------- 
   */
  static Can det1, det2;
		
  bool	beamON,
    trg,
    data;
	
  float	X;
  int	multi;
	
  ifstream fp[16];
			 
  string 	line;

  /// Tracelength increased from 100 to 300 (1200us) for liquid cell tests
  int i,j, pposition,
    Tracelength = 300;
			
  float pulse[500],
    SG_pulse[500],
    SGderv2_pulse[500],
    baseline[500];
			
  Float_t amplitude,
    risetime,
    falltime,
    width,
    CFD,
    tac,
    paraL,
    paraS,
    runtime;
	
  // For SG filtered pulse
  Float_t trace_min, trace_max;
  int trace_min_loc, trace_max_loc;
  bool zero_cross;
  
  char 	filename[250],		
    prompt[10],
    openfile[250],
    prefix[100],
    interrputPrompt;
			
  Float_t trgtime, prevtime, difftime;
  Float_t prevtrgtime[10];
  long	TEvt = 0;

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
	15900,
	15900,
	15900,
	15900,
	15900,
	15900,
	15900,
	15900
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
    
  cout << "Run file directory name: ";
  cin >> prefix;

  TFile *ff = new TFile(filename,"RECREATE");

  TTree *tt = new TTree("T","Liquid Scintillator Tests");

  TH1F *trace0 = new TH1F("trace0","Trace for channel 0",300,0,1200);
  TH1F *trace1 = new TH1F("trace1","Trace for channel 1",300,0,1200);
  TH1F *trace2 = new TH1F("trace2","Trace for channel 2",300,0,1200);

  TH1F *trace0_SG = new TH1F("trace0_SG","SG filtered trace for channel 0",300,0,1200);
  TH1F *trace0_SG_derv2 = new TH1F("trace0_SG_derv2","2nd derivative of SG filtered trace for channel 0",300,0,1200);
  TH1F *trace1_SG = new TH1F("trace1_SG","SG filtered trace for channel 1",300,0,1200);
  TH1F *trace1_SG_derv2 = new TH1F("trace1_SG_derv2","2nd derivative of SG filtered trace for channel 1",300,0,1200);

	
  tt->Branch("det1",&det1,"l:s:amp:cfd:psd:trg:tac:pp");
  tt->Branch("det2",&det2,"l:s:amp:cfd:psd:trg:tac:pp");
  tt->Branch("runtime",&runtime,"runtime/F");     // Runtime in ms

  tt->Branch("trace0","TH1F",&trace0);
  tt->Branch("trace1","TH1F",&trace1);
  tt->Branch("trace2","TH1F",&trace2);

  tt->Branch("trace0_SG","TH1F",&trace0_SG);
  tt->Branch("trace0_SG_derv2","TH1F",&trace0_SG_derv2);
  tt->Branch("trace1_SG","TH1F",&trace1_SG);
  tt->Branch("trace1_SG_derv2","TH1F",&trace1_SG_derv2);
 
    
  const int numFiles = 2;	
  // Open files
  for (i = 0; i < numFiles; i++)
  {
    sprintf(openfile, "./%s/run0_wave%d.txt", prefix, i);
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
      beamON = 0;
      for (j = 0; j < numFiles; j++)
	{
	  if(!fp[j].is_open()){data = 0; cout << "Could not open file!!" << endl;}
	  if(fp[j].is_open())
	    {
	      trace_min = 0;
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
	      tac = 0;
	      pposition = -1;
	                  
	      // Get traces
	      for (i = 0; i < Tracelength; i++)
		{
		  if (!getline(fp[j], line)) {data = 0; break;}
		  pulse[i] = 16383 - atof(line.c_str());
		}

	      /** Liquid can processing **/
	      if(Tracelength > 1)
		{
		  // Process trace
		  Analysis->Baseline_restore(pulse, baseline, Tracelength, 5, 3); 
		  Analysis->Parameters(pulse, Tracelength, 3, &CFD, &amplitude, &risetime, &falltime, &width);
                

		  // Experimental routines
		  trace_min = 0;
		  trace_min_loc = 0;
		  for (i = 2; i< Tracelength - 2; i++) {
                    SG_pulse[i] = (-3.*pulse[i-2] + 12.*pulse[i-1] + 17.*pulse[i] + 12.*pulse[i+1] -3.*pulse[i+2])/35.;
                    SGderv2_pulse[i] = (2.*pulse[i-2] - pulse[i-1] - 2.*pulse[i] - pulse[i+1] + 2.*pulse[i+2])/7.;
                    trace0->SetBinContent(i, pulse[i]);
                    trace0_SG->SetBinContent(i, SG_pulse[i]);
                    trace0_SG_derv2->SetBinContent(i, SGderv2_pulse[i]);
                    
                    if(SGderv2_pulse[i] < trace_min) {
		      trace_min = SGderv2_pulse[i];
		      trace_min_loc = i;
                    }
		  }
		  
		  zero_cross = 0;
		  if (trace_min_loc > 1 && trace_min_loc < Tracelength - 2)  {
                    for (i = trace_min_loc; i>1; i--) {
		      if(SGderv2_pulse[i] > 0) {
			tac = -1.*(SGderv2_pulse[i] - (SGderv2_pulse[i] - SGderv2_pulse[i + 1])*(float)i)/(SGderv2_pulse[i] - SGderv2_pulse[i + 1]);
			break;
		      }
                    }
                    
                    for (i = trace_min_loc; i<Tracelength - 2; i++) {
		      if(SGderv2_pulse[i] > 0) {zero_cross = 1;}
		      if(SGderv2_pulse[i] < 0 && zero_cross) {
			tac = -1.*(SGderv2_pulse[i] - (SGderv2_pulse[i] - SGderv2_pulse[i - 1])*(float)i)/(SGderv2_pulse[i] - SGderv2_pulse[i - 1]) - tac;
			break;
		      }
                    }
		  }
		  		  
		  SG_pulse[0] = 0;
		  SG_pulse[1] = 0;
		  		  
		  Analysis->PSD_Integration(pulse, Tracelength, 12, 50, 10, 3, &paraL, &paraS);
		}
	      
	    }
	  

	  switch(j) {
	  case 0 : 
	    det1.amp = amplitude;
	    det1.l = paraL;
	    det1.s = paraS;
	    det1.cfd = CFD;
	    det1.trg = trg;
	    det1.tac = tac;
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
                    
	    break;

	  case 1 : 
	    det2.amp = amplitude;
	    det2.l = paraL;
	    det2.s = paraS;
	    det2.cfd = CFD;
	    det2.trg = trg;
	    det2.tac = tac;
	    det2.psd = det2.s / det2.l;
	    if (det2.trg){ det2.trg = 1;}
	    else {det2.trg = 0;}
	                     
	    break;
	  
	  }
                        
        }
      	
      tt->Fill();
      TEvt++;
      if (TEvt%1000==0) {cout << "\rEvent counter: " << TEvt << flush;}
	
    }
	
  
  for (i = 0; i < numFiles; i++)
    {
      if(fp[j].is_open()) fp[j].close();
    }
  ff->cd();
  tt->Write(); 
  ff->Close();
  
  cout << "\nFinsihed! " << TEvt << " events" << endl; 
    
  return 0;
}

int main() {
  return NewScan();
}
