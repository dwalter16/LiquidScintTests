/************************************************************************
 *
 *  Filename: NewScan.cpp
 *
 *  Description: 
 *
 *	Author(s):
 *     Michael T. Febbraro
 *     
 *
 *  Creation Date: 9/25/2016
 *  Last modified: 1/17/2017
 * 
 *  To compile: g++ -O3 -pedantic -o Scan.exe `root-config --cflags --libs` -lSpectrum NewScan.cpp
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
  Float_t ll;              // Long integral left (port)
  Float_t lr;              // Long integral right
  Float_t s;               // Short integral sqrt(sl*sr)
  Float_t sl;              // Short integral left (port)
  Float_t sr;              // Short integral right
  Float_t al;              // Amplitude left (port)
  Float_t ar;              // Amplitude right 
  Float_t time;            // Trigger time 0.5 * (tl + tr)
  Float_t tl;              // Trigger time left (port)
  Float_t tr;              // Trigger time right
  Float_t position;        // Position ln(ll/lr)
  Float_t psd;             // PSD parameter s/l
  Float_t trgl;            // Left side trigger
  Float_t trgr;            // Right side triger
  Float_t trg;             // Bar left and right side trigger
  Float_t dpsd;            // Bar left psd - bar right psd
  Float_t apl;             // Bar left after pulsing detected
  Float_t apr;             // Bar right after pulsing detected
  Float_t tdl;             // time since last trigger left (us)
  Float_t tdr;             // Time since last trigger right (us)
  Float_t ppl;             // Position of max amplitude left
  Float_t ppr;             // Position of max amplitude right
  Float_t npl;             // Number of peaks left
  Float_t npr;             // Number of peaks right
} Bar;

int NewScan (){
	
	
  /** ---------------------------------------------------- 
   *	Variable declairation
   *   ---------------------------------------------------- 
   */
  static Bar bar1,
    bar2,
    bar3,
    bar4,
    bar5;
		
  bool	beamON,
    trg,
    data;
	
  float	X, Y, beamSweeper;
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
    runtime,
    FC,
    vll, vlr,
    val, var, hagA, hagL;
			
  char 	filename[250],		
    prompt[10],
    openfile[250],
    prefix[100],
    interrputPrompt;
			
  Float_t trgtime, prevtime, difftime;
  Float_t prevtrgtime[10];
  long	TEvt = 0;

  TH1F *phistl = new TH1F("hl","hl",Tracelength,0,Tracelength-1);
  TH1F *phistr = new TH1F("hr","hr",Tracelength,0,Tracelength-1);

  TSpectrum *s = new TSpectrum();
	
  /** ---------------------------------------------------- 
   *	Calibrations and threshold
   *   ---------------------------------------------------- 
   */
        
  float cal[5] = 
    {   0.0833,
	0.0815,
	0.0753,
	0.1574,
	0.1227
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
    
  cout << "Run file prefix: ";
  cin >> prefix;

  TFile *ff = new TFile(filename,"RECREATE");

  TTree *tt = new TTree("T","Liquid Scintillator Tests");
	
  tt->Branch("bar1",&bar1,"l:ll:lr:s:sl:sr:al:ar:time:tl:tr:position:psd:trgl:trgr:trg:dpsd:apl:apr:tdl:tdr:ppl:ppr:npl:npr");
  tt->Branch("bar2",&bar2,"l:ll:lr:s:sl:sr:al:ar:time:tl:tr:position:psd:trgl:trgr:trg:dpsd:apl:apr:tdl:tdr:ppl:ppr:npl:npr");
  tt->Branch("bar3",&bar3,"l:ll:lr:s:sl:sr:al:ar:time:tl:tr:position:psd:trgl:trgr:trg:dpsd:apl:apr:tdl:tdr:ppl:ppr:npl:npr");
  tt->Branch("bar4",&bar4,"l:ll:lr:s:sl:sr:al:ar:time:tl:tr:position:psd:trgl:trgr:trg:dpsd:apl:apr:tdl:tdr:ppl:ppr:npl:npr");
  tt->Branch("bar5",&bar5,"l:ll:lr:s:sl:sr:al:ar:time:tl:tr:position:psd:trgl:trgr:trg:dpsd:apl:apr:tdl:tdr:ppl:ppr:npl:npr");
  tt->Branch("FC",&FC,"FC/F");
  tt->Branch("y",&Y,"y/F");
  tt->Branch("m",&multi,"m/I");
  tt->Branch("beamON",&beamON,"beamON/O");
  tt->Branch("beamSweeper",&beamSweeper,"beamSweeper/F");
  tt->Branch("runtime",&runtime,"runtime/F");     // Runtime in ms
    
  tt->Branch("vll",&vll,"vll/F");
  tt->Branch("vlr",&vlr,"vlr/F");
  tt->Branch("val",&val,"val/F");
  tt->Branch("var",&var,"var/F");
	
  tt->Branch("hagA",&hagA,"hagA/F");
  tt->Branch("hagL",&hagL,"hagL/F");
  //cout << "File to read: ";
  //cin >> openfile;up and higher accuracy.
    
  // Open files
  for (i = 0; i < 1; i++)
    {
      sprintf(openfile, "%s_wave%d.txt", prefix, i);
      fp[i].open(openfile, std::ifstream::in);
    }
    
  data = 1;
  runtime = 0;
  prevtime = 0;
    
  /** ---------------------------------------------------- 
   *	Process liquid scintillator bar events
   *   ----------------------------------------------------  
   */
	
  while (data)
    {
      multi = 0;
      beamSweeper = -1;
      X = -1;
      Y = -1;
      ap = 0;
      beamON = 0;
      for (j = 0; j < 1; j++)
	{

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
            
	      // Get trace
	      for (i = 0; i < Tracelength; i++)
		{
		  if (!getline(fp[j], line)) {data = 0; break;}
		  pulse[i] = 16383 - atof(line.c_str());

		  switch(j) {
		  case 0 : phistl->SetBinContent(i, pulse[i]); break;
		  case 1 : phistr->SetBinContent(i, pulse[i]); break;
		  case 2 : phistl->SetBinContent(i, pulse[i]); break;
		  case 3 : phistr->SetBinContent(i, pulse[i]); break;
		  case 4 : phistl->SetBinContent(i, pulse[i]); break;
		  case 5 : phistr->SetBinContent(i, pulse[i]); break;
		  case 6 : phistl->SetBinContent(i, pulse[i]); break;
		  case 7 : phistr->SetBinContent(i, pulse[i]); break;
		  case 8 : phistl->SetBinContent(i, pulse[i]); break;
		  case 9 : phistr->SetBinContent(i, pulse[i]); break;
		  }
                
		  if (j < 10)
                    if((16838 - pulse[i]) <= threshold[j]) {trg = 1;}
		}

	      /** Liquid bar processing **/
	      if(Tracelength > 1 && j < 10)
		{
		  // Process trace
		  Analysis->Baseline_restore(pulse, baseline, Tracelength, 10, 3); 
		  Analysis->Parameters(pulse, Tracelength, 3, &CFD, &amplitude, &risetime, &falltime, &width);
                
		  // Note this method has been updated to including after pulsing detection
		  Analysis->PSD_Integration_Afterpulsing(pulse, Tracelength, 12, 50, 5, 15, 50, &pposition, &paraL, &paraS, &ap);

		}
            
	      /** Veto panel processing **/
	      if(Tracelength > 1 && j > 12)
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
	    bar1.al = amplitude;
	    bar1.ll = paraL;
	    bar1.sl = paraS;
	    bar1.tl = CFD;
	    bar1.trgl = trg;
	    bar1.ppl = pposition;
		                    
	    bar1.apl = ap;
   
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
	    bar1.tdl = trgtime - prevtrgtime[0];
	    if(bar1.tdl < 0) { bar1.tdl = (8*(bar1.tdl + 2147483647))/1.0E3; prevtrgtime[0] = trgtime;}
	    else { bar1.tdl = (8*(bar1.tdl))/1.0E3; prevtrgtime[0] = trgtime;}
                    
	    break;
	  case 1 : 
	    bar1.ar = amplitude;
	    bar1.lr = paraL;
	    bar1.sr = paraS;
	    bar1.tr = CFD;
	    bar1.trgr = trg;
	    bar1.ppr = pposition;
		    
	    bar1.apr = ap;
                    
	    bar1.l = sqrt(bar1.ll * bar1.lr)*cal[0];
	    bar1.s = sqrt(bar1.sl * bar1.sr)*cal[0];
	    bar1.psd = bar1.s / bar1.l;
	    bar1.time = (4.0/2.0) * (bar1.tl + bar1.tr);
	    bar1.position = (bar1.ll - bar1.lr)/(bar1.ll + bar1.lr)*19.56;
	    if (bar1.trgl && bar1.trgr)
	      {
		bar1.trg = 1; 
		multi++;
	      }
	    else {bar1.trg = 0;}
	    bar1.dpsd = (bar1.sl/bar1.ll) - (bar1.sr/bar1.lr);
                    
	    //Time trigger difference
	    bar1.tdr = trgtime - prevtrgtime[1];
	    if(bar1.tdr < 0) { bar1.tdr = (8*(bar1.tdr + 2147483647))/1.0E3; prevtrgtime[1] = trgtime;}
	    else { bar1.tdr = (8*(bar1.tdr))/1.0E3; prevtrgtime[1] = trgtime;}

	    // Noise and after pulsing detection
	    if (bar1.trg && bar1.al<15500 && bar1.ar<15500 && bar1.psd > (egate[0][0]/sqrt(bar1.l) + egate[0][1]*bar1.l + egate[0][2]))
	      {
		bar1.npl = s->Search(phistl, 2.0, "", 0.05);
		bar1.npr = s->Search(phistr, 2.0, "", 0.05);
	      }
	    else { bar1.npl = -1; bar1.npr = -1; }

                    
	    break;
	  case 2 : 
	    bar2.al = amplitude;
	    bar2.ll = paraL;
	    bar2.sl = paraS;
	    bar2.tl = CFD;
	    bar2.trgl = trg;
	    bar2.ppl = pposition;
                    
	    bar2.apl = ap;
                    
	    //Time trigger difference
	    bar2.tdl = trgtime - prevtrgtime[2];
	    if(bar2.tdl < 0) { bar2.tdl = (8*(bar2.tdl + 2147483647))/1.0E3; prevtrgtime[2] = trgtime;}
	    else { bar2.tdl = (8*(bar2.tdl))/1.0E3; prevtrgtime[2] = trgtime;}
                    
	    break;
	  case 3 : 
	    bar2.ar = amplitude;
	    bar2.lr = paraL;
	    bar2.sr = paraS;
	    bar2.tr = CFD;
	    bar2.trgr = trg;
	    bar2.ppr = pposition;
                    
	    bar2.apr = ap;
                    
	    bar2.l = sqrt(bar2.ll * bar2.lr)*cal[1];
	    bar2.s = sqrt(bar2.sl * bar2.sr)*cal[1];
	    bar2.psd = bar2.s / bar2.l;
	    bar2.time = (4.0/2.0) * (bar2.tl + bar2.tr);
	    bar2.position = (bar2.ll - bar2.lr)/(bar2.ll + bar2.lr)*19.56;
	    if (bar2.trgl && bar2.trgr)
	      {
		bar2.trg = 1; 
		multi++;
	      }
	    else {bar2.trg = 0;}
	    bar3.dpsd = (bar3.sl/bar3.ll) - (bar3.sr/bar3.lr);
                    
	    //Time trigger difference
	    bar2.tdr = trgtime - prevtrgtime[3];
	    if(bar2.tdr < 0) { bar2.tdr = (8*(bar2.tdr + 2147483647))/1.0E3; prevtrgtime[3] = trgtime;}
	    else { bar2.tdr = (8*(bar2.tdr))/1.0E3; prevtrgtime[3] = trgtime;}

	    // Noise and after pulsing detection
	    if (bar2.trg && bar2.al<15500 && bar2.ar<15500 && bar2.psd > (egate[1][0]/sqrt(bar2.l) + egate[1][1]*bar2.l + egate[1][2]))
	      {
		bar2.npl = s->Search(phistl, 2.0, "", 0.05);
		bar2.npr = s->Search(phistr, 2.0, "", 0.05);
	      }
	    else { bar2.npl = -1; bar2.npr = -1; }
                    
	    break;
	  case 4 : 
	    bar3.al = amplitude;
	    bar3.ll = paraL;
	    bar3.sl = paraS;
	    bar3.tl = CFD;
	    bar3.trgl = trg;
	    bar3.ppl = pposition;
                    
	    bar3.apl = ap;
                    
	    //Time trigger difference
	    bar3.tdl = trgtime - prevtrgtime[4];
	    if(bar3.tdl < 0) { bar3.tdl = (8*(bar3.tdl + 2147483647))/1.0E3; prevtrgtime[4] = trgtime;}
	    else { bar3.tdl = (8*(bar3.tdl))/1.0E3; prevtrgtime[4] = trgtime;}
                    
	    break;
	  case 5 : 
	    bar3.ar = amplitude;
	    bar3.lr = paraL;
	    bar3.sr = paraS;
	    bar3.tr = CFD;
	    bar3.trgr = trg;
	    bar3.ppr = pposition;

	    bar3.apr = ap;
                    
	    bar3.l = sqrt(bar3.ll * bar3.lr)*cal[2];
	    bar3.s = sqrt(bar3.sl * bar3.sr)*cal[2];
	    bar3.psd = bar3.s / bar3.l;
	    bar3.time = (4.0/2.0) * (bar3.tl + bar3.tr);
	    bar3.position = (bar3.ll - bar3.lr)/(bar3.ll + bar3.lr)*19.56;
	    if (bar3.trgl && bar3.trgr)
	      {
		bar3.trg = 1; 
		multi++;
	      }
	    else {bar3.trg = 0;}
	    bar4.dpsd = (bar4.sl/bar4.ll) - (bar4.sr/bar4.lr);
                    
	    //Time trigger difference
	    bar3.tdr = trgtime - prevtrgtime[5];
	    if(bar3.tdr < 0) { bar3.tdr = (8*(bar3.tdr + 2147483647))/1.0E3; prevtrgtime[5] = trgtime;}
	    else { bar3.tdr = (8*(bar3.tdr))/1.0E3; prevtrgtime[5] = trgtime;}

	    // Noise and after pulsing detection
	    if (bar3.trg && bar3.al<15500 && bar3.ar<15500 && bar3.psd > (egate[2][0]/sqrt(bar3.l) + egate[2][1]*bar3.l + egate[2][2]))
	      {
		bar3.npl = s->Search(phistl, 2.0, "", 0.05);
		bar3.npr = s->Search(phistr, 2.0, "", 0.05);
	      }
	    else { bar3.npl = -1; bar3.npr = -1; }
                    
	    break;
	  case 6 : 
	    bar4.al = amplitude;
	    bar4.ll = paraL;
	    bar4.sl = paraS;
	    bar4.tl = CFD;
	    bar4.trgl = trg;
	    bar4.ppl = pposition;
                    
	    bar4.apl = ap;
                    
	    //Time trigger difference
	    bar4.tdl = trgtime - prevtrgtime[6];
	    if(bar4.tdl < 0) { bar4.tdl = (8*(bar4.tdl + 2147483647))/1.0E3; prevtrgtime[6] = trgtime;}
	    else { bar4.tdl = (8*(bar4.tdl))/1.0E3; prevtrgtime[6] = trgtime;}
                    
	    break;
	  case 7 : 
	    bar4.ar = amplitude;
	    bar4.lr = paraL;
	    bar4.sr = paraS;
	    bar4.tr = CFD;
	    bar4.trgr = trg;
	    bar4.ppr = pposition;
		    
	    bar4.apr = ap;
                    
	    bar4.l = sqrt(bar4.ll * bar4.lr)*cal[3];
	    bar4.s = sqrt(bar4.sl * bar4.sr)*cal[3];
	    bar4.psd = bar4.s / bar4.l;
	    bar4.time = (4.0/2.0) * (bar4.tl + bar4.tr);
	    bar4.position = (bar4.ll - bar4.lr)/(bar4.ll + bar4.lr)*19.56;
	    if (bar4.trgl && bar4.trgr)
	      {
		bar4.trg = 1; 
		multi++;
	      }
	    else {bar4.trg = 0;}
                    
	    //Time trigger difference
	    bar4.tdr = trgtime - prevtrgtime[7];
	    if(bar4.tdr < 0) { bar4.tdr = (8*(bar4.tdr + 2147483647))/1.0E3; prevtrgtime[7] = trgtime;}
	    else { bar4.tdr = (8*(bar4.tdr))/1.0E3; prevtrgtime[7] = trgtime;}

	    // Noise and after pulsing detection
	    if (bar4.trg && bar4.al<15500 && bar4.ar<15500 && bar4.psd > (egate[3][0]/sqrt(bar4.l) + egate[3][1]*bar4.l + egate[3][2]))
	      {
		bar4.npl = s->Search(phistl, 2.0, "", 0.05);
		bar4.npr = s->Search(phistr, 2.0, "", 0.05);
	      }
	    else { bar4.npl = -1; bar4.npr = -1; }
                    
	    break;
	  case 8 : 
	    bar5.al = amplitude;
	    bar5.ll = paraL;
	    bar5.sl = paraS;
	    bar5.tl = CFD;
	    bar5.trgl = trg;
	    bar5.ppl = pposition;
                    
	    bar5.apl = ap;
                    
	    //Time trigger difference
	    bar5.tdl = trgtime - prevtrgtime[8];
	    if(bar5.tdl < 0) { bar5.tdl = (8*(bar5.tdl + 2147483647))/1.0E3; prevtrgtime[8] = trgtime;}
	    else { bar5.tdl = (8*(bar5.tdl))/1.0E3; prevtrgtime[8] = trgtime;}
                    
	    break;
	  case 9 : 
	    bar5.ar = amplitude;
	    bar5.lr = paraL;
	    bar5.sr = paraS;
	    bar5.tr = CFD;
	    bar5.trgr = trg;
	    bar5.ppr = pposition;
                    
	    bar5.apr = ap;
                    
	    bar5.l = sqrt(bar5.ll * bar5.lr)*cal[4];
	    bar5.s = sqrt(bar5.sl * bar5.sr)*cal[4];
	    bar5.psd = bar5.s / bar5.l;
	    bar5.time = (4.0/2.0) * (bar5.tl + bar5.tr);
	    bar5.position = (bar5.ll - bar5.lr)/(bar5.ll + bar5.lr)*19.56;
	    if (bar5.trgl && bar5.trgr)
	      {
		bar5.trg = 1; 
		multi++;
	      }
	    else {bar5.trg = 0;}
	    bar5.dpsd = (bar5.sl/bar5.ll) - (bar5.sr/bar5.lr);
                    
	    //Time trigger difference
	    bar5.tdr = trgtime - prevtrgtime[9];
	    if(bar5.tdr < 0) { bar5.tdr = (8*(bar5.tdr + 2147483647))/1.0E3; prevtrgtime[9] = trgtime;}
	    else { bar5.tdr = (8*(bar5.tdr))/1.0E3; prevtrgtime[9] = trgtime;}

	    // Noise and after pulsing detection
	    if (bar5.trg && bar5.al<15500 && bar5.ar<15500 && bar5.psd > (egate[4][0]/sqrt(bar5.l) + egate[4][1]*bar5.l + egate[4][2]))
	      {
		bar5.npl = s->Search(phistl, 2.0, "", 0.05);
		bar5.npr = s->Search(phistr, 2.0, "", 0.05);
	      }
	    else { bar5.npl = -1; bar5.npr = -1; }
                    
	    break;
	  case 10 : 
	    beamSweeper = pulse[0];
	    break;
	  case 11 :
	    FC = pulse[0];
	    break;
	  case 12 :
	    Y = pulse[0];
	    break;
	  case 13 :
	    vll = paraL;
	    val = amplitude;
	    break;
	  case 14 :
	    vlr = paraL;
	    var = amplitude;
	    break;
	  case 15 :
	    hagA = amplitude;
	    hagL = paraL;
	    break;            }
               
            
            
        }
    
      //if (TEvt > 100) {data = 0;}
	
      tt->Fill();
      TEvt++;
      if (TEvt%1000==0) {cout << "\rEvent counter: " << TEvt << flush;}
	
    }
	
	
  for (i = 0; i < 15; i++)
    {
      if(fp[i].is_open())
	fp[i].close();
    }
  ff->cd(); tt->Write(); 
  ff->Close();
	
  cout << "\nFinsihed! " << TEvt << " events" << endl; 
    
  return 0;
}

int main() {
  return NewScan();
}
