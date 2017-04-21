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
  Float_t amp;              // Amplitude
  Float_t time;            // Trigger time 
  Float_t psd;             // PSD parameter s/l
  Float_t trg;             // Det trigger
  Float_t ap;             // Bar right after pulsing detected
  Float_t td;             // Time since last trigger  (us)
  Float_t pp;             // Position of max amplitude 
  Float_t np;             // Number of peaks 

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
	
  tt->Branch("can1",&can1,"l:s:amp::time:psd:trg::ap:td:pp:np");
  tt->Branch("can2",&can2,"l:s:amp::time:psd:trg::ap:td:pp:np");
  
  tt->Branch("runtime",&runtime,"runtime/F");     // Runtime in ms
    
  tt->Branch("vl",&vl,"vl/F");
  tt->Branch("val",&val,"val/F");
  tt->Branch("var",&var,"var/F");
	
  tt->Branch("hagA",&hagA,"hagA/F");
  tt->Branch("hagL",&hagL,"hagL/F");
  //cout << "File to read: ";
  //cin >> openfile;up and higher accuracy.
    
  // Open files
  for (i = 0; i < 1; i++)
    {
      //sprintf(openfile, "%s_wave%d.txt", prefix, i);
      sprintf(openfile, "BillsDet_1500V_cf%s.txt", prefix);
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

	      /** Liquid det processing **/
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
	    det1.al = amplitude;
	    det1.ll = paraL;
	    det1.sl = paraS;
	    det1.tl = CFD;
	    det1.trgl = trg;
	    det1.ppl = pposition;
		                    
	    det1.apl = ap;
   
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
	    det1.tdl = trgtime - prevtrgtime[0];
	    if(det1.tdl < 0) { det1.tdl = (8*(det1.tdl + 2147483647))/1.0E3; prevtrgtime[0] = trgtime;}
	    else { det1.tdl = (8*(det1.tdl))/1.0E3; prevtrgtime[0] = trgtime;}
                    
	    break;
	  case 1 : 
	    det1.ar = amplitude;
	    det1.lr = paraL;
	    det1.sr = paraS;
	    det1.tr = CFD;
	    det1.trgr = trg;
	    det1.ppr = pposition;
		    
	    det1.apr = ap;
                    
	    det1.l = sqrt(det1.ll * det1.lr)*cal[0];
	    det1.s = sqrt(det1.sl * det1.sr)*cal[0];
	    det1.psd = det1.s / det1.l;
	    det1.time = (4.0/2.0) * (det1.tl + det1.tr);
	    det1.position = (det1.ll - det1.lr)/(det1.ll + det1.lr)*19.56;
	    if (det1.trgl && det1.trgr)
	      {
		det1.trg = 1; 
		multi++;
	      }
	    else {det1.trg = 0;}
	    det1.dpsd = (det1.sl/det1.ll) - (det1.sr/det1.lr);
                    
	    //Time trigger difference
	    det1.tdr = trgtime - prevtrgtime[1];
	    if(det1.tdr < 0) { det1.tdr = (8*(det1.tdr + 2147483647))/1.0E3; prevtrgtime[1] = trgtime;}
	    else { det1.tdr = (8*(det1.tdr))/1.0E3; prevtrgtime[1] = trgtime;}

	    // Noise and after pulsing detection
	    if (det1.trg && det1.al<15500 && det1.ar<15500 && det1.psd > (egate[0][0]/sqrt(det1.l) + egate[0][1]*det1.l + egate[0][2]))
	      {
		det1.npl = s->Search(phistl, 2.0, "", 0.05);
		det1.npr = s->Search(phistr, 2.0, "", 0.05);
	      }
	    else { det1.npl = -1; det1.npr = -1; }

                    
	    break;
	  case 2 : 
	    det2.al = amplitude;
	    det2.ll = paraL;
	    det2.sl = paraS;
	    det2.tl = CFD;
	    det2.trgl = trg;
	    det2.ppl = pposition;
                    
	    det2.apl = ap;
                    
	    //Time trigger difference
	    det2.tdl = trgtime - prevtrgtime[2];
	    if(det2.tdl < 0) { det2.tdl = (8*(det2.tdl + 2147483647))/1.0E3; prevtrgtime[2] = trgtime;}
	    else { det2.tdl = (8*(det2.tdl))/1.0E3; prevtrgtime[2] = trgtime;}
                    
	    break;
	  case 3 : 
	    det2.ar = amplitude;
	    det2.lr = paraL;
	    det2.sr = paraS;
	    det2.tr = CFD;
	    det2.trgr = trg;
	    det2.ppr = pposition;
                    
	    det2.apr = ap;
                    
	    det2.l = sqrt(det2.ll * det2.lr)*cal[1];
	    det2.s = sqrt(det2.sl * det2.sr)*cal[1];
	    det2.psd = det2.s / det2.l;
	    det2.time = (4.0/2.0) * (det2.tl + det2.tr);
	    det2.position = (det2.ll - det2.lr)/(det2.ll + det2.lr)*19.56;
	    if (det2.trgl && det2.trgr)
	      {
		det2.trg = 1; 
		multi++;
	      }
	    else {det2.trg = 0;}
	    det3.dpsd = (det3.sl/det3.ll) - (det3.sr/det3.lr);
                    
	    //Time trigger difference
	    det2.tdr = trgtime - prevtrgtime[3];
	    if(det2.tdr < 0) { det2.tdr = (8*(det2.tdr + 2147483647))/1.0E3; prevtrgtime[3] = trgtime;}
	    else { det2.tdr = (8*(det2.tdr))/1.0E3; prevtrgtime[3] = trgtime;}

	    // Noise and after pulsing detection
	    if (det2.trg && det2.al<15500 && det2.ar<15500 && det2.psd > (egate[1][0]/sqrt(det2.l) + egate[1][1]*det2.l + egate[1][2]))
	      {
		det2.npl = s->Search(phistl, 2.0, "", 0.05);
		det2.npr = s->Search(phistr, 2.0, "", 0.05);
	      }
	    else { det2.npl = -1; det2.npr = -1; }
                    
	    break;
	  case 4 : 
	    det3.al = amplitude;
	    det3.ll = paraL;
	    det3.sl = paraS;
	    det3.tl = CFD;
	    det3.trgl = trg;
	    det3.ppl = pposition;
                    
	    det3.apl = ap;
                    
	    //Time trigger difference
	    det3.tdl = trgtime - prevtrgtime[4];
	    if(det3.tdl < 0) { det3.tdl = (8*(det3.tdl + 2147483647))/1.0E3; prevtrgtime[4] = trgtime;}
	    else { det3.tdl = (8*(det3.tdl))/1.0E3; prevtrgtime[4] = trgtime;}
                    
	    break;
	  case 5 : 
	    det3.ar = amplitude;
	    det3.lr = paraL;
	    det3.sr = paraS;
	    det3.tr = CFD;
	    det3.trgr = trg;
	    det3.ppr = pposition;

	    det3.apr = ap;
                    
	    det3.l = sqrt(det3.ll * det3.lr)*cal[2];
	    det3.s = sqrt(det3.sl * det3.sr)*cal[2];
	    det3.psd = det3.s / det3.l;
	    det3.time = (4.0/2.0) * (det3.tl + det3.tr);
	    det3.position = (det3.ll - det3.lr)/(det3.ll + det3.lr)*19.56;
	    if (det3.trgl && det3.trgr)
	      {
		det3.trg = 1; 
		multi++;
	      }
	    else {det3.trg = 0;}
	    det4.dpsd = (det4.sl/det4.ll) - (det4.sr/det4.lr);
                    
	    //Time trigger difference
	    det3.tdr = trgtime - prevtrgtime[5];
	    if(det3.tdr < 0) { det3.tdr = (8*(det3.tdr + 2147483647))/1.0E3; prevtrgtime[5] = trgtime;}
	    else { det3.tdr = (8*(det3.tdr))/1.0E3; prevtrgtime[5] = trgtime;}

	    // Noise and after pulsing detection
	    if (det3.trg && det3.al<15500 && det3.ar<15500 && det3.psd > (egate[2][0]/sqrt(det3.l) + egate[2][1]*det3.l + egate[2][2]))
	      {
		det3.npl = s->Search(phistl, 2.0, "", 0.05);
		det3.npr = s->Search(phistr, 2.0, "", 0.05);
	      }
	    else { det3.npl = -1; det3.npr = -1; }
                    
	    break;
	  case 6 : 
	    det4.al = amplitude;
	    det4.ll = paraL;
	    det4.sl = paraS;
	    det4.tl = CFD;
	    det4.trgl = trg;
	    det4.ppl = pposition;
                    
	    det4.apl = ap;
                    
	    //Time trigger difference
	    det4.tdl = trgtime - prevtrgtime[6];
	    if(det4.tdl < 0) { det4.tdl = (8*(det4.tdl + 2147483647))/1.0E3; prevtrgtime[6] = trgtime;}
	    else { det4.tdl = (8*(det4.tdl))/1.0E3; prevtrgtime[6] = trgtime;}
                    
	    break;
	  case 7 : 
	    det4.ar = amplitude;
	    det4.lr = paraL;
	    det4.sr = paraS;
	    det4.tr = CFD;
	    det4.trgr = trg;
	    det4.ppr = pposition;
		    
	    det4.apr = ap;
                    
	    det4.l = sqrt(det4.ll * det4.lr)*cal[3];
	    det4.s = sqrt(det4.sl * det4.sr)*cal[3];
	    det4.psd = det4.s / det4.l;
	    det4.time = (4.0/2.0) * (det4.tl + det4.tr);
	    det4.position = (det4.ll - det4.lr)/(det4.ll + det4.lr)*19.56;
	    if (det4.trgl && det4.trgr)
	      {
		det4.trg = 1; 
		multi++;
	      }
	    else {det4.trg = 0;}
                    
	    //Time trigger difference
	    det4.tdr = trgtime - prevtrgtime[7];
	    if(det4.tdr < 0) { det4.tdr = (8*(det4.tdr + 2147483647))/1.0E3; prevtrgtime[7] = trgtime;}
	    else { det4.tdr = (8*(det4.tdr))/1.0E3; prevtrgtime[7] = trgtime;}

	    // Noise and after pulsing detection
	    if (det4.trg && det4.al<15500 && det4.ar<15500 && det4.psd > (egate[3][0]/sqrt(det4.l) + egate[3][1]*det4.l + egate[3][2]))
	      {
		det4.npl = s->Search(phistl, 2.0, "", 0.05);
		det4.npr = s->Search(phistr, 2.0, "", 0.05);
	      }
	    else { det4.npl = -1; det4.npr = -1; }
                    
	    break;
	  case 8 : 
	    det5.al = amplitude;
	    det5.ll = paraL;
	    det5.sl = paraS;
	    det5.tl = CFD;
	    det5.trgl = trg;
	    det5.ppl = pposition;
                    
	    det5.apl = ap;
                    
	    //Time trigger difference
	    det5.tdl = trgtime - prevtrgtime[8];
	    if(det5.tdl < 0) { det5.tdl = (8*(det5.tdl + 2147483647))/1.0E3; prevtrgtime[8] = trgtime;}
	    else { det5.tdl = (8*(det5.tdl))/1.0E3; prevtrgtime[8] = trgtime;}
                    
	    break;
	  case 9 : 
	    det5.ar = amplitude;
	    det5.lr = paraL;
	    det5.sr = paraS;
	    det5.tr = CFD;
	    det5.trgr = trg;
	    det5.ppr = pposition;
                    
	    det5.apr = ap;
                    
	    det5.l = sqrt(det5.ll * det5.lr)*cal[4];
	    det5.s = sqrt(det5.sl * det5.sr)*cal[4];
	    det5.psd = det5.s / det5.l;
	    det5.time = (4.0/2.0) * (det5.tl + det5.tr);
	    det5.position = (det5.ll - det5.lr)/(det5.ll + det5.lr)*19.56;
	    if (det5.trgl && det5.trgr)
	      {
		det5.trg = 1; 
		multi++;
	      }
	    else {det5.trg = 0;}
	    det5.dpsd = (det5.sl/det5.ll) - (det5.sr/det5.lr);
                    
	    //Time trigger difference
	    det5.tdr = trgtime - prevtrgtime[9];
	    if(det5.tdr < 0) { det5.tdr = (8*(det5.tdr + 2147483647))/1.0E3; prevtrgtime[9] = trgtime;}
	    else { det5.tdr = (8*(det5.tdr))/1.0E3; prevtrgtime[9] = trgtime;}

	    // Noise and after pulsing detection
	    if (det5.trg && det5.al<15500 && det5.ar<15500 && det5.psd > (egate[4][0]/sqrt(det5.l) + egate[4][1]*det5.l + egate[4][2]))
	      {
		det5.npl = s->Search(phistl, 2.0, "", 0.05);
		det5.npr = s->Search(phistr, 2.0, "", 0.05);
	      }
	    else { det5.npl = -1; det5.npr = -1; }
                    
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
