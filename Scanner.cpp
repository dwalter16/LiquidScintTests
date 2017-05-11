/************************************************************************
 *
 *  Filename: Scanner.cpp
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
 *  To compile: g++ -O3 -pedantic -o Scan.exe `root-config --cflags --libs` -lSpectrum Scanner.cpp
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
    Float_t apl;             // Bar left after pulsing detected
    Float_t apr;             // Bar right after pulsing detected
    Float_t tdl;             // time since last trigger left (us)
    Float_t tdr;             // Time since last trigger right (us)
    Float_t ppl;             // Position of max amplitude left
    Float_t ppr;             // Position of max amplitude right
    Float_t npl;             // Number of peaks left
    Float_t npr;             // Number of peaks right
} Bar;  

typedef struct 
{
	Float_t l;               // Long integral
    Float_t s;               // Short integral
    Float_t a;               // Amplitude 
	Float_t cfd;             // Constant fraction timing
	Float_t psd;             // PSD parameter s/l
	Float_t tac;             // Zero-crossing timing used for PSD
} Can;


int Scanner (){
	
	
	/** ---------------------------------------------------- 
	*	Variable declairation
	*   ---------------------------------------------------- 
	*/
	Can     can_start,
            can_stop;
		
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
            SG_pulse[500],
            SGderv2_pulse[500],
            baseline[500];
			
	Float_t	amplitude,
			risetime,
			falltime,
			width,
			CFD,
            tac,
			paraL,
			paraS,
            runtime,
            FC,
            vll, vlr,
            val, var, hagA, hagL;
            
            
    // Temporary for testing
    Float_t trace_min, trace_max;
    int     trace_min_loc, trace_max_loc;
    bool    zero_cross;
			
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
	cout << " | Scanner.cpp                                   |" << endl; 
	cout << " |   Experiment: 13Can @ ORNL - MIRF             |" << endl; 
	cout << " |   Date: January 17th 2017                     |" << endl; 
	cout << " |   Calibration used: January 5th 2017          |" << endl; 
	cout << " |   ORNL Nuclear Astrophysics                   |" << endl;
	cout << " ------------------------------------------------ " << endl;
	
	cout << "Root file name to be created: ";
	cin >> filename;
    
    //cout << "Run file prefix: ";
    //cin >> prefix;

	TFile *ff = new TFile(filename,"RECREATE");

	TTree *tt = new TTree("T","n-ToF");
    TH1F *trace0 = new TH1F("trace0","Trace for channel 0",100,0,400);
    TH1F *trace1 = new TH1F("trace1","Trace for channel 1",100,0,400);
    TH1F *trace2 = new TH1F("trace2","Trace for channel 2",100,0,400);
    TH1F *trace3 = new TH1F("trace3","Trace for channel 3",100,0,400);
    TH1F *trace4 = new TH1F("trace4","Trace for channel 4",100,0,400);
    TH1F *trace5 = new TH1F("trace5","Trace for channel 5",100,0,400);
    TH1F *trace6 = new TH1F("trace6","Trace for channel 6",100,0,400);
    TH1F *trace7 = new TH1F("trace7","Trace for channel 7",100,0,400);
    TH1F *trace8 = new TH1F("trace8","Trace for channel 8",100,0,400);
    TH1F *trace9 = new TH1F("trace9","Trace for channel 9",100,0,400);
    
    TH1F *trace0_SG = new TH1F("trace0_SG","Trace for channel 0",100,0,400);
    TH1F *trace0_SG_derv2 = new TH1F("trace0_SG_derv2","Trace for channel 0",100,0,400);
    
	tt->Branch("can_start",&can_start,"l:s:a:cfd:psd:tac");
    tt->Branch("can_stop",&can_stop,"l:s:a:cfd:psd:tac");
    tt->Branch("runtime",&runtime,"runtime/F");     // Runtime in ms
    
    tt->Branch("trace0","TH1F",&trace0);
    tt->Branch("trace1","TH1F",&trace1);
    tt->Branch("trace2","TH1F",&trace2);
    tt->Branch("trace3","TH1F",&trace3);
    tt->Branch("trace4","TH1F",&trace4);
    tt->Branch("trace5","TH1F",&trace5);
    tt->Branch("trace6","TH1F",&trace6);
    tt->Branch("trace7","TH1F",&trace7);
    tt->Branch("trace8","TH1F",&trace8);
    tt->Branch("trace9","TH1F",&trace9);
    
    tt->Branch("trace0_SG","TH1F",&trace0_SG);
    tt->Branch("trace0_SG_derv2","TH1F",&trace0_SG_derv2);
    
	cout << "File to read: ";
	cin >> openfile;
    
    // Open files
    for (i = 0; i < 1; i++)
    {
        //sprintf(openfile, "%s_wave%d.txt
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
        {trace_min = 0;
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
            npeaks = 0;
            
            // Get trace
            for (i = 0; i < Tracelength; i++)
            {
                if (!getline(fp[j], line)) {data = 0; break;}
                pulse[i] = 16383 - atof(line.c_str());
                
            }

            /** Liquid can processing **/
            if(Tracelength > 1 && j < 1)
            {
                // Process trace
                Analysis->Baseline_restore(pulse, baseline, Tracelength, 10, 3); 
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
                Analysis->PSD_Integration(SG_pulse, Tracelength - 2, 5, 50, 10, 3, &paraL, &paraS);
                

            }
            
        }
            switch(j) {
                case 0 : 
                    can_stop.a = amplitude;
                    can_stop.l = paraL;
                    can_stop.s = paraS;
                    can_stop.psd = paraS/paraL;
                    can_stop.cfd = CFD;
                    can_stop.tac = tac;

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
        if(fp[i].is_open())
            fp[i].close();
    }
	ff->cd(); tt->Write(); 
	ff->Close();
	
    cout << "\nFinsihed! " << TEvt << " events" << endl; 
    
	return 0;
}

int main() {
    return Scanner();
}
