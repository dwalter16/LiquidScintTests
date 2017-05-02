# LiquidScintTests

Scanner.cpp was used for scanning the liquid bar data from MIRF 13C(a,n) experiment with the CAEN digitizer and program Wavedump.

NewScan.cpp is being modified from the existing scan code (Scanner.cpp). NewScan.cpp is to be used for time of flight characterizations of the new liquid cells. Handle (at first) two liquid cans for a start and stop signal.

To compile:

   g++ -03 -pedantic -o Scan.exe 'root-config --cflags --libs' -lSpectrum NewScan.cpp

- if errors compiling in Mac OSX
  - remove -03 option
  - clang++ instead of g++
  - $(root-config --clfags --libs) instead of 'root-config --cflags --libs

- if "error while loading shared libriaries: libcore.so: ..." occurs, type "source 'root-config --prefix'/bin/thisroot.sh"