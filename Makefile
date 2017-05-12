#For Mac OSX
CXX = clang++
CXXFLAGS = -pedantic 

#For linux
#CXX = g++
#CXXFLAGS = -03 -pedantic

ROOTFLAGS = $(shell root-config --cflags --libs) -lSpectrum

TARGET = NewScan.cpp

OBJS = NewScan.o

OBJ = newscan

$(OBJS): $(TARGET)
	$(CXX) $(CXXFLAGS) -o $(OBJ) $(ROOTFLAGS) $(TARGET) 

clean:
	-rm -f $(OBJ)