CXX=g++
CXXFLAGS=-O2 -I/home/soplex-4.0.1//src/ -std=c++11
LIBS= /home/soplex-4.0.1//build/lib/libsoplex.a -lgmp -lz -lm -lmpfr 
TARGET1=hlog2
TARGET2=hexpm1
TARGET3=flog2



all: $(TARGET1) $(TARGET2) $(TARGET3)

.PHONY: $(TARGET1) $(TARGET2) $(TARGET3)

$(TARGET3): flog2.cpp fhelper.cpp
	$(CXX) $(CXXFLAGS) $< fhelper.cpp $(LIBS) -o $@
	./$(TARGET3)


clean:
	rm -f $(TARGET1) $(TARGET2) $(TARGET3)

