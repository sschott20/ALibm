CXX=g++
CXXFLAGS=-O2 -I/home/soplex-4.0.1//src/ -std=c++11
LIBS= /home/soplex-4.0.1//build/lib/libsoplex.a -lgmp -lz -lm -lmpfr 
TARGET1=hlog2
TARGET2=hexpm1
TARGET3=float


all: $(TARGET1) $(TARGET2)

.PHONY: $(TARGET1) $(TARGET2) $(TARGET3)


$(TARGET1): halflog2.cpp 
	$(CXX) $(CXXFLAGS) $< $(LIBS) -o $@
	./$(TARGET1)

$(TARGET2): halfexpm1.cpp 
	$(CXX) $(CXXFLAGS) $< $(LIBS) -o $@
	./$(TARGET2)

$(TARGET3): float.cpp
	$(CXX) $(CXXFLAGS) $< $(LIBS) -o $@
	./$(TARGET2)


clean:
	rm -f $(TARGET1) $(TARGET2)

