CXX=g++
CXXFLAGS=-O2 -std=c++11
LIBS= -lmpfr 
TARGET1=flog2lib




all: $(TARGET1) 

.PHONY: $(TARGET1) 

$(TARGET1): flog2lib.cpp 
	$(CXX) $(CXXFLAGS) $<  $(LIBS) -o $@
	./$(TARGET1)
	


clean:
	rm -f $(TARGET1)

