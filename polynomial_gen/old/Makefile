SOPLEX_INCLUDE=$(SOPLEXPATH)/src/
SOPLEX_LIB=$(SOPLEXPATH)/build/lib/libsoplex.a

all: progressive

progressive: rlibm-fast.h progressive.cpp
	g++ -O2 -I$(SOPLEX_INCLUDE) -std=c++11 progressive.cpp  $(SOPLEX_LIB) -o polygen -lgmp -lz -lm

clean:
	rm -f *~ polygen filter *.txt
 
