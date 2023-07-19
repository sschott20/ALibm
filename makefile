float: FORCE
	g++ -O0 -Wall --std=c++14 -L/usr/lib -o float float.cpp -lmpfr -lgmp ; ./float

comp: FORCE
	g++ -O0 -Wall --std=c++14 -L/usr/lib -o comp compute_correct.cpp -lmpfr -lgmp ; ./comp
test: FORCE
	g++ -O0 -Wall --std=c++14 -L/usr/lib -o test test.cpp -lmpfr -lgmp ; ./test	
FORCE: ;