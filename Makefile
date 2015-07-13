all:
	g++ -Wall -o sbfc -g -O2 sbfc.cpp -lgsl -lgslcblas
	
debug:	
	g++ -Wall -o sbfcd -g -DDEBUG sbfc.cpp -lgsl -lgslcblas
