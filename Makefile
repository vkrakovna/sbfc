all:
	g++ -Wall -o sbfc -g -O2 sbfc.cpp -lgsl -lgslcblas
	
debug:	
	g++ -Wall -o sbfcd -g -DDEBUG sbfc.cpp -lgsl -lgslcblas
	
test:
	g++ -Wall -o sbfcu -g -DTEST unit_tests.cpp -lgsl -lgslcblas
