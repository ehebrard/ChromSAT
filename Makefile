CPOHOME = /Users/ehebrard/Applications/IBM/ILOG/CPLEX_Studio1271
ARCHI = x86-64_osx

CCC = g++
CONCERTDIR = $(CPOHOME)/concert
CPLEXDIR = $(CPOHOME)/cplex
CFLAGS = -std=c++17 -DIL_STD -DDEBUG -g -I. -Isrc -I$(CPOHOME)/cpoptimizer/include -I$(CONCERTDIR)/include -I$(CPLEXDIR)/include -fPIC -pedantic -Wall -Wno-long-long -fexceptions -m64 -fno-strict-aliasing -DILOUSEMT -D_REENTRANT -DILM_REENTRANT
LDFLAGS = -L$(CPOHOME)/cpoptimizer/lib/$(ARCHI)/static_pic -lcp -L$(CPLEXDIR)/lib/$(ARCHI)/static_pic -lcplex -lilocplex -L$(CONCERTDIR)/lib/$(ARCHI)/static_pic -lconcert  -lpthread -framework CoreFoundation -framework IOKit
LIBRARYPATH = $(CPOHOME)/cpoptimizer/bin/$(ARCHI):$(CPLEXDIR)/bin/$(ARCHI)
CPPEXDIR = $(CPOHOME)/cpoptimizer/examples/src/cpp

cplex: src/cplex.cpp options.o intstack.o basic_graph.o
	$(CCC) -o cplex $(CFLAGS) src/cplex.cpp $(LDFLAGS) options.o basic_graph.o intstack.o 

options.o: src/options.cpp
	$(CCC) -o options.o -c $(CFLAGS) src/options.cpp
	
basic_graph.o: src/basic_graph.cpp 
	$(CCC) -o basic_graph.o -c $(CFLAGS) src/basic_graph.cpp 
	
intstack.o: src/intstack.cpp
	$(CCC) -o intstack.o -c $(CFLAGS) src/intstack.cpp

clean:
	rm -f cplex *.o
