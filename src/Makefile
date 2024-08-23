LDFLAGS="-L./lib -lpthread -lgzstream -lz -Wl,-rpath ${PWD}/lib"

default: count_tnf fakmer

count_tnf : count_tnf.cpp
	g++ -O3 -o ../bin/count_tnf count_tnf.cpp -I./include ${LDFLAGS}

fakmer: fakmer.cpp
	g++ -O3 -o ../bin/fakmer fakmer.cpp -I./include ${LDFLAGS}

faKmtx: faKmtx.cpp
	g++ -O3 -o ../bin/faKmtx faKmtx.cpp -I./include ${LDFLAGS}
	
fqKdist: fqKdist.cpp
	g++ -O3 -o ../bin/fqKdist fqKdist.cpp -I./include ${LDFLAGS}

clean:
rm -f ../bin/count_tnf
