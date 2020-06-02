LIB_DIR=/usr/local/lib
NTL_INC=/usr/local/include/NTL
GMPL_INC=/usr/local/include

LOCAL_LIB_DIR=lib

GF2: clean clean_lib
	cd lib && ${MAKE} lib
	g++ -Wall -Wpedantic -g -O5 -I${LOCAL_LIB_DIR}/include hanmat_GF2.cpp -o hanmat_GF2 -std=c++11 -L${LOCAL_LIB_DIR} -ldetfct

multithread: clean
	g++ -Wall -Wpedantic -g -O5 -I${NTL_INC} -I${GMPL_INC} hanmat_mt.cpp -o hanmat_mt -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread

hanmat: clean
	g++ -Wall -Wpedantic -g -O5 -I${NTL_INC} -I${GMPL_INC} hanmat.cpp -o hanmat -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread

hanmat_sampling: clean
	g++ -Wall -Wpedantic -g -O5 -I${NTL_INC} -I${GMPL_INC} hanmat_sampling.cpp -o hanmat_sampling -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread

debug: clean
	g++ -Wall -Wpedantic -g -I${NTL_INC} -I${GMPL_INC} hanmat_mt.cpp -o hanmat_mt -std=c++11 -L${LIB_DIR} -lntl -lgmp -lm -lpthread
	gdb hanmat_mt

clean:
	rm -rf hanmat hanmat_sampling hanmat_mt hanmat_GF2

clean_lib:
	cd lib && ${MAKE} clean
