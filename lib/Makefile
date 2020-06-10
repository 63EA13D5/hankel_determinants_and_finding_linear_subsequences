LIB_DIR=/usr/local/lib
NTL_INC=/usr/local/include/NTL
GMPL_INC=/usr/local/include

MAT_MAX_SIZE_OPTION=MAT_MAX_SIZE
MAT_MAX_SIZE=512

lib:
	g++ -Wall -Wpedantic -c -g -O5 -Iinclude src/detfct.cpp -std=c++11 -D${MAT_MAX_SIZE_OPTION}=${MAT_MAX_SIZE}
	ar rvs libdetfct.a detfct.o
	rm -rf detfct.o

lib_debug:
	g++ -Wall -Wpedantic -c -g -Iinclude src/detfct.cpp -std=c++11 -D${MAT_MAX_SIZE_OPTION}=${MAT_MAX_SIZE}
	ar rvs libdetfct.a detfct.o
	rm -rf detfct.o


test: lib
	g++ -Wall -Wpedantic -g -O5 -Iinclude -I${NTL_INC} -I${GMPL_INC} test/detfct_test.cpp -o detfct_test -L${LIB_DIR} -L. -ldetfct -lntl -lgmp -lm -lpthread -std=c++11 -D${MAT_MAX_SIZE_OPTION}=${MAT_MAX_SIZE}
	./detfct_test


clean:
	rm -rf libdetfct.a detfct_test
