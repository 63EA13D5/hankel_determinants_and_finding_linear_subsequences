# INSTALLATION

## Requirements
The application requires the following libraries to be installed :
- GNU Multiple Precision Arithmetic Library, GMP 6.2.0 ([link](https://gmplib.org/))
- Number Theory Library, NTL 11.4.3 ([link](https://www.shoup.net/ntl/))

## Installation of GMP library
IMPORTANT, GMP should be installed before NTL. The library has to be configured to be thread-safe, this is done by adding the flag `--enable-alloca=malloc-notreentrant` during configuration.
Requires root access:
```
wget https://gmplib.org/download/gmp/gmp-6.2.0.tar.lz
tar -xf gmp-6.2.0
cd gmp-6.2.0
./configure --enable-alloca=malloc-notreentrant
make
make check
sudo make install
```

## Installation of NTL library
NTL also needs to be configured to be thread-safe, this is done by adding the following flags `NTL_THREADS=on NTL_GMP_LIP=on` during configuration.
Requires root access:
```
wget https://www.shoup.net/ntl/ntl-11.4.3.tar.gz
tar -xf ntl-11.4.3.tar.gz
cd ntl-11.4.3.tar.gz
./configure NTL_THREADS=on NTL_GMP_LIP=on
make
make check
sudo make install
```

# Usage
The application can be built using the provided makefile:

## Build with NTL library
```
make multithread
./hanmat_mt
```

Build with NTL library for matrice representation and determinant computation. The data vector length and maximum of thread can be adjusted by directly modifying the `hanmat_mt.cpp` file.

## Build with custom GF2 library
```
make GF2 MAT_MAX_SIZE=1024
./hanmat_GF2
```

Build with vector of bitset for matrice representation and custom optimized determinant computation in GF2. Due to the static size allocation of `std::bitset` the data vector length has to be defined at compilation time using the compilation variable `MAT_MAX_SIZE`. The number of thread can be adjusted by directly modyfing the `hanmat_GF2.cpp` file.

# File of interest
- `hanmat_mt.cpp` - Features the NTL library based approach, single and multithreaded.
- `hanmat_GF2.cpp` - Features the custom GF2 library based approach, single and multithread.
- `lib/src/detfct.cpp` - Features the custom GF2 library, contains the optimized determinant function.
