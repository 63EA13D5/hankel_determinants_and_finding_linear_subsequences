#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <bitset>

#include <NTL/matrix.h>
#include <NTL/GF2.h>
#include <NTL/mat_GF2.h>

#include "detfct.h"

void test(long N)
{
	std::random_device r_dev;
	std::mt19937_64 mt(r_dev());
	std::bernoulli_distribution dist(0.5);

	double time_ntl=0.0;
	double time_artisanal=0.0;

	std::vector<std::bitset<MAT_MAX_SIZE>> A;
	for(long r = 0 ; r < N; r++)
	{
		std::bitset<MAT_MAX_SIZE> tmp;
		A.push_back(tmp);
	}

	NTL::Mat<NTL::GF2> M;
	M.SetDims(N,N);

	int chk_sum = 0;
	int nb_tests = 1000;
	for(int test_ct = 0 ; test_ct < nb_tests; test_ct++)
	{
		for(int r = 0 ; r < N ; r++)
		{
			for(int c = 0 ; c < N ; c++)
			{
				if( (bool)dist(mt) )
				{
					A[r].set(c);
				}
				M[r][c] = NTL::GF2(A[r][c]);
			}
		}

		auto t_start = std::chrono::high_resolution_clock::now();
		NTL::GF2 detNTL = NTL::determinant(M);
		auto t_end = std::chrono::high_resolution_clock::now();
		time_ntl += std::chrono::duration<double, std::milli>(t_end-t_start).count();

		t_start = std::chrono::high_resolution_clock::now();
		bool det_artisanal = GF2_Utils::det_b(A, N);
		t_end = std::chrono::high_resolution_clock::now();
		time_artisanal += std::chrono::duration<double, std::milli>(t_end-t_start).count();

		if(detNTL==NTL::GF2(det_artisanal))
		{
			chk_sum++;
		}
	}

	std::cout << "\n\n" << chk_sum << " / " << nb_tests << "\n\n";
	std::cout << "NTL : " << time_ntl/(double)nb_tests;
	std::cout << "\nartisanal : " << time_artisanal/(double)nb_tests << "\n\n";
}

int main(void)
{
	std::cout << "TEST FOR N=512, MAT_MAX_SIZE=" << MAT_MAX_SIZE << "\n";
	test(512);

	std::cout << "TEST FOR N=256, MAT_MAX_SIZE=" << MAT_MAX_SIZE << "\n";
	test(256);

	return 0;
}
