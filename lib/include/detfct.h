/*
 * Library to compute determinant for GF2 matrices
 *
 * After importing the library please redefine N with the used size
 */

#ifndef DETFCT_H
#define DETFCT_H

#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <bitset>

#ifndef MAT_MAX_SIZE
#define MAT_MAX_SIZE 512
#endif

namespace GF2_Utils {
	/*
	   Triangularise M sur place donc M originale detruite.
	   Determinant est egale aux produits des elements de la diagonale de la triangularisation.
	   Des qu'un element nul est rencontre sur la diagonale alors retourne zero.
	   */
	bool det_b(std::vector<std::bitset<MAT_MAX_SIZE>> M, size_t msize);

	bool chk_triangular_tables_not_equal(std::vector<std::bitset<MAT_MAX_SIZE>> & M1, std::vector<std::bitset<MAT_MAX_SIZE>> & M2, long maxdim, size_t msize);

	int how_many_differences(std::vector<std::bitset<MAT_MAX_SIZE>> & M1, std::vector<std::bitset<MAT_MAX_SIZE>> & M2, long maxdim, size_t msize);
}

#endif
