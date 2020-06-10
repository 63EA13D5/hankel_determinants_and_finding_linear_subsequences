/***** BEGINNING OF DISCLAIMER *****

  This code is provided without warranty of any kind, either express or implied. Use at your own risk.

  The use of this is done at your own discretion and risk and with agreement that you will be solely responsible for any damage to your computer system or loss of data that results from the use of this code. You are solely responsible for adequate protection and backup of the data and equipment used in connection with this code, and we will not be liable for any damages that you may suffer in connection with using, modifying or distributing any of this code. No advice or information, whether oral or written, obtained by you from us or from this website shall create any warranty for this code.

  We make makes no warranty that:

  1) the code will meet your requirements
  2) the code will be uninterrupted, timely, secure or error-free
  3) the results that may be obtained from the use of the software will be effective, accurate or reliable
  4) the quality of the code will meet your expectations
  5) any errors in the code obtained from us will be corrected. 

  The code and its documentation made available on this website:

  1) could include technical or other mistakes, inaccuracies or typographical errors. We may make changes to the code or documentation made available on its web site at any time without prior-notice.
  2) may be out of date, and we make no commitment to update such materials. 

  We assume no responsibility for errors or omissions in the code or documentation available from its web site.

  In no event shall we be liable to you or any third parties for any special, punitive, incidental, indirect or consequential damages of any kind, or any damages whatsoever, including, without limitation, those resulting from loss of use, data or profits, and on any theory of liability, arising out of or in connection with the use of this code. 

 ***** END OF DISCLAIMER *****/

/*

AUTHORS: Claude Gravel and Bastien Rigault

DATE: April 7th, 2020 

TECHNICAL REFERENCE: (arxiv) "Finding linearly generated subsequences" by Claude Gravel and Daniel Panario

NOTES:

This code relies on NTL library which may (or may not) depend on GMPL.
For information about NTL (Number Theory Library) see: https://www.shoup.net/ntl/
For information about GMPL (GNU Multiple Precision Arithmetic Library) see: https://gmplib.org/

*/

#include <iomanip>
#include <iostream>
#include <fstream>
#include <random>
#include <ctime>
#include <chrono>
#include <thread>
#include <mutex>
#include <vector>
#include <bitset>

#include "detfct.h"


/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

const int N=MAT_MAX_SIZE; //length data vector
const int nb_trials=100;//sample size for timing (and debugging)
const int max_dim=(N/2)*(N%2 == 0) + ((N+1)/2)*(N%2 == 1);
const double ratio_size=1.0/16.0;
const int len_C_gen = (int)floor(N*ratio_size);//length of generating vector
const int number_of_threads = 4; // Maximum number of thread that can be created for each j-th row

/*
   Possible locations of the hidden linear subsequence.
   Case [1] is representative of a hard case.
   Case [2] is representative of an easy case.
   */

const int left_index = (int)(floor(7*len_C_gen));//where the linear subsequence begins, case [1]
const int right_index = (int)(floor(9*len_C_gen));//where the linear subsequence ends, case [1]
//const int left_index = len_C_gen;//where the linear subsequence begins, case [2]
//const int right_index = n;//where the linear subsequence ends, case [2]

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

const int max_len_side_grid=7;//Research required to increase value here...

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

const bool RUN_NAIVE_FLAG = true;
//If RUN_NAIVE_FLAG is false, then skip naive (thus resulting in nan, warning of mismatch between naive and new methods unless the last part is modified appropriately.)
//RUN_NAIVE_FLAG can be false when the new method is wanted solely and no need to compare methods.

bool MULTI_THREAD_ON=false;

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

std::bitset<len_C_gen> C_gen; //generation vector
std::bitset<MAT_MAX_SIZE> V0; //data vector to be filled with a linear subsequence

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

std::mutex determinant_mutex;
std::mutex asptl_mutex;

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

struct experiment {
	std::vector<std::bitset<MAT_MAX_SIZE>> M;
	std::vector<std::bitset<MAT_MAX_SIZE>> flags_M; // Not used by trivial approach
	std::vector<std::bitset<MAT_MAX_SIZE>> AspTL[max_dim+1]; // Not used by trivial approach
};

struct determinant_result {
	int j, i;
	bool determinant;
}; //Data structure used to make the multithread version of our fast algorithm outputs correct result

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

/* PRINT FUNCTIONS */

void print_stats(double time_triv_st, double time_triv_mt, double time_fast_st, double time_fast_mt)// for timing
{
	/*
	   To output timing results for a single case.
	   */

	int number_of_entries = ((N%2)==0)*((N/2)-1)*(N/2)+((N%2)==1)*(N/2)*(N-1)/2;
	int out_width0 = (int)floor(1.0+(log(N)/log(10)));

	std::cout << "\n[A." << std::setw(out_width0) <<0<< "]  length of a sequence                  = " << N; std::cout << "\n[A." << std::setw(out_width0) <<1<< "]  length of generating vector           = " << len_C_gen;
	std::cout << "\n[A." << std::setw(out_width0) <<2<< "]  left index linear subsequence         = " << left_index;
	std::cout << "\n[A." << std::setw(out_width0) <<3<< "]  right index linear subsequene         = " << right_index;
	std::cout << "\n[A." << std::setw(out_width0) <<4<< "]  number of entries of a table          = " << number_of_entries << " (analytic count)\n";


	std::cout << "\n[B." << std::setw(out_width0) <<0<< "]  Wall time trivial algo single thread  = " << time_triv_st << " ms";
	std::cout << "\n[B." << std::setw(out_width0) <<1<< "]  Wall time trivial algo multi thread   = " << time_triv_mt << " ms";
	std::cout << "\n[B." << std::setw(out_width0) <<2<< "]  Wall time faster algo single thread   = " << time_fast_st << " ms";
	std::cout << "\n[B." << std::setw(out_width0) <<3<< "]  Wall time faster algo multi thread    = " << time_fast_mt << " ms\n";

}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void print_initial_data(void)
{
	/*
	   To output the generating vector used to hide the linear subsequence,
	   and the input sequence itself.
	   */

	//std::cout << "\n\n";
	std::cout << "Original sequence (is identical to row 1 of triangular table), length = " << V0.size() << ", sequence:\n     ";
	for(unsigned int l1 = 0; l1<V0.size(); l1++)
		std::cout << V0[l1] << " ";

	std::cout << "\n\n";
	std::cout << "Generating vector, length = " << C_gen.size() << ", vector:\n     ";
	for(unsigned int l1 = 0; l1<C_gen.size(); l1++)
		std::cout << C_gen[l1] << " ";

	std::cout << "\n\n";
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void print_mat(std::vector<std::bitset<MAT_MAX_SIZE>> M,int max_dim,int n,int out_width1)
{
	/*
	   DESCRIPTION
	   print_ntl_gf2_mat prints on the standard output M with max_dim rows, and n columns.
	   out_width1 is to control the number of caracters for the index of the row.
	   */

	std::cout << "\n";
	std::cout << std::setw(out_width1) << 0 << " : ";
	for(int c0=0;c0 < n;c0++)
		std::cout << M[0][c0] << " ";

	std::cout << "\n";
	for(int r0=1;r0<max_dim+1;r0++)
	{
		std::cout << std::setw(out_width1) << r0 << " : ";
		for(int c0=0;c0 < r0-1;c0++)
			std::cout << "  ";

		for(int c0=r0-1;c0 <= n-r0;c0++)
			std::cout << M[r0][c0] << " ";
		std::cout << "\n";
	}
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

/* UTILS FUNCTIONS */
void generate_initial_data()
{
	/*
	   Create an input sequence and a generating vector for the LFSR used to create
	   the linear subsequence. The linear subsequence begins at position left_index and 
	   ends at position right_index.

	   Thus 0 < left_index < right_index <= n-1, where n is the length of the input sequence.
	   */

	std::random_device r_dev;
	std::mt19937_64 mt(r_dev());
	std::bernoulli_distribution dist(0.5);

	for(unsigned int i1=0;i1<C_gen.size();i1++)
	{
		C_gen[i1] = dist(mt);
	}
	C_gen[C_gen.size()-1] = true;//make sure the generating vector is nonzero

	for(int x1=0;x1<N;x1++)
	{
		V0[x1] =  dist(mt);
	}

	for(int x1=left_index ;  x1 < right_index; x1++)
	{
		bool tmpsum = false;
		for(size_t y1 = 0; y1 < C_gen.size()-1; y1++)
		{
			tmpsum = tmpsum || (C_gen[C_gen.size()-2-y1 ] && V0[x1 -1 - y1]);//scalar product computation --- linear subsequence
		}
		V0[x1] = tmpsum;//scalar product value
	}
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void assign_GF2(std::vector<std::bitset<MAT_MAX_SIZE>> & M, bool value, int j, int i)
{
	//A very nice idea found by Bastien Rigault! 

	if (MULTI_THREAD_ON) {
		determinant_mutex.lock();
		M[j][i] = value;
		determinant_mutex.unlock();
	}
	else
		M[j][i] = value;
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void init_experiment(struct experiment & experiment)
{
	//Initialize trivial or fast methods. For the trivial, some data is created even though not required. Those data are marked by *fast*

	experiment.M = std::vector<std::bitset<MAT_MAX_SIZE>>(max_dim+1, 0); //number of rows = max_dim + 1 (an extra one required for initializing the top part)
	experiment.flags_M = std::vector<std::bitset<MAT_MAX_SIZE>>(max_dim+1, 0); //*fast*
	experiment.AspTL[0] = std::vector<std::bitset<MAT_MAX_SIZE>>(1, 0); //*fast*

	experiment.AspTL[0][0][0]=false;//*fast*, initialized however never needed

	experiment.M[0].set();
	experiment.M[1] |= V0;
	for(int c1=0;c1<N;c1++)
	{
		experiment.flags_M[0][c1]=true; //*fast*
		experiment.flags_M[1][c1]=true; //*fast*
	}

	for(int t=1;t<max_dim+1;t++)
	{
		experiment.AspTL[t] = std::vector<std::bitset<MAT_MAX_SIZE>>(t); //*fast*
	}
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

/* ALGORITHM FUNCTIONS */

/*
   DESCRIPTION
   partition_thread_load for a given j-th row, determines each thread load (start index and range)
   Returned value is a list of pair, pair.first corresponds to start index and pair.second the range.
   The size of the list might be smaller than the number_of_threads value in cases when the number of
   determinant to compute is too small.
   Returns an empty list if the computation of the j-th row should be handle by one (current) thread
   */

std::vector<std::pair<int,int>> partition_thread_load(int j)
{
	/*
	   Distribute the number of tasks uniformly (or almost uniformly) across processes.
	   A task is computing a determinant.
	   */

	std::vector<std::pair<int, int>> result;

	if(N-2*j+1 > number_of_threads)//then uniformly distribute number of tasks on the maximal number of allowed threads
	{
		// Uniform share strategy
		int quo = (N-2*j+1) / number_of_threads;
		int rem = (N-2*j+1) % number_of_threads;

		for(long a = 0 ; a < number_of_threads; a++)
		{
			result.push_back(std::pair<int,int>(j-1+a*quo,quo+1));//for thread 0 to rem-1 inclusively, give them 1 more task then the quotient as the number of tasks
		}

		if (rem != 0)
			result.push_back(std::pair<int,int>(j-1+(number_of_threads*quo), rem+1));
	}
	else//then distribute one task per thread on n-2*j+2 threads < number_of_threads
	{
		for(long a = (j-1) ; a <= N-j ; a++)//safe no out of bound since n-2*j+2 <= number_of_threads
		{
			result.push_back(std::pair<int,int>(a,1));
		}
	}


	if(result.empty()){std::cerr << "List for the distribution of workloads is empty. EXIT\n";exit(-1);}//This should never happen. The only this could happen is if n-2j+2 = 0 but this impossible.

	return result;

}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void s_f(struct experiment & experiment, int j)
{
	/*
DESCRIPTION:
Fill the triangular table with zeros when the top part of a square is detected.
By top part, for position indices i_0, i_1 such that j-1 < i_0<i_1 <n-j, and j>=2,
If a run of non-zero elements is detected from i_0-1 to to i_1+1 on level j-2,
and a run of zero elements is detected from i_0 to i_1 on level j-1, then fill the table
below accordingly.
*/

	/*
	   Use knowledge of table with rows 0-th to (j-1)-th included to build square of zeros

	   StartAtLeft & EndInMiddle
	   StartAtLeft & EndAtRight
	   StartInMiddle & EndInMiddle
	   StartInMiddle & EndAtRight
	   */

	/*
	   (...)                                   ...                 ...                                         (...)

	   (j-2,j-3) |*| (j-2,j-2) (j-2,j-1) ....................... (j-2, n-j) (j-2,n-j+1) |*| (j,n-j+2) < n-(j-2)+1=n-j+3

	   (j-1,j-2) (j-1,j-1) ....................... (j-1, n-j) (j-1,n-j+1)               < n-(j-1)+1=n-j+2

	   --------------------------->  (j  ,j-1) ....................... (j  , n-j)                           < n-(j  )+1=n-j+1

*/

	int i = j-1;
	//struct experiment & experiment = const_cast<struct experiment &>(exp);

	while(i<N-j+1)//using knowledge of (j-2)th, (j-1)th rows, stop when right position = n-j+1 excluded
	{
		int lp = 0;//index of the left non-zero element
		int rp = 0;//index of the right non-zero element

		bool chkL1 = (i==j-1) && (experiment.M[j-2][i-1]) && (!experiment.M[j-1][i-1]) && (experiment.M[j-2][i]) && (!experiment.M[j-1][i]);

		bool chkL2 = (experiment.M[j-2][i-1]) && (experiment.M[j-2][i]) && (experiment.M[j-1][i-1]) && (!experiment.M[j-1][i]);

		if(chkL1||chkL2)
		{
			lp = i;
			int t1 = lp+1;

			while(t1 < N-j+1)
			{
				bool chkR1 = (experiment.M[j-2][t1-1]) && (experiment.M[j-2][t1]) && (!experiment.M[j-1][t1-1]) && (experiment.M[j-1][t1]);
				bool chkR2 = (t1==N-j) && (experiment.M[j-2][t1]) && (!experiment.M[j-1][t1]) && (experiment.M[j-2][t1+1]) && (!experiment.M[j-1][t1+1]);

				if(chkR1||chkR2)
				{
					rp = t1-1;
					break;
				}
				else
					t1 = t1 + 1;
			}
			if(rp>lp)
			{
				int stop_dim;
				if(max_dim<=(j+rp-lp-1))
					stop_dim = max_dim;
				else
					stop_dim = (j+rp-lp-1);

				for(int x_j=j;x_j<stop_dim;x_j++)//inside square excluding boundaries, square of zeros top row = j-1 which has been completed, thus begin at j
				{
					int x_i_l,x_i_r;//extreme possible values of i as the filling process goes from top part to bottom part

					if(lp>=x_j-1)
						x_i_l = lp;
					else
						x_i_l = x_j-1;

					if(rp<=N-x_j)
						x_i_r = rp;
					else
						x_i_r = N-x_j;

					for(int x_i=x_i_l;x_i<=x_i_r;x_i++)
					{
						experiment.M[x_j][x_i] = false;
						experiment.flags_M[x_j][x_i] = true;
					}
				}
				i = rp;
			}
			else
				i = i + 1;
		}
		else
			i = i+1;
	}
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

bool solve_eq_for_lower_corner(struct experiment & experiment, int j, int i, int q)
{

	/*
	   q is the minimum length of the grid required to solve for d_{i,j}.
	   The Main matrix has size therefore (q+1) X (q+1).
	   The Top/Bottom & Left/Right matrices used to solve for d_{i,j} has size q X q.
	   */


	bool ret = false;

	std::vector<std::bitset<MAT_MAX_SIZE>> T(q, 0); //temporary matrix of size q X q for a minor obtained from the (q+1) X (q+1) Main matrix
	for(int g=0;g<q;g++)//index mineure
	{
		for(int h=0;h<q+1;h++)//index rangee de la mineure
		{
			if(g<h)
			{
				for(int k=0;k<q;k++)//index colonne de la mineure
					T[h-1][k]=experiment.M[j-2*q+h+k][i-h+k];
			}
			else if (g>h)
			{
				for(int k=0;k<q;k++)
					T[h][k]=experiment.M[j-2*q+h+k][i-h+k];
			}
			//else it is not part of the expansion of the determinant
		}
		ret = ret ^ (experiment.M[j-q+g][i+q-g] & GF2_Utils::det_b(std::vector<std::bitset<MAT_MAX_SIZE>>(T), q));
	}
	return ret;
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

bool chk_conds_for_solvability(struct experiment & experiment, int j, int i, int & effective_length)
{
	/*
	   Checks the hypotheses of the conjecture in technical reference given above.
	   */

	if (MULTI_THREAD_ON)
		asptl_mutex.lock();
	bool ret=false;

	for(int w1=2;w1<=max_len_side_grid;w1++)//w1 is the length of the side which is growing so that the main matrix has size (w1+1)x(w1+1)
	{
		for(int rx=0;rx<w1;rx++)
		{
			for(int cx=0;cx<w1;cx++)
			{
				experiment.AspTL[w1][rx][cx] = experiment.M[j-2*w1+rx+cx][i-rx+cx];
			}
		}

		if ( (GF2_Utils::det_b(std::vector<std::bitset<MAT_MAX_SIZE>>(experiment.AspTL[w1]), w1)) && (!experiment.M[j-w1][i]) )/******/
		{
			ret=true;
			effective_length=w1;//length of a side of the grid so the main square matrix of size (w1+1) X (w1+1)	  
			break;
		}
	}

	if (MULTI_THREAD_ON)
		asptl_mutex.unlock();

	return ret;
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void d_c_s(struct experiment & experiment, int j, int start, int range)
{
	/*
	   This function computes determinants using results from the technical reference above.
	   This function must be called after the filling of zero-square is completed.
	   */

	int i = start;

	while(i < start+range)
	{
		int effective_len;

		if(!experiment.flags_M[j][i])
		{
			/************/

			if( (j>=2) && (experiment.M[j-2][i]) )//North-South-East-West
			{
				experiment.M[j][i] = (experiment.M[j-1][i-1] & experiment.M[j-1][i+1]) ^ experiment.M[j-1][i];
				experiment.flags_M[j][i]=true;
			}
			else if ( (j>=2*max_len_side_grid) && chk_conds_for_solvability(experiment, j, i, effective_len) )
			{
				experiment.M[j][i] = solve_eq_for_lower_corner(experiment, j, i, effective_len);
				experiment.flags_M[j][i] = true;
			}
			else
			{
				int t_x=1;

				while(t_x<=j)
				{
					if(!experiment.M[j-t_x][i])//start check from j-1 to 0 (experiment.M[0][i] == 1 by definition) ( level of sequence )
						t_x+=1;
					else
						break;//when here, experiment.M[j-t_x][i]!=0 AND experiment.M[j-1, j-2, ..., j-(t_x-1)][i]==0
				}
				//Note that t_x can never be 2 actually for it is handled by the North-South-East-West 
				if(t_x==1)//explicit computation
				{
					std::vector<std::bitset<MAT_MAX_SIZE>> tmp(j, 0);

					for(int r = 0;r<j ; r++)
					{
						for(int c = 0; c<j ; c++)
						{
							//tmp[c][r] = experiment.M[t_x][i-r+c];//==V0[i-j+r+c+1];
							tmp[c][r] = V0[i-j+r+c+1];
						}
					}

					experiment.M[j][i] = GF2_Utils::det_b(tmp, j);
					experiment.flags_M[j][i]=true;

				}
				else//computation by cross identities
				{
					std::vector<std::bitset<MAT_MAX_SIZE>> tmp(t_x, 0);

					for(int r=0;r<t_x;r++)
					{
						for(int c=0;c<t_x;c++)
							tmp[r][c]=experiment.M[j-(t_x-1)][i+c-r];
					}

					experiment.M[j][i] = GF2_Utils::det_b(tmp, t_x);
					experiment.flags_M[j][i]=true;

				}
			}
		}
		i=i+1;
	}
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void d_c_s_mt(struct experiment & experiment, int j, int start, int range)
{
	//Adaptation of d_c_s for multithread thread purpose by Bastien Rigault! Absolutely wonderful!

	int i = start;
	std::vector<struct determinant_result> results;

	while(i < start+range)
	{
		int effective_len;

		if(!experiment.flags_M[j][i])
		{
			/************/

			if( (j>=2) && (experiment.M[j-2][i]) )//North-South-East-West
			{
				bool determinant = (experiment.M[j-1][i-1] & experiment.M[j-1][i+1]) ^ experiment.M[j-1][i];
				assign_GF2(experiment.M, determinant, j, i);
				experiment.flags_M[j][i] = true;
			}
			else if ( (j>=2*max_len_side_grid) && chk_conds_for_solvability(experiment, j, i, effective_len) )
			{
				bool determinant = solve_eq_for_lower_corner(experiment, j, i, effective_len);
				assign_GF2(experiment.M, determinant, j, i);
				experiment.flags_M[j][i] = true;

			}
			else
			{
				int t_x=1;

				while(t_x<=j)
				{
					if(!experiment.M[j-t_x][i])//start check from j-1 to 0 (experiment.M[0][i] == 1 by definition) ( level of sequence )
						t_x+=1;
					else
						break;//when here, experiment.M[j-t_x][i]!=0 AND experiment.M[j-1, j-2, ..., j-(t_x-1)][i]==0
				}
				//Note that t_x can never be 2 actually for it is handled by the North-South-East-West 
				if(t_x==1)//explicit computation
				{
					std::vector<std::bitset<MAT_MAX_SIZE>> tmp(j, 0);

					for(int r = 0;r<j ; r++)
					{
						for(int c = 0; c<j ; c++)
						{
							//tmp[c][r] = experiment.M[t_x][i-r+c];//==V0[i-j+r+c+1];
							tmp[c][r] = V0[i-j+r+c+1];
						}
					}
					bool determinant = GF2_Utils::det_b(tmp, j);
					assign_GF2(experiment.M, determinant, j, i);
					experiment.flags_M[j][i] = true;

				}
				else//computation by cross identities
				{
					std::vector<std::bitset<MAT_MAX_SIZE>> tmp(t_x, 0);

					for(int r=0;r<t_x;r++)
					{
						for(int c=0;c<t_x;c++)
							tmp[r][c]=experiment.M[j-(t_x-1)][i+c-r];
					}

					bool determinant = GF2_Utils::det_b(tmp, t_x);
					assign_GF2(experiment.M, determinant, j, i);
					experiment.flags_M[j][i] = true;
				}
			}
		}
		i=i+1;
	}
}


/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void solve_j_trivial(struct experiment & experiment, int j, int start, int range)
{
	/*
	   Compute determinants for a given row of the triangular difference table.
	   The given row is the j-th row.
	   */

	for(int i = start; i<start+range; i++)//pour single thread est equivalent a i=j-1 ... i < j-1+n-2*j+2=n-j+1
	{
		std::vector<std::bitset<MAT_MAX_SIZE>> tmp(j, 0);

		for( int r = 0; r<j ; r++)
		{
			for( int c = 0; c<j ; c++)
			{
				tmp[r][c] = V0[i+1-j+r+c];
			}
		}
		bool determinant = GF2_Utils::det_b(tmp, j);
		assign_GF2(experiment.M, determinant, j, i);
	}
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void solve_trivial(struct experiment & experiment)
{

	if (!MULTI_THREAD_ON)
	{
		// Single thread solving
		for (int j=2; j<max_dim+1; j++)
			solve_j_trivial(experiment, j, j-1, N-2*j+2);
	}
	else
	{
		// Multi thread solving

		for (int j=2; j<max_dim+1; j++)
		{
			std::vector<std::pair<int, int>> thread_load = partition_thread_load(j);
			std::thread * thread_pool = new std::thread[thread_load.size()];

			//Uncomment lines followed by PRINT_WORKLOAD to see work loads and left inclusive boundaries for each sub-interval. Suggestion set nb_trials=1.
			//std::cout << "row # " << j << " : " << j-1 << " <= i < " << n-j+1 << " :: ";//PRINT_WORKLOAD

			for (unsigned int k=0; k<thread_load.size(); k++)
			{
				//std::cout << thread_load[k].first << ", " << thread_load[k].second << " | ";//PRINT_WORKLOAD
				thread_pool[k] = std::thread(&solve_j_trivial, std::ref(experiment), j, thread_load[k].first, thread_load[k].second);
			}
			//std::cout << "\n";//PRINT_WORKLOAD

			for (unsigned int k=0; k<thread_load.size(); k++)
				thread_pool[k].join();

			delete [] thread_pool;
		}
	}
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

void solve_fast(struct experiment & experiment)
{
	if (!MULTI_THREAD_ON)
	{
		// Single thread solving
		for (int j=2; j<max_dim+1; j++)
		{
			//Zero-square filling 
			s_f(experiment, j);
			//Use other tricks NSEW, cross shape (Bareiss) identities, conjecture, etc.
			d_c_s(experiment, j, j-1, N-2*j+2);
		}
	}
	else
	{
		// Multi thread solving


		for (int j=2; j<max_dim+1; j++)
		{
			std::vector<std::pair<int, int>> thread_load = partition_thread_load(j);//We can't have optimized strategy for the workloads that uses the flags before the square filling is over. For the uniform stategy it does not matter.
			std::thread* thread_pool = new std::thread[thread_load.size()];

			//Zero-square filling (Note that s_f cannot be done in parallel.) 
			s_f(experiment, j);

			//Use other tricks NSEW, cross shape (Bareiss) identities, conjecture, etc in parallel.
			for (size_t k=0; k<thread_load.size(); k++)
			{
				//Make sure that inside solve_j_fast, the square filling is done mono thread only as it cannot be multi threaded.
				//This explains why there were lots of discrepancies for the optimized strategy which I removed here.
				thread_pool[k] = std::thread(&d_c_s_mt, std::ref(experiment), j, thread_load[k].first, thread_load[k].second);
			}

			for (size_t k=0; k<thread_load.size(); k++)
				thread_pool[k].join();

			delete [] thread_pool;
		}
	}
}

/**** ***** ***** ***** ***** ***** ***** ***** ***** *****/

int main(void)
{
	std::cout.precision(6);
	std::cout.setf( std::ios::fixed, std::ios::floatfield );

	if(right_index>N){std::cerr << "\n\nOut of bound\n\n"; exit(-1);}

	// Experiments declarations
	struct experiment 
		trivial_experiment, 
		trivial_experiment_mt, 
		fast_experiment,
		fast_experiment_mt;


	double time_triv_st,time_triv_mt,time_fast_st,time_fast_mt;// for timing
	double sum_of_time_differences[4]={0.0, 0.0, 0.0, 0.0};//[0] trivial mono, [1] fast mono, [2] trivial multi, [3], fast multi 

	for(long ct_trial = 0; ct_trial<nb_trials;ct_trial++)
	{ 
		std::cout << "Test # " << ct_trial+1 << " / " << nb_trials << " ";
		std::cout.flush();

		// Initialize data
		generate_initial_data();
		//print_initial_data();

		init_experiment(trivial_experiment);
		init_experiment(trivial_experiment_mt);
		init_experiment(fast_experiment);
		init_experiment(fast_experiment_mt);

		// Run experiments

		MULTI_THREAD_ON = false;

		auto t_start = std::chrono::high_resolution_clock::now();
		solve_trivial(trivial_experiment);
		auto t_end = std::chrono::high_resolution_clock::now();
		time_triv_st = std::chrono::duration<double, std::milli>(t_end-t_start).count();
		sum_of_time_differences[0] += time_triv_st;

		t_start = std::chrono::high_resolution_clock::now();
		solve_fast(fast_experiment);
		t_end = std::chrono::high_resolution_clock::now();
		time_fast_st = std::chrono::duration<double, std::milli>(t_end-t_start).count();
		sum_of_time_differences[1] += time_fast_st;

		MULTI_THREAD_ON = true;

		t_start = std::chrono::high_resolution_clock::now();
		solve_trivial(trivial_experiment_mt);
		t_end = std::chrono::high_resolution_clock::now();
		time_triv_mt = std::chrono::duration<double, std::milli>(t_end-t_start).count();
		sum_of_time_differences[2] += time_triv_mt;

		t_start = std::chrono::high_resolution_clock::now();
		solve_fast(fast_experiment_mt);
		t_end = std::chrono::high_resolution_clock::now();
		time_fast_mt = std::chrono::duration<double, std::milli>(t_end-t_start).count();
		sum_of_time_differences[3] += time_fast_mt;

		// Prints
		//print_stats(time_triv_st,time_triv_mt,time_fast_st,time_fast_mt);

		int nb_differences;
		// Check results against trivial mono thread
		if (GF2_Utils::chk_triangular_tables_not_equal(trivial_experiment.M, fast_experiment.M, max_dim, N))
		{
			nb_differences = GF2_Utils::how_many_differences(trivial_experiment.M, fast_experiment.M, max_dim, N);
			std::cerr << "\n\n*****Mismatches between trivial single threaded and fast single threaded.*****\n";
			std::cerr << "#####Number of mismatches = " << nb_differences << "\n";
			std::cerr << "EXIT - error - mismatch between mono thread trivial and mono thread fast\n";

			std::cerr << "\n";
			print_mat(trivial_experiment.M, max_dim, N, 1);
			std::cerr << "\n";
			print_mat(fast_experiment.M, max_dim, N, 1);
			exit(-1);
		}

		if (GF2_Utils::chk_triangular_tables_not_equal(trivial_experiment.M, trivial_experiment_mt.M, max_dim, N))
		{
			nb_differences = GF2_Utils::how_many_differences(trivial_experiment.M, trivial_experiment_mt.M, max_dim, N);
			std::cerr << "\n\n*****Mismatches trivial single threaded and trivial multi threaded.*****\n";
			std::cerr << "#####Number of mismatches = " << nb_differences << "\n";
			std::cerr << "EXIT - error - mismatch between mono thread trivial and multi thread trivial\n";

			std::cerr << "\n";
			print_mat(trivial_experiment.M, max_dim, N, 1);
			std::cerr << "\n";
			print_mat(trivial_experiment_mt.M, max_dim, N, 1);
			exit(-1);
		}

		if (GF2_Utils::chk_triangular_tables_not_equal(trivial_experiment.M, fast_experiment_mt.M, max_dim, N))
		{
			nb_differences = GF2_Utils::how_many_differences(trivial_experiment.M, fast_experiment_mt.M, max_dim, N);
			std::cerr << "\n\n*****Mismatches trivial single threaded and fast multi threaded.*****\n";
			std::cerr << "#####Number of mismatches = " << nb_differences << "\n";
			std::cerr << "EXIT - error - mismatch between mono thread trivial and multi thread fast\n";

			std::cerr << "\n";
			print_mat(trivial_experiment.M, max_dim, N, 1);
			std::cerr << "\n";
			print_mat(fast_experiment_mt.M, max_dim, N, 1);
			exit(-1);
		}


		std::cout << "  ***** " << sum_of_time_differences[0]/(double)(ct_trial+1) << " " << sum_of_time_differences[2]/(double)(ct_trial+1) << " " << sum_of_time_differences[1]/(double)(ct_trial+1) << " " << sum_of_time_differences[3]/(double)(ct_trial+1) << "\n"; 
	}

	std::cout << "\n\n***** ***** ***** ***** ***** ***** ***** ***** ***** *****\n"
		<< "\nNumber of threads = " << number_of_threads
		<< "\nLength of sequences = " << N
		<< "\nLength of generating vector = " <<  len_C_gen
		<< "\nNumber of entries = " << ((N%2)==0)*((N/2)-1)*(N/2)+((N%2)==1)*(N/2)*(N-1)/2
		<< "\nHidden linear subsequence starting index = " << left_index
		<< "\nHidden linear subsequence ending index = " << right_index << "\n";

	std::cout << "\nThe following averages of running times are computed from a sample of size " << nb_trials << ".\n"; 
	std::cout << "\n   Average time trivial method single thread = " << sum_of_time_differences[0]/(double)nb_trials;
	std::cout << "\n   Average time trivial method multi thread  = " << sum_of_time_differences[2]/(double)nb_trials;
	std::cout << "\n   Average time fast method single thread    = " << sum_of_time_differences[1]/(double)nb_trials;
	std::cout << "\n   Average time fast method multi thread     = " << sum_of_time_differences[3]/(double)nb_trials;

	std::cout << "\n\n";
	return 0;
}
