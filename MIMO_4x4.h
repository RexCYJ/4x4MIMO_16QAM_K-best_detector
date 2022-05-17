#ifndef MIMO_H
#define MIMO_H

#include <complex>
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::complex;
using std::norm;
using std::ofstream;
using std::vector;
using std::pair;
using std::make_pair;

#define N 4
#define M 4

ofstream out("matrix.csv");
const int CORDIC_ITER = 15;					// Number of iteration
const double arctan2[CORDIC_ITER] 
		= {45.0, 26.565051, 14.036243, 7.125016, 3.576334, 
			1.789910, 0.895174, 0.447614, 0.223811, 0.111906, 
			0.055953, 0.027976, 0.013988, 0.006994, 0.003497};
const double An = 1.6467602;				// Rotation gain
const int QAM16_val[4] = {-3, -1, 1, 3};	// 16-QAM value
const double QAM16_normval[4] = {-0.94868, -0.31623, 0.31623, 0.94868};
const double KMOD = 3.16228;				// sqrt(10)

// An 4x4 (NxM) MIMO system
// Modulation: 16-QAM
class MIMO_4x4
{
	public:
		MIMO_4x4();
		void detect();
		void setH(complex<double>** );
		void setY(complex<double>* );
		void getX(complex<int>*);
		void output_X_CSV(complex<int>*);
		void output_R_tmn();
	private:
		complex<double> **H;
		double Hr[2 * N][2 * M];
		double R[2 * N][2 * M];
		double Yr[2 * N];
		double Yexp[2 * N];
		int xIndex[2 * M];
		int x[2 * M];
		complex<int> *xComplex;
		complex<double> *y;
		int decomposeType;
		
		void ML();
		void Kbest();
		void Kbest3();
		void Kbest_optiK();
		void decompositionHalf();
		void decompositionFull();
		double CORDIC(double &, double &, double &, int);
		void outputMatrix_CSV();
		void CORDIC_ROW(double* , double* , int, double);
		void xRearrange();
};

#include "MIMO_4x4.cpp"

#endif