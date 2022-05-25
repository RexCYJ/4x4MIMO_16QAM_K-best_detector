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

#define FRAC 12
#define INTE 2
#define ANGLE_FRAC 15

ofstream out("matrix.csv");

const int CORDIC_ITER = 9;					// Number of iteration
const double arctan2[15] 
		= {45.0, 26.565051, 14.036243, 7.125016, 3.576334, 
			1.789910, 0.895174, 0.447614, 0.223811, 0.111906, 
			0.055953, 0.027976, 0.013988, 0.006994, 0.003497};

const double An = 1.6467602;				// Rotation gain
const int q_Const = (1 / An) * (1 << FRAC);

const int QAM16_val[4] = {-3, -1, 1, 3};	// 16-QAM value
const double QAM16_normval[4] = {-0.94868, -0.31623, 0.31623, 0.94868};
const int32_t QAM16_q_normval[4] = {-3886, -1295, 1295, 3886};
const double KMOD = 3.16228;				// sqrt(10)

// An 4x4 (NxM) MIMO system
// Modulation: 16-QAM
class MIMO_4x4
{
	public:
		MIMO_4x4();
		void detect();
		void q_detect(complex<double> **, complex<double> *);	// Input H and Y
		void setH(complex<double>** );
		void setY(complex<double>* );
		void getX(complex<int>*);
		void output_X_CSV(complex<int>*);
		void output_RY_tmn();

	private:
		complex<double> **H;		// input H
		complex<double> *y;			// input y
		complex<int> *xComplex;		// output x
		double Hr[2 * N][2 * M];	// for ML
		double Yr[2 * N];			// 
		double R[2 * N][2 * M];		// for floating-point-decomposed result
		double Yexp[2 * N];			// 
		int xIndex[2 * M];			// indicate the rearranged column index
		int x[2 * M];				
		int16_t qH[2 * N][2 * N], qY[2 * N];
		int16_t qR[2 * N][2 * N], qYtrans[2 * N];		// for fixed-point-decomposed result

		double qR_flp[2 * N][2 * N], qYtrans_flp[2 * N];

		int decomposeType;

		void q_decompositionFull();
		void q_Kbest_optiK();
		void DecomposeError();
		
		void ML();
		void Kbest();
		void Kbest3();
		void Kbest_optiK();

		void decompositionHalf();
		void decompositionFull();
		double CORDIC(double &, double &, double &, int);
		void CORDIC_ROW(double* , double* , int, double);
		
		void q_Input(complex<double> **, complex<double> *);
		
		void xRearrange();
		void output_Matrix_CSV();
		void output_HY_tmn();
		void output_qRY_tmn();
		void output_qRY_flp_tmn();
};

MIMO_4x4::MIMO_4x4()
{
	H = new complex<double>*[N];
	H[0] = new complex<double>[M];
	H[1] = new complex<double>[M];
	H[2] = new complex<double>[M];
	H[3] = new complex<double>[M];
	y = new complex<double>[M];
	xComplex = new complex<int>[M];
}

#include "Detector_FLP.cpp"
#include "Detector_FXP.cpp"
#include "Decomposition_FLP.cpp"
#include "Decomposition_FXP.cpp"
#include "ports.cpp"

#endif