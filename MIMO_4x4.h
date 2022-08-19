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
using std::ifstream;
using std::ofstream;
using std::vector;
using std::pair;
using std::make_pair;

#define N 4
#define M 4

#define FRAC 12
#define COLLEN_FRAC 7
#define ERR_FRAC 12
#define INTE 2

ofstream out(".\\data\\matrix.csv");
ofstream out_verilog(".\\data\\dat.txt");
ofstream Rout_verilog(".\\data\\Rdat.txt");
ofstream out0_verilog(".\\data\\data0.txt");
ofstream out1_verilog(".\\data\\data1.txt");
ofstream out2_verilog(".\\data\\data2.txt");
ofstream out3_verilog(".\\data\\data3.txt");
ofstream xout_verilog(".\\data\\xdat.txt");
ifstream file_HYX(".\\data\\testHY.txt");

const int CORDIC_ITER = 9;					// Number of iteration
const double arctan2[15] 
		= {45.0, 26.565051, 14.036243, 7.125016, 3.576334, 
			1.789910, 0.895174, 0.447614, 0.223811, 0.111906, 
			0.055953, 0.027976, 0.013988, 0.006994, 0.003497};

const double An = 1.6467602;				// Rotation gain
const int q_Const = (1 / An) * (1 << FRAC);

const int QAM16_val[4] = {-3, -1, 1, 3};	// 16-QAM value
const double QAM16_normval[4] = {-0.94868, -0.31623, 0.31623, 0.94868};
int32_t QAM16_q_normval[4]; 	//  = {-3886, -1295, 1295, 3886};
const double KMOD = 3.16228;				// sqrt(10)

const int zigzag_path[8][3] = {
	{1, 2, 3}, {1, 2, 3}, {0, 2, 3}, {2, 0, 3},
	{1, 3, 0}, {3, 1, 0}, {2, 1, 0}, {2, 1, 0}
};

// An 4x4 (NxM) MIMO system
// Modulation: 16-QAM
class MIMO_4x4
{
	public:
		MIMO_4x4();
		void detect();
		void q_detect();	// Input H and Y
		void setH(complex<double>** );
		void setY(complex<double>* );
		void getX(complex<int>*);
		void output_X_CSV(complex<int>*);
		void output_RY_tmn();
		void q_input_complex2int(complex<double> **, complex<double> *);
		void q_input_verifyHYX(unsigned char *);

		int *columnOrder;

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
		int32_t qH[2 * N][2 * N], qY[2 * N];
		int32_t qR[2 * N][2 * N], qYtrans[2 * N];		// for fixed-point-decomposed result
		uint32_t colpower[2 * M];
		double qR_flp[2 * N][2 * N], qYtrans_flp[2 * N];

		int decomposeType;

		void q_decompositionFull();
		void q_decompositionFull_test();
		void q_Kbest_optiK();
		void q_Kbest4();
		void DecomposeError();
		
		void ML();
		void Kbest();
		void Kbest3();
		void Kbest4();
		void Kbest_optiK();

		void decompositionHalf();
		void decompositionFull();
		double CORDIC(double &, double &, double &, int);
		void CORDIC_ROW(double* , double* , int, double);
		
		void xRearrange();
		void output_Matrix_CSV();
		void output_qHY_tmn();
		void output_qRY_tmn();
		void output_qRY_flp_tmn();
		void output_qRY_verilog();
		void output_qHY_verilog();
		void output_qHY_verilog_batch();
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
	for (int i = 0; i < 4; i++)
		QAM16_q_normval[i] = (QAM16_val[i] / sqrt(10)) * (1 << FRAC);
}

#include "Detector_FLP.cpp"
#include "Detector_FXP.cpp"
#include "Decomposition_FLP.cpp"
#include "Decomposition_FXP.cpp"
#include "ports.cpp"

#endif