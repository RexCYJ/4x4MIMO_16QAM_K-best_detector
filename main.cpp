// MIMO Detector Simulation
// Date: 2022-04-12
// Author: CYJ

#include <ctime>
#include <bitset>
#include <random>
#include "MIMO_4x4.h"

using std::norm;

const int64_t TEST_NUM = 50000;

const double NOISE_SIGMA = 0.0224;

const double NOISE_R_SIGMA = sqrt(NOISE_SIGMA * NOISE_SIGMA * 2);
const double H_X_SIGMA = 0.35355;		// 1/sqrt(8)

complex<int> bit2int_dec(unsigned char);
unsigned char int2bit_dec(complex<int>);
void channel(complex<double>** , complex<double>* , complex<double>*);
int check(unsigned char, unsigned char);

// generate uniform distribution random variable served as random seed
std::random_device rd;
// give to random number generator
std::mt19937 rand_gen(time(NULL));
// normal distribution with 0 mean and standard deviation = NOISE_SIGMA
std::normal_distribution<double> gaussN(0.0, NOISE_SIGMA);
std::normal_distribution<double> gaussCH(0.0, H_X_SIGMA);
std::uniform_int_distribution<int> xrand(0, 15);

int main(void)
{
	MIMO_4x4 module1;
	complex<double> **H, *y_normal, *x_normal;
	complex<int> *x_symbol, *x_dtc;
	double y_total = 0, x_total = 0, H_total = 0;
	int Bit_err[N] = {0}, total_Bit_err = 0, total_Sym_err = 0;
	unsigned char *x_bit, *x_bit_dtc;

	H = new complex<double>*[N];
	H[0] = new complex<double>[M];
	H[1] = new complex<double>[M];
	H[2] = new complex<double>[M];
	H[3] = new complex<double>[M];
	y_normal = new complex<double>[N];
	x_normal = new complex<double>[M];
	x_symbol = new complex<int>[M];
	x_dtc = new complex<int>[M];
	x_bit = new unsigned char[M];
	x_bit_dtc = new unsigned char[M];

	for (int64_t i = 0; i < TEST_NUM; i++) {
		// cout << endl << "ITER: " << i << " ___________________________\n" << endl;

		for (int j = 0; j < N; j++) {
			x_bit[j] = xrand(rand_gen);
			x_symbol[j] = bit2int_dec(x_bit[j]);
			x_normal[j].real(x_symbol[j].real() / KMOD);
			x_normal[j].imag(x_symbol[j].imag() / KMOD);
		}
		
		channel(H, y_normal, x_normal);						// pass x through the simulated channel
		
		// for (int l = 0; l < N; l++)
		// 	printf("%3d %3d ", x_symbol[l].real(), x_symbol[l].imag());
		// cout << endl;

		// -------- floating point ---------
		module1.setH(H);									// Feed the H and Y to detector
		module1.setY(y_normal);
		module1.detect();									// detect	
		module1.getX(x_dtc);								// get the solved x symbol
		// module1.output_X_CSV(x_symbol);
		// module1.output_R_tmn();

		// -------- fixed point ---------
		// module1.q_detect(H, y_normal);
		// module1.getX(x_dtc);								// get the solved x symbol

		// check result
		for (int j = 0; j < N; j++) {
			x_bit_dtc[j] = int2bit_dec(x_dtc[j]);
			Bit_err[j] = check(x_bit[j], x_bit_dtc[j]);
			total_Bit_err += Bit_err[j];
			total_Sym_err += (Bit_err[j] > 0)? 1:0;
		}

		// ERROR checking
		// if (Bit_err[0] || Bit_err[1]) {
		// 	cout << endl << "ERROR: " << endl;
			// std::bitset<4> bs1(x_bit[0]);
			// std::bitset<4> bs2(x_bit[1]);
			// cout << "bit1: " << bs1 << "  --->   symbol1: " << x_symbol[0] << endl;
			// cout << "bit2: " << bs2 << "  --->   symbol2: " << x_symbol[1] << endl;
		// 	cout << "X symbol detect: " << x_dtc[0] << endl;
		// 	cout << "X symbol detect: " << x_dtc[1] << endl;
		// 	std::bitset<4> bs3(x_bit_dtc[0]);
		// 	std::bitset<4> bs4(x_bit_dtc[1]);
		// 	cout << "X bit detect: " << bs3 << endl;
		// 	cout << "X bit detect: " << bs4 << endl;
		// }

		// y_total += sqrt(norm(y_normal[0]));
		// x_total += sqrt(norm(x_normal[0]));
		// H_total += sqrt(norm(H[0][0]) + norm(H[0][1]) + norm(H[0][2]) + norm(H[0][3]));

		if (i % (TEST_NUM / 10) == 0) printf("%2d\%\n", i / (TEST_NUM / 100));
	}

	// cout << "E[|x|] = " << x_total / TEST_NUM << endl;
	// cout << "E[|y|] = " << y_total / TEST_NUM << endl;
	// cout << "E[|H|] = " << H_total / TEST_NUM << endl;
	cout << endl;
	cout << "BER: " << total_Bit_err * 1.0 / (4 * M * TEST_NUM) << endl;
	cout << "Number of error bits " << total_Bit_err << endl;
	cout << "SER: " << total_Sym_err * 1.0 / (M * TEST_NUM) << endl;

	return 0;
}

void channel(complex<double>** H, complex<double>* y, complex<double>* x)
{
	complex<double> n[4];

	for (int i = 0; i < N; i++) {
		n[i] = complex<double>(gaussN(rand_gen), gaussN(rand_gen));
		y[i] = n[i];
		for (int j = 0; j < M; j++) {
			H[i][j] = complex<double>(gaussCH(rand_gen), gaussCH(rand_gen));
			y[i] += H[i][j] * x[j];
		}
	}
}

complex<int> bit2int_dec(unsigned char n)
{
	int x, y;

	x = n & 0x03;
	y = n >> 2;

	switch (x) {
		case 0: x = 3; break;
		case 1: x = 1; break;
		case 3: x = -1; break;
		case 2: x = -3; break;
	}
	
	switch (y) {
		case 0: y = 3; break;
		case 1: y = 1; break;
		case 3: y = -1; break;
		case 2: y = -3; break;
	}

	return complex<int>(x, y);
}

unsigned char int2bit_dec(complex<int> S)
{
	unsigned code;

	switch (S.imag()) {
		case 3: code = 0; break;
		case 1: code = 1; break;
		case -1: code = 3; break;
		case -3: code = 2; break;
	}
	
	code <<= 2;
	
	switch (S.real()) {
		case 3: code += 0; break;
		case 1: code += 1; break;
		case -1: code += 3; break;
		case -3: code += 2; break;
	}

	return code;	
}

int check(unsigned char orig, unsigned char solve)
{
	int err = 0;
	unsigned char errbit;

	errbit = orig ^ solve;
	err = (errbit & 0x01) + ((errbit & 0x02) >> 1) + ((errbit & 0x04) >> 2) + ((errbit & 0x08) >> 3);

	return err;		// return error bit
}