// MIMO Detector Simulation
// Date: 2022-04-12
// Author: CYJ

#include <ctime>
#include <bitset>
#include <random>
#include "MIMO_4x4.h"

using std::norm;

const double NOISE_SIGMA[12] = {0.5006, 0.3976, 0.3159, 0.2236, 0.1583, 0.1257,
							    0.0999, 0.0707, 0.0562, 0.0398, 0.0282, 0.0224};
const int SNR[12] = {3, 5, 7, 10, 13, 15, 17, 20, 22, 25, 28, 30};
const double H_X_SIGMA = 0.35355;		// 1/sqrt(8)

complex<int> bit2int_dec(unsigned char);
unsigned char int2bit_dec(complex<int>);
void channel(complex<double>** , complex<double>* , complex<double>*, std::normal_distribution<double>);
int check(unsigned char, unsigned char);

// generate uniform distribution random variable served as random seed
std::random_device rd;
// give to random number generator
std::mt19937 rand_gen(time(NULL));
// normal distribution with 0 mean and standard deviation = NOISE_SIGMA
std::normal_distribution<double> gaussCH(0.0, H_X_SIGMA);
std::uniform_int_distribution<int> xrand(0, 15);

int main(void)
{
	MIMO_4x4 module1;
	complex<double> **H, *y_normal, *x_normal;
	complex<int> *x_symbol, *x_dtc;
	double y_total = 0, x_total = 0, H_total = 0;
	int Bit_err = 0, total_Bit_err = 0, total_Sym_err = 0;
	int patherr = 0;
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
	
	int64_t TEST_NUM = 10000;
	
	auto tmn_cout_buf = cout.rdbuf(xout_verilog.rdbuf());
	cout.rdbuf(tmn_cout_buf);

	int SNR_idx = 11;

	// for (SNR_idx = 0; SNR_idx < 12; SNR_idx++) {
		// if (SNR_idx > 6) TEST_NUM = 200000;
		// else TEST_NUM = 5000;
		std::normal_distribution<double> gaussN(0.0, NOISE_SIGMA[SNR_idx]);
		printf("\nSNR: %2d \n", SNR[SNR_idx]);

		total_Bit_err = total_Sym_err = 0;

		ofstream FILE_xdtc(".\\data\\xdtc_cpp.txt");

		for (int64_t i = 0; i < TEST_NUM; i++) {
			for (int j = 0; j < N; j++) {
				x_bit[j] = xrand(rand_gen);
				x_symbol[j] = bit2int_dec(x_bit[j]);
				// std::bitset<4> bs(x_bit[j]);
				// cout << "bit: " << bs << "  --->   symbol: " << x_symbol[j] << endl;
				x_normal[j].real(x_symbol[j].real() / KMOD);
				x_normal[j].imag(x_symbol[j].imag() / KMOD);
			}

			// cout << endl;

			// -------- floating point ---------
			// channel(H, y_normal, x_normal, gaussN);			// pass x through the simulated channel
			// module1.setH(H);									// Feed the H and Y to detector
			// module1.setY(y_normal);
			// module1.detect();									// detect	
			// module1.getX(x_dtc);								// get the solved x symbol
			// module1.output_X_CSV(x_symbol);
			// module1.output_R_tmn();

			// -------- fixed point ---------
			channel(H, y_normal, x_normal, gaussN);				// pass x through the simulated channel
			module1.q_input_complex2int(H, y_normal);
			module1.q_detect();
			module1.getX(x_dtc);								// get the solved x symbol
			
			// check result
			Bit_err = 0;
			for (int j = N - 1; j >= 0; j--) {
				x_bit_dtc[j] = int2bit_dec(x_dtc[j]);			// convert to binary representation
				Bit_err += check(x_bit[j], x_bit_dtc[j]);		// compare the result with original input x
				// std::bitset<4> bs(x_bit[j]);
				// cout  << bs << ' ';
				xout_verilog << int(x_bit[j]) << ' ';		// export the x
				FILE_xdtc << int(x_bit_dtc[j]) << ' ';
			}
			FILE_xdtc << endl;
			if (Bit_err > 0) patherr++;
			total_Bit_err += Bit_err;						// sum up the number of error bits
			xout_verilog << endl;

			if (i % (TEST_NUM / 10) == 0) printf("-");
		}
		// cout << "E[|x|] = " << x_total / TEST_NUM << endl;
		// cout << "E[|y|] = " << y_total / TEST_NUM << endl;
		// cout << "E[|H|] = " << H_total / TEST_NUM << endl;

		cout << endl;
		printf("BER: %f (%d/%d)\n", 
				total_Bit_err * 1.0 / (4 * M * TEST_NUM), total_Bit_err, 4 * M * TEST_NUM);
		// printf("SER: %9f\n", total_Sym_err * 1.0 / (M * TEST_NUM));
		printf("PER: %f (%d/%d)\n", patherr * 1.0 / TEST_NUM, patherr, TEST_NUM);

	// }


	return 0;
}

void read_testHY(complex<double>** Hcpx, complex<double>* Ycpx)
{	
	int H[8][8], Y[8];
	ifstream file_HY("testHY.txt");

	for (int i = 0; i < 2 * N; i++)
		for (int j = 0; j < 2 * M; j++)
			file_HY >> H[i][j];
	
	for (int i = 0; i < 2 * M; i++) 
		file_HY >> Y[i];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			Hcpx[i][j].real(H[i * 2][j * 2]);
			Hcpx[i][j].imag(H[i * 2 + 1][j * 2]);
		}
		Ycpx[i].real(Y[i * 2]);
		Ycpx[i].imag(Y[i * 2 + 1]);
	}
}

void channel(complex<double>** H, complex<double>* y, complex<double>* x, std::normal_distribution<double> gaussN)
{
	complex<double> n[4];
	// std::normal_distribution<double> gaussN(0.0, noise_sigma);

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