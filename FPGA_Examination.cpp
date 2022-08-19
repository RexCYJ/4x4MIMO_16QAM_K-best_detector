#include <fstream>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

// #define TESTNUM 200000
// #define SNR 25

int check(int, int);

int main (void)
{
	int in;
	int x_orig, x_fpga, x_fxp;
	int errbit, errbit_total = 0;
	int errpath_total = 0;
	int SNR, TESTNUM;
	double BER, PER;

	stringstream FPGA_datapath;
	stringstream ORIG_datapath;
	stringstream BER_datapath;

	cout << "Input SNR: ";
	cin >> SNR;
	cout << "Input number of data: ";
	cin >> TESTNUM;
	cout << "SNR: " << SNR << endl;
	FPGA_datapath << ".\\data\\FPGA_result\\" << SNR << "dB\\X.out";
	ORIG_datapath << ".\\data\\FPGA_result\\" << SNR << "dB\\xdat.txt";
	BER_datapath << ".\\data\\FPGA_result\\" << SNR << "dB\\BER.csv";

	ifstream FILE_fpga(FPGA_datapath.str());
	ifstream FILE_orig(ORIG_datapath.str());
	ofstream FILE_BER(BER_datapath.str());
	ifstream FILE_fxp(".\\data\\xdtc_cpp.txt");
	
	cout << "Source path: " << FPGA_datapath.str() << endl;

	for (int i = 0; i < 8; i++)
		FILE_fpga >> x_fpga;

	for (int i = 0; i < TESTNUM; i++) {
		errbit = 0;
		x_orig = 0;
		// x_fxp = 0;
		FILE_fpga >> x_fpga;
		for (int j = 0; j < 4; j++) {
			FILE_orig >> in;
			x_orig = (x_orig << 4) | in;
		}
		// for (int j = 0; j < 4; j++) {
		// 	FILE_fxp >> in;
		// 	x_fxp = (x_fxp << 4) | in;
		// }
		errbit += check(x_orig, x_fpga);
		errbit_total += errbit;
		// if (check(x_fxp, x_fpga)) {
		// 	printf("mismatch: %6d FPGA error bits: %2d FXP error bit: %2d\n",
		// 			i + 1, errbit, check(x_orig, x_fxp));
		// }
		if (errbit > 0)
			errpath_total++;
		
		FILE_BER << (errbit_total / ((i + 1) * 16.0)) << ',' << ' ';
	}

	BER = errbit_total / (TESTNUM * 16.0);
	PER = errpath_total * 1.0 / TESTNUM;
	printf("BER: %f (%d/%d)\n", BER, errbit_total, TESTNUM * 16);
	printf("PER: %f (%d/%d)\n", PER, errpath_total, TESTNUM);
	
	return 0;
}

int check(int orig, int solve)
{
	int err = 0;
	int errbit;

	errbit = orig ^ solve;
	for (int i = 0; i < 16; i++) {
		err += errbit & 1;
		errbit = errbit >> 1;
	}

	return err;		// return error bit
}