#include <fstream>
#include <iostream>

using std::ifstream;
using std::cout;
using std::endl;

#define TESTNUM 20

int check(unsigned char, unsigned char);

int main (void)
{
	ifstream FILE_veri(".\\data\\xdtc_vivado_dat.txt");
	ifstream FILE_orig(".\\data\\xdat.txt");
	// ifstream FILE_orig(".\\data\\xdtc_cpp.txt");
	int in;
	unsigned char x_veri, x_orig;
	int errbit, errbit_total = 0;
	int errpath_total = 0;
	double BER, PER;

	for (int i = 0; i < TESTNUM; i++) {
		errbit = 0;
		for (int j = 0; j < 4; j++) {
			FILE_veri >> in;
			x_veri = in;
			FILE_orig >> in;
			x_orig = in;
			errbit += check(x_orig, x_veri);
		}
		errbit_total += errbit;
		if (errbit > 0) {
			errpath_total++;
			cout << "err: " << i + 1 << endl;
		}
	}

	BER = errbit_total / (TESTNUM * 16.0);
	PER = errpath_total * 1.0 / TESTNUM;
	printf("BER: %f (%d/%d)\n", BER, errbit_total, TESTNUM * 16);
	printf("PER: %f (%d/%d)\n", PER, errpath_total, TESTNUM);
	
	return 0;
}

int check(unsigned char orig, unsigned char solve)
{
	int err = 0;
	unsigned char errbit;

	errbit = orig ^ solve;
	err = (errbit & 0x01) + ((errbit & 0x02) >> 1) + ((errbit & 0x04) >> 2) + ((errbit & 0x08) >> 3);

	return err;		// return error bit
}