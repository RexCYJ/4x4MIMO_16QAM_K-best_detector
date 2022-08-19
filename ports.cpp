// I/O functions for MIMO_4x4.h

#include "MIMO_4x4.h"

void MIMO_4x4::output_qHY_tmn()
{
	cout << "  Y[] |   H[][]\n";
	for (int i = 0; i < 2 * N; i++) {
		printf("%5d | ", qY[i]);
		for (int j = 0; j < 2 * N; j++)
			printf("%5d ", qH[i][j]);
		cout << endl;
	}
	cout << endl;
}

void MIMO_4x4::output_qHY_verilog()
{
	// cout << "\n>> Output verilog format\n\n";
	auto cout_buf = cout.rdbuf(out_verilog.rdbuf());

	for (int i = 0; i < 2 * N; i++) {
		for (int j = 0; j < 2 * N; j++)
			cout << qH[i][j] << ' ';
		cout << endl;
	}

	for (int i = 0; i < 2 * N; i++) 
		cout << qY[i] << ' ';
	cout << endl << endl;

	cout.rdbuf(cout_buf);
}

void MIMO_4x4::output_qHY_verilog_batch()
{
	static bool is_first = true;

	// if (is_first) {
	// 	out3_verilog << '0' << endl;
	// 	out2_verilog << '0' << endl;
	// 	out1_verilog << '0' << endl;
	// 	out0_verilog << '0' << endl;
	// 	is_first = false;
	// } 
	// for (int i = 7; i >= 0; i -= 2) {
	// 	for (int j = 7; j >=0; j -= 4) {
	// 		out3_verilog << qH[i][j] << endl;
	// 		out2_verilog << qH[i][j - 1] << endl;
	// 		out1_verilog << qH[i][j - 2] << endl;
	// 		out0_verilog << qH[i][j - 3] << endl;
	// 	}
	// }

	// for (int i = 7; i >= 0; i -= 4) {
	// 	out3_verilog << qY[i] << endl;
	// 	out2_verilog << qY[i - 1] << endl;
	// 	out1_verilog << qY[i - 2] << endl;
	// 	out0_verilog << qY[i - 3] << endl;
	// }

	if (is_first) {
		out3_verilog << '0' << endl;
		out2_verilog << '0' << endl;
		out1_verilog << '0' << endl;
		out0_verilog << '0' << endl;
		is_first = false;
	} 

	int i = 31;
	while (i > 2) {
		out2_verilog << qH[1 + 2 * (i / 8)][i-- % 8] << endl;
		out1_verilog << qH[1 + 2 * (i / 8)][i-- % 8] << endl;
		out0_verilog << qH[1 + 2 * (i / 8)][i-- % 8] << endl;
	}
	out2_verilog << qH[1][1] << endl;
	out1_verilog << qH[1][0] << endl;
	out0_verilog << qY[7] << endl;

	out2_verilog << qY[6] << endl;
	out1_verilog << qY[5] << endl;
	out0_verilog << qY[4] << endl;

	out2_verilog << qY[3] << endl;
	out1_verilog << qY[2] << endl;
	out0_verilog << qY[1] << endl;

	out2_verilog << qY[0] << endl;
	out1_verilog << 0 << endl;
	out0_verilog << 0 << endl;
}

void MIMO_4x4::output_qRY_verilog()
{
	for (int i = 0; i < 2 * N; i++) {
		for (int j = 0; j < 2 * N; j++)
			Rout_verilog << qR[i][j] << ' ';
		Rout_verilog << endl;
	}

	for (int i = 0; i < 2 * N; i++) 
		Rout_verilog << qYtrans[i] << ' ';
	Rout_verilog << endl << endl;
		
	// for (int i = 0; i < 2 * N; i++)
	// 	Rout_verilog << xIndex[i] << ' ';
	// Rout_verilog << endl;
	
	// for (int i = 0; i < 2 * N; i++)
	// 	Rout_verilog << colpower[i] << ' ';
	// Rout_verilog << endl << endl;
}

void MIMO_4x4::output_qRY_tmn()
{
	// auto cout_buf = cout.rdbuf(out.rdbuf());

	cout << " qY[] |  qR[][]\n";
	for (int i = 0; i < 2 * N; i++) {
		printf("%5d | ", qYtrans[i]);
		// if (qYtrans[i] < 0)
		// 	printf("%4X | ", qYtrans[i]);
		// else
		// 	printf("%4X | ", qYtrans[i]);
		for (int j = 0; j < 2 * N; j++)
			printf("%5d ", qR[i][j]);
			// if (qR[i][j] < 0)
			// 	printf("%4X ", qR[i][j]);
			// else 
			// 	printf("%4X ", qR[i][j]);
		cout << endl;
	}

	cout << "  col | ";
	for (int i = 0; i < 8; i++)
		printf("%-5d ", xIndex[i]);

	cout <<endl;
	cout <<endl;
	
	// cout.rdbuf(cout_buf);
}

void MIMO_4x4::output_qRY_flp_tmn()
{
	int i, j;
	cout << " qY[]    |  qR[][]\n";
	for (i = 0; i < 2*N; i++) {
		printf("%8.4f | ", qYtrans_flp[i]);
		for (j = 0; j < 2*N; j++)
			printf("%8.4f ", qR_flp[i][j]);
		cout <<endl;
	}

	for (int i = 0; i < 8; i++)
		printf("%2d ", xIndex[i]);
	cout <<endl;
}

void MIMO_4x4::output_RY_tmn()
{
	int i, j;
	cout << "  Y[]    |   R[][]\n";
	for (i = 0; i < 2*N; i++) {
		printf("%8.4f | ", Yexp[i]);
		for (j = 0; j < 2*N; j++)
			printf("%8.4f ", R[i][j]);
		cout <<endl;
	}

	for (int i = 0; i < 8; i++)
		printf("%2d ", xIndex[i]);
	cout <<endl;
}

void MIMO_4x4::q_input_complex2int(complex<double> **Hin, complex<double> *Yin)
{
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			qH[2 * i][2 * j] = Hin[i][j].real() * (1 << FRAC);
			qH[2 * i][2 * j + 1] = -Hin[i][j].imag() * (1 << FRAC);
			qH[2 * i + 1][2 * j] = Hin[i][j].imag() * (1 << FRAC);
			qH[2 * i + 1][2 * j + 1] = Hin[i][j].real() * (1 << FRAC);			
		}
	}

	for (int i = 0; i < N; i++) {
		qY[2 * i] = Yin[i].real() * (1 << FRAC);
		qY[2 * i + 1] = Yin[i].imag() * (1 << FRAC);
	}
}

void MIMO_4x4::q_input_verifyHYX(unsigned char *x_bit)
{
	for (int i = 0; i < 2 * N; i++) 
		for (int j = 0; j < 2 * N; j++) 
			file_HYX >> qH[i][j];

	for (int i = 0; i < 2 * N; i++) 
		file_HYX >> qY[i];

	int in;	
	for (int i = 3; i >= 0; i--) {
		file_HYX >> in;
		x_bit[i] = in;
		xout_verilog << in << ' ';
	}
}

void MIMO_4x4::xRearrange()
{
	if (decomposeType == 0) 
		for (int i = 0; i < 4; i++) {
			xComplex[xIndex[i]].real(x[2 * i]);
			xComplex[xIndex[i]].imag(x[2 * i + 1]);
		}
	else if (decomposeType == 1) 
		for (int i = 0; i < 8; i++) {
			if (xIndex[i] % 2 == 0)
				xComplex[xIndex[i] / 2].real(x[i]);
			else 
				xComplex[xIndex[i] / 2].imag(x[i]);
		}
}

void MIMO_4x4::getX(complex<int>* out)
{
	for (int i = 0; i < M; i++)
		out[i] = xComplex[i];
}

void MIMO_4x4::setH(complex<double>** Hin)
{
	int i, j;
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			H[i][j] = Hin[i][j];

	return;
}

void MIMO_4x4::setY(complex<double>* Yin)
{
	for (int j = 0; j < N; j++)
		y[j] = Yin[j];
	return;
}

void MIMO_4x4::output_Matrix_CSV()
{
	auto cout_buf = cout.rdbuf(out.rdbuf());
	
	// output H
	for (int j = 0; j < 2 * M; j++) {
		for (int i = 0; i < N; i++) {
			if (j % 2 == 0)
				cout << H[i][j/2].real() << ", " << H[i][j/2].imag();
			else 
				cout << -H[i][j/2].imag() << ", " << H[i][j/2].real();	
			if (j != 2 * M - 1 || i != N - 1)
				cout << ", ";
		}
	}
	cout << endl;
	
	// Output R
	for (int j = 0; j < 2 * M; j++) {
		for (int i = 0; i < 2 * N; i++) {
			cout << R[i][j];
			if (i != 2 * N - 1 || j != 2 * M - 1)
				cout << ", ";
		}
	}
	cout << endl;

	// output y
	for (int i = 0; i < N; i++) {
		cout << y[i].real() << ", " << y[i].imag();
		if (i != N - 1)
			cout << ", ";
	}
	cout << endl;

	// output yh
	for (int i = 0; i < 2 * N; i++) {
		cout << Yexp[i];
		if (i != 2 * N - 1)
			cout << ", ";
	}
	cout << endl;

	cout.rdbuf(cout_buf);

}

void MIMO_4x4::output_X_CSV(complex<int>* xsymbol)
{
	auto cout_buf = cout.rdbuf(out.rdbuf());
	
	for (int i = 0; i < N; i++) {
		cout << xsymbol[i].real() << ", " << xsymbol[i].imag();
		if (i != N - 1)
			cout << ", ";
	}
	cout << endl;

	cout.rdbuf(cout_buf);
}

