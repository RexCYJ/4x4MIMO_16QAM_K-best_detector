#include "MIMO_4x4.h"

void MIMO_4x4::decompositionHalf()
{
	double Hexp[2 * N][N];
	double z = 0;
	double colpower[4] = {0}, power, minpower, temp;
	int mincol = 0, tmp;

	decomposeType = 0;

	// expand the complex matrix H into real number matrix Hexp
	Hexp[0][0] = H[0][0].real(); Hexp[0][1] = H[0][1].real(); Hexp[0][2] = H[0][2].real(); Hexp[0][3] = H[0][3].real();
	Hexp[1][0] = H[0][0].imag(); Hexp[1][1] = H[0][1].imag(); Hexp[1][2] = H[0][2].imag(); Hexp[1][3] = H[0][3].imag();
	Hexp[2][0] = H[1][0].real(); Hexp[2][1] = H[1][1].real(); Hexp[2][2] = H[1][2].real(); Hexp[2][3] = H[1][3].real();
	Hexp[3][0] = H[1][0].imag(); Hexp[3][1] = H[1][1].imag(); Hexp[3][2] = H[1][2].imag(); Hexp[3][3] = H[1][3].imag();
	Hexp[4][0] = H[2][0].real(); Hexp[4][1] = H[2][1].real(); Hexp[4][2] = H[2][2].real(); Hexp[4][3] = H[2][3].real();
	Hexp[5][0] = H[2][0].imag(); Hexp[5][1] = H[2][1].imag(); Hexp[5][2] = H[2][2].imag(); Hexp[5][3] = H[2][3].imag();
	Hexp[6][0] = H[3][0].real(); Hexp[6][1] = H[3][1].real(); Hexp[6][2] = H[3][2].real(); Hexp[6][3] = H[3][3].real();
	Hexp[7][0] = H[3][0].imag(); Hexp[7][1] = H[3][1].imag(); Hexp[7][2] = H[3][2].imag(); Hexp[7][3] = H[3][3].imag();
	 
	Yexp[0] = y[0].real(); Yexp[1] = y[0].imag(); 
	Yexp[2] = y[1].real(); Yexp[3] = y[1].imag(); 
	Yexp[4] = y[2].real(); Yexp[5] = y[2].imag(); 
	Yexp[6] = y[3].real(); Yexp[7] = y[3].imag(); 


	// Calculate the column power
	minpower = 10000;
	for (int j = 0; j < 4; j++) {
		power = 0;
		xIndex[j] = j;		// initialize xIndex
		for (int i = 0; i < 8; i++) 
			power += Hexp[i][j] * Hexp[i][j];
		colpower[j] = power;
		if (power < minpower) {
			minpower = power;
			mincol = j;
		}
	}
	// // move the minimum power column to first column
	// if (mincol != 0) {
	// 	for (int i = 0; i < 8; i++) {
	// 		temp = Hexp[i][0]; Hexp[i][0] = Hexp[i][mincol]; Hexp[i][mincol] = temp;
	// 	}
	// 	temp = colpower[0]; colpower[0] = colpower[mincol]; colpower[mincol] = temp;
	// 	tmp = xIndex[0]; xIndex[0] = xIndex[mincol]; xIndex[mincol] = tmp;
	// }

	// COLUMN 1
	CORDIC(Hexp[0][0], Hexp[1][0], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[0], Hexp[1], 1, z);
	CORDIC(Yexp[0],    Yexp[1], z, 1);		// y Rotation-mode CORDIC
	z = 0;
	CORDIC(Hexp[2][0], Hexp[3][0], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[2], Hexp[3], 1, z);
	CORDIC(Yexp[2],    Yexp[3], z, 1);		// y Rotation-mode CORDIC
	z = 0;
	CORDIC(Hexp[4][0], Hexp[5][0], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[4], Hexp[5], 1, z);
	CORDIC(Yexp[4],    Yexp[5], z, 1);		// y Rotation-mode CORDIC
	z = 0;
	CORDIC(Hexp[6][0], Hexp[7][0], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[6], Hexp[7], 1, z);
	CORDIC(Yexp[6],    Yexp[7], z, 1);		// y Rotation-mode CORDIC
	z = 0;
	CORDIC(Hexp[0][0], Hexp[2][0], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[0], Hexp[2], 1, z);
	CORDIC_ROW(Hexp[1], Hexp[3], 1, z);
	CORDIC(Yexp[0],    Yexp[2], z, 1);		// y Rotation-mode CORDIC
	CORDIC(Yexp[1],    Yexp[3], z, 1);		// y Rotation-mode CORDIC
	z = 0;
	CORDIC(Hexp[4][0], Hexp[6][0], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[4], Hexp[6], 1, z);
	CORDIC_ROW(Hexp[5], Hexp[7], 1, z);
	CORDIC(Yexp[4],    Yexp[6], z, 1);		// y Rotation-mode CORDIC
	CORDIC(Yexp[5],    Yexp[7], z, 1);		// y Rotation-mode CORDIC
	z = 0;
	CORDIC(Hexp[0][0], Hexp[4][0], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[0], Hexp[4], 1, z);
	CORDIC_ROW(Hexp[1], Hexp[5], 1, z);
	CORDIC(Yexp[0],    Yexp[4], z, 1);		// y Rotation-mode CORDIC
	CORDIC(Yexp[1],    Yexp[5], z, 1);		// y Rotation-mode CORDIC
	z = 0;

	// minpower = 10000;
	// for (int j = 1; j < 4; j++) {
	// 	power = Hexp[0][j] * Hexp[0][j] + Hexp[1][j] * Hexp[1][j];
	// 	colpower[j] -= power;
	// 	if (colpower[j] < minpower) {
	// 		minpower = power;
	// 		mincol = j;
	// 	}
	// }
	// // move the minimum power column to first column
	// if (mincol != 1) {
	// 	for (int i = 0; i < 8; i++) {
	// 		temp = Hexp[i][1]; Hexp[i][1] = Hexp[i][mincol]; Hexp[i][mincol] = temp;
	// 	}
	// 	temp = colpower[1]; colpower[1] = colpower[mincol]; colpower[mincol] = temp;
	// 	tmp = xIndex[1]; xIndex[1] = xIndex[mincol]; xIndex[mincol] = tmp;
	// }

	// COLUMN 2
	CORDIC(Hexp[2][1], Hexp[3][1], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[2], Hexp[3], 2, z);
	CORDIC(Yexp[2],    Yexp[3], z, 1);		// y Rotation-mode CORDIC
	z = 0;
	CORDIC(Hexp[4][1], Hexp[5][1], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[4], Hexp[5], 2, z);
	CORDIC(Yexp[4],    Yexp[5], z, 1);		// y Rotation-mode CORDIC
	z = 0;
	CORDIC(Hexp[6][1], Hexp[7][1], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[6], Hexp[7], 2, z);
	CORDIC(Yexp[6],    Yexp[7], z, 1);		// y Rotation-mode CORDIC
	z = 0;
	CORDIC(Hexp[2][1], Hexp[4][1], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[2], Hexp[4], 2, z);
	CORDIC_ROW(Hexp[3], Hexp[5], 2, z);
	CORDIC(Yexp[2],    Yexp[4], z, 1);		// y Rotation-mode CORDIC
	CORDIC(Yexp[3],    Yexp[5], z, 1);		// y Rotation-mode CORDIC
	z = 0;
	CORDIC(Hexp[2][1], Hexp[6][1], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[2], Hexp[6], 2, z);
	CORDIC_ROW(Hexp[3], Hexp[7], 2, z);
	CORDIC(Yexp[2],    Yexp[6], z, 1);		// y Rotation-mode CORDIC
	CORDIC(Yexp[3],    Yexp[7], z, 1);		// y Rotation-mode CORDIC
	z = 0;

	// minpower = 10000;
	// for (int j = 2; j < 4; j++) {
	// 	power = Hexp[2][j] * Hexp[2][j] + Hexp[3][j] * Hexp[3][j];
	// 	colpower[j] -= power;
	// 	if (colpower[j] < minpower) {
	// 		minpower = power;
	// 		mincol = j;
	// 	}
	// }
	// // move the minimum power column to first column
	// if (mincol != 2) {
	// 	for (int i = 0; i < 8; i++) {
	// 		temp = Hexp[i][2]; Hexp[i][2] = Hexp[i][mincol]; Hexp[i][mincol] = temp;
	// 	}
	// 	temp = colpower[2]; colpower[2] = colpower[mincol]; colpower[mincol] = temp;
	// 	tmp = xIndex[2]; xIndex[2] = xIndex[mincol]; xIndex[mincol] = tmp;
	// }

	// COLUMN 3
	CORDIC(Hexp[4][2], Hexp[5][2], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[4], Hexp[5], 3, z);
	CORDIC(Yexp[4],    Yexp[5], z, 1);		// y Rotation-mode CORDIC
	z = 0;
	CORDIC(Hexp[6][2], Hexp[7][2], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[6], Hexp[7], 3, z);
	CORDIC(Yexp[6],    Yexp[7], z, 1);		// y Rotation-mode CORDIC
	z = 0;
	CORDIC(Hexp[4][2], Hexp[6][2], z, 0);	// Vector-mode CORDIC
	CORDIC_ROW(Hexp[4], Hexp[6], 3, z);
	CORDIC_ROW(Hexp[5], Hexp[7], 3, z);
	CORDIC(Yexp[4],    Yexp[6], z, 1);		// y Rotation-mode CORDIC
	CORDIC(Yexp[5],    Yexp[7], z, 1);		// y Rotation-mode CORDIC
	z = 0;

	// COLUMN 4
	CORDIC(Hexp[6][3], Hexp[7][3], z, 0);	// Vector-mode CORDIC
	CORDIC(Yexp[6],    Yexp[7], z, 1);		// y Rotation-mode CORDIC
	z = 0;

	R[0][0] = Hexp[0][0];	R[0][1] = 0;			R[0][2] = Hexp[0][1];	R[0][3] = -Hexp[1][1];
	R[1][0] = 0;			R[1][1] = Hexp[0][0];	R[1][2] = Hexp[1][1];	R[1][3] =  Hexp[0][1];
	R[2][0] = 0;			R[2][1] = 0;			R[2][2] = Hexp[2][1];	R[2][3] =  0;
	R[3][0] = 0;			R[3][1] = 0;			R[3][2] = 0;			R[3][3] =  Hexp[2][1];
	R[4][0] = 0;			R[4][1] = 0;			R[4][2] = 0;			R[4][3] =  0;
	R[5][0] = 0;			R[5][1] = 0;			R[5][2] = 0;			R[5][3] =  0;
	R[6][0] = 0;			R[6][1] = 0;			R[6][2] = 0;			R[6][3] =  0;
	R[7][0] = 0;			R[7][1] = 0;			R[7][2] = 0;			R[7][3] =  0;

	R[0][4] = Hexp[0][2];	R[0][5] = -Hexp[1][2];	R[0][6] = Hexp[0][3];	R[0][7] = -Hexp[1][3];
	R[1][4] = Hexp[1][2];	R[1][5] =  Hexp[0][2];	R[1][6] = Hexp[1][3];	R[1][7] =  Hexp[0][3];
	R[2][4] = Hexp[2][2];	R[2][5] = -Hexp[3][2];	R[2][6] = Hexp[2][3];	R[2][7] = -Hexp[3][3];
	R[3][4] = Hexp[3][2];	R[3][5] =  Hexp[2][2];	R[3][6] = Hexp[3][3];	R[3][7] =  Hexp[2][3];
	R[4][4] = Hexp[4][2];	R[4][5] =  0;			R[4][6] = Hexp[4][3];	R[4][7] = -Hexp[5][3];
	R[5][4] = 0;			R[5][5] =  Hexp[4][2];	R[5][6] = Hexp[5][3];	R[5][7] =  Hexp[4][3];
	R[6][4] = 0;			R[6][5] = 0;			R[6][6] = Hexp[6][3];	R[6][7] =  0;
	R[7][4] = 0;			R[7][5] = 0;			R[7][6] = 0;			R[7][7] =  Hexp[6][3];

	return;
}

void MIMO_4x4::decompositionFull()
{
	double z = 0;
	double colpower[2 * M] = {0}, minpower, temp;
	int mincol = 0, tmp;

	decomposeType = 1;

	// expand H matrix, and y vector
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			R[2 * i][2 * j] = H[i][j].real();
			R[2 * i][2 * j + 1] = -H[i][j].imag();
			R[2 * i + 1][2 * j] = H[i][j].imag();
			R[2 * i + 1][2 * j + 1] = H[i][j].real();
		}
		Yexp[2 * i] = y[i].real();
		Yexp[2 * i + 1] = y[i].imag();
	}

	// Initiate the column power
	minpower = 10000;
	for (int j = 0; j < 8; j++) {
		colpower[j] = 0;
		xIndex[j] = j;							// initialize xIndex
		for (int i = 0; i < 8; i++) 
			colpower[j] += R[i][j] * R[i][j];
		if (colpower[j] < minpower) {
			minpower = colpower[j];
			mincol = j;
		}
		// printf("[%d]=%.4f ", j, colpower[j]);
	}
	// cout << endl;

	for (int j = 0; j < 7; j++) {
		// Move the minimum power column to first column
		if (mincol != j) {
			for (int i = 0; i < 8; i++) {
				temp = R[i][j]; 					// swap 
				R[i][j] = R[i][mincol]; 
				R[i][mincol] = temp;
			}
			temp = colpower[j]; 
			colpower[j] = colpower[mincol]; 
			colpower[mincol] = temp;
			tmp = xIndex[j]; 
			xIndex[j] = xIndex[mincol]; 
			xIndex[mincol] = tmp;
		}
		for (int i = j + 1; i < 8; i++) {
			z = 0;
			CORDIC(R[j][j], R[i][j], z, 0);		// Vector-mode CORDIC
			CORDIC_ROW(R[j], R[i], j + 1, z);
			CORDIC(Yexp[j], Yexp[i], z, 1);		// y Rotation-mode CORDIC
		}
		minpower = 10000;
		for (int k = j + 1; k < 8; k++) {
			colpower[k] -= R[j][k] * R[j][k];
			if (colpower[k] < minpower) {
				minpower = colpower[k];
				mincol = k;
			}
			// printf("[%d]=%.4f ", k, colpower[k]);
		}
		// cout << endl;
	}
}

void MIMO_4x4::CORDIC_ROW(double* r1, double* r2, int startCol, double z)
{
	if (decomposeType == 0)
		for (int j = startCol; j < N; j++)
			CORDIC(r1[j], r2[j], z, 1);
	else if (decomposeType == 1)
		for (int j = startCol; j < 2 * N; j++)
			CORDIC(r1[j], r2[j], z, 1);

}

// CORDIC function, Mode: Vector:0, Rotation:1.
// vector mode will modify z
double MIMO_4x4::CORDIC(double &x, double &y, double &z, int mode)
{
	double x0, xi, y0, yi, z0, zi;	// variables
	double di;						// direction-adjusting factor
	double temp;					// temperary var
	double A = 0;					// rotation gain
	double power2 = 1;				// power of 2

	// Rotation Limit Correction
	if (x < 0 && mode == 0) {	// Vector mode rotate over 90degree, x < 0
		di = (y < 0)? 1: -1;
		x0 = xi = -di * y;
		y0 = yi = di * x;
		z0 = zi = z - di * 90;
	} else if ((z > 90 || z < -90) && mode == 1) {	// Rotation mode ratote over 90degree, |z| > 90 
		di = (z < 0)? 1: -1;
		x0 = xi = di * y;
		y0 = yi = -di * x;
		z0 = zi = z + di * 90;
	} else {
		x0 = xi = x;
		y0 = yi = y;
		z0 = zi = z;
	}

	// Iteration to rotate the vector
	for (int iter = 0; iter < CORDIC_ITER; iter++) {
		if (mode == 0)
			di = (yi < 0)? 1: -1;				// update d
		else if (mode == 1)
			di = (zi < 0)? -1: 1;

		temp = xi - yi * di * power2;					// update x, y, z
		yi = yi + xi * di * power2;
		xi = temp;
		zi = zi - di * arctan2[iter];
		power2 *= 0.5;
	}

	if (mode == 0) z = -zi;
	
	// return the new value
	x = xi / An;	// Normalize
	y = yi / An;

	return A;
}
