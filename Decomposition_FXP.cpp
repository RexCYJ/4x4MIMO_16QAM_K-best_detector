
#include "MIMO_4x4.h"

void MIMO_4x4::q_decompositionFull()
{
	int power2 = 0;
	int32_t colpower[2 * M] = {0}, minpower, temp;
	int32_t mincol = 0, tmp;
	int32_t di = 0;
	int i, j, col, iter;

	decomposeType = 1;

	// expand H matrix, and y vector
	for (i = 0; i < 8; i++) {
		for (j = 0; j < 8; j++)
			qR[i][j] = qH[i][j];
		qYtrans[i] = qY[i];
	}

	// Initiate the column power
	minpower = (1 << 25);
	for (j = 0; j < 8; j++) {
		colpower[j] = 0;
		xIndex[j] = j;							// initialize xIndex
		for (i = 0; i < 8; i++) 
			colpower[j] += (((int32_t)qR[i][j] >> (FRAC - COLLEN_FRAC)) 
							* (((int32_t)qR[i][j]) >> (FRAC - COLLEN_FRAC))) >> COLLEN_FRAC;
		// printf("[%1d]=%4d ", j, colpower[j]);
		if (colpower[j] < minpower) {
			minpower = colpower[j];
			mincol = j;
		}
	}
	// cout << endl;

	// Givens Rotation
	for (j = 0; j < 7; j++) {
		// output_qRY_tmn();

		// Move the minimum power column to first column
		if (mincol != j) {
			for (i = 0; i < 8; i++) {
				temp = qR[i][j]; 					// swap 
				qR[i][j] = qR[i][mincol]; 
				qR[i][mincol] = temp;
			}
			temp = colpower[j]; 					// swap column power
			colpower[j] = colpower[mincol]; 
			colpower[mincol] = temp;
			
			tmp = xIndex[j]; 						// swap index
			xIndex[j] = xIndex[mincol]; 
			xIndex[mincol] = tmp;
		}

		// Run CORDIC for each row (Parallelizable)
		for (i = j + 1; i < 8; i++) {
			// z = 0;
			// Adjust the angel to 90 ~ -90 degree
			if (qR[j][j] < 0) {
				di = (qR[i][j] < 0)? 1: -1;
				for (col = j; col < 2 * N; col++) {
					tmp = -di * qR[i][col];
					qR[i][col] = di * qR[j][col];
					qR[j][col] = tmp;
				}
				tmp = -di * qYtrans[i];
				qYtrans[i] = di * qYtrans[j];
				qYtrans[j] = tmp;

			}

			// CORDIC core
			for (iter = 0; iter < CORDIC_ITER; iter++) {
				di = (qR[i][j] < 0)? 1: -1;
				for (col = j; col < 2 * N; col++) {
					tmp 	   = qR[j][col] - ((di * qR[i][col]) >> iter);
					qR[i][col] = ((di * qR[j][col]) >> iter) + qR[i][col];
					qR[j][col] = tmp;
				}
				tmp 	   = qYtrans[j] - ((di * qYtrans[i]) >> iter);
				qYtrans[i] = ((di * qYtrans[j]) >> iter) + qYtrans[i];
				qYtrans[j] = tmp;
			}

			// Normalize 
			for (col = j; col < 2 * N; col++) {
				qR[j][col] = ((int64_t)qR[j][col] * (int32_t)q_Const) >> FRAC;
				qR[i][col] = ((int64_t)qR[i][col] * (int32_t)q_Const) >> FRAC;
			}
			qYtrans[j] = ((int64_t)qYtrans[j] * (int32_t)q_Const) >> FRAC;
			qYtrans[i] = ((int64_t)qYtrans[i] * (int32_t)q_Const) >> FRAC;
		}	// End CORDIC

		// Recalculate the column power
		minpower = 1 << 25;
		for (int col = j + 1; col < 8; col++) {
			colpower[col] -= ((qR[j][col] >> (FRAC - COLLEN_FRAC)) 
							* (qR[j][col] >> (FRAC - COLLEN_FRAC))) >> COLLEN_FRAC;
			if (colpower[col] < minpower) {
				minpower = colpower[col];
				mincol = col;
			}
			// printf("[%1d]=%4d ", col, colpower[col]);
		}
		// cout << endl;
	}	// Go to eliminate next column
		// End Givens Rotation

	// convert to floating point
	power2 = 1 << FRAC;
	for (i = 0; i < 8; i++) {
		for (j = 0; j < 8; j++) {
			qR_flp[i][j] = R[i][j] = (double)qR[i][j] / power2;
		}
		qYtrans_flp[i] = Yexp[i] = (double)qYtrans[i] / power2;
	}
}


void MIMO_4x4::DecomposeError()
{
	double err = 0, errsum = 0;

	cout << " error: \n";
	for (int i = 0; i < 2 * N; i++) {
		err = Yexp[i] - qYtrans_flp[i];
		printf("%8.4f | ", err * err);
		for (int j = 0; j < 2 * N; j++) {
			err = R[i][j] - qR_flp[i][j];
			err = err * err;
			printf("%8.4f ", err);
			errsum += err;
		}
		cout << endl;
	}
	cout << endl;
}
