#include "MIMO_4x4.h"

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

void MIMO_4x4::detect()
{
	// decompositionHalf();
	decompositionFull();

	// outputMatrix_CSV();
	
	// Kbest();
	Kbest3();
	// Kbest_optiK();
	// ML();

	xRearrange();
}

void MIMO_4x4::Kbest_optiK()
{
	const int K = 3;					// the number of candidates
	vector<int> path(8, 0);				// selected X
	vector<vector<int>> temp(K, path);	// store the next level candidate
	vector<vector<int>> Kpath(K, path);	// the K solutions
	double T[K] = {0}, tempT[K] = {0};	// the path error accumulator
	double CurLevelErr[K] = {0};
	pair<int, double> kcand[K];			// store the possible candidate
										// (index of QAM16_val[], error)
	int ZZstep[K];

	double err = 0, minerr = 0, accumErr = 0,hx = 0;
	double Err[3];
	int i, j, k, bestpath = 0;

	for (i = 0; i < K; i++)
		kcand[i] = make_pair(0, 10000);	// initialize the candidate container

	// run through all possible combination in the lowest two level
	// search the best path to the third lowest level
	for (i = 0; i < 4; i++) {
		err = Yexp[7] - (R[7][7] * QAM16_normval[i]);
		accumErr = err * err;
		for (j = 0; j < 4; j++) {
			err = Yexp[6] - (R[6][6] * QAM16_normval[j] + R[6][7] * QAM16_normval[i]);
			accumErr += Err[1] = err * err;
			for (int l = 0; l < 4; l++) {
				err = Yexp[5] - (R[5][5] * QAM16_normval[l] + R[5][6] * QAM16_normval[j] + R[5][7] * QAM16_normval[i]);
				accumErr += Err[2] = err * err;
				if (accumErr < kcand[0].second) {
					kcand[2] = kcand[1];
					kcand[1] = kcand[0];
					kcand[0] = make_pair(l, accumErr);
					Kpath[2].assign(Kpath[1].begin(), Kpath[1].end());
					Kpath[1].assign(Kpath[0].begin(), Kpath[0].end());
					Kpath[0][7] = i; Kpath[0][6] = j; Kpath[0][5] = l;
				} else if (accumErr < kcand[1].second) {
					kcand[2] = kcand[1];
					kcand[1] = make_pair(l, accumErr);
					Kpath[2].assign(Kpath[1].begin(), Kpath[1].end());
					Kpath[1][7] = i; Kpath[1][6] = j; Kpath[1][5] = l;
				} else if (accumErr < kcand[2].second) {
					kcand[2] = make_pair(l, accumErr);
					Kpath[2][7] = i; Kpath[2][6] = j; Kpath[2][5] = l;
				}
				accumErr -= Err[2];
			}
			accumErr -= Err[1];
		}
	}

	// Kpath.at(0).at(5) = kcand[0].first;			// lowest error path
	// Kpath.at(1).at(5) = kcand[1].first;
	// Kpath.at(2).at(5) = kcand[2].first;
	T[0] = kcand[0].second;
	T[1] = kcand[1].second;
	T[2] = kcand[2].second;

	for (int l = 4; l >= 0; l--) {				// search each level
		// Find the first best child of each path

		// for (k = 0; k < K; k++) {
		// 	for (i = 7; i > l; i--)
		// 		cout << Kpath[k][i] << " -> ";
		// 	cout << endl;
		// }

		for (k = 0; k < K; k++) {
			hx = 0;
			// update the accumulated error
			for (j = 7; j > l; j--) 
				hx += R[l][j] * QAM16_normval[Kpath.at(k).at(j)];
			CurLevelErr[k] = Yexp[l] - hx;					// the rest y value
			// find the best first child
			kcand[k].second = 100000;
			for (i = 0; i < 4; i++) {
				err = CurLevelErr[k] - R[l][l] * QAM16_normval[i];
				accumErr = T[k] + err * err;
				if (kcand[k].second > accumErr) 
					kcand[k] = make_pair(i, accumErr);
			}
			Kpath.at(k).at(l) = kcand[k].first;
			if (kcand[k].first == 3)
				ZZstep[k] = -1;
			else if (kcand[k].first == 0)
				ZZstep[k] = 1;
			else if (err * R[l][l] > 0)
				ZZstep[k] = 1;
			else if (err * R[l][l] < 0)
				ZZstep[k] = -1;
		}
		
		// find the best path
		for (k = 0; k < K; k++) {
			minerr = 10000;
			bestpath = 0;
			for (j = 0; j < K; j++)					// Find the current best path
				if (kcand[j].second < minerr) {
					bestpath = j;
					minerr = kcand[j].second;
				}
			// store the selected path to temp
			temp.at(k).assign(Kpath[bestpath].begin(), Kpath[bestpath].end());
			tempT[k] = kcand[bestpath].second;

			// replace the selected path with its sibling, select the neighbor child
			// Zig-Zag searching
			kcand[bestpath].first += ZZstep[bestpath];
			if (kcand[bestpath].first < 0)
				kcand[bestpath].first = 2;
			else if (kcand[bestpath].first > 3)
				kcand[bestpath].first = 1;
			if (ZZstep[bestpath] > 0)
				ZZstep[bestpath] = -(ZZstep[bestpath] + 1);
			else if (ZZstep[bestpath] < 0)
				ZZstep[bestpath] = -(ZZstep[bestpath] - 1);

			err = CurLevelErr[bestpath] - R[l][l] * QAM16_normval[kcand[bestpath].first];
			accumErr = T[bestpath] + err * err;
			kcand[bestpath].second = accumErr;
			Kpath[bestpath][l] = kcand[bestpath].first;
		}

		Kpath.at(0).assign(temp.at(0).begin(), temp.at(0).end());	// the best path
		Kpath.at(1).assign(temp.at(1).begin(), temp.at(1).end());
		Kpath.at(2).assign(temp.at(2).begin(), temp.at(2).end());
		T[0] = tempT[0]; T[1] = tempT[1]; T[2] = tempT[2];
	}
	
	for (int j = 0; j < 8; j++)
		x[j] = QAM16_val[Kpath[0][j]];
}

void MIMO_4x4::Kbest3()
{
	const int K = 3;					// the number of candidates
	vector<int> path(8, 0);				// selected X
	vector<vector<int>> temp(K, path);	// store the next level candidate
	vector<vector<int>> Kpath(K, path);	// the K solutions
	double T[K] = {0}, tempT[K] = {0};	// the path error accumulator
	double CurLevelErr[K] = {0};
	pair<int, double> kcand[K];			// store the possible candidate
										// (index of QAM16_val[], error)
	int ZZstep[K];

	double err = 0, minerr = 0, accumErr = 0,hx = 0;
	int i, j, k, bestpath = 0;

	for (i = 0; i < K; i++)
		kcand[i] = make_pair(0, 100);	// initialize the candidate container

	// initialize the first level of Kpath
	for (i = 0; i < 4; i++) {
		err = Yexp[7] - (R[7][7] * QAM16_normval[i]);
		err = err * err;								// square
		if (err < kcand[0].second) {
			kcand[2] = kcand[1];
			kcand[1] = kcand[0];
			kcand[0] = make_pair(i, err);
		} else if (err < kcand[1].second) {
			kcand[2] = kcand[1];
			kcand[1] = make_pair(i, err);
		} else if (err < kcand[2].second) 
			kcand[2] = make_pair(i, err);
	}

	Kpath.at(0).at(7) = kcand[0].first;			// lowest error path
	Kpath.at(1).at(7) = kcand[1].first;
	Kpath.at(2).at(7) = kcand[2].first;
	T[0] = kcand[0].second;
	T[1] = kcand[1].second;
	T[2] = kcand[2].second;

	for (int l = 6; l >= 0; l--) {				// search each level
		// Find the first best child of each path

		// for (k = 0; k < K; k++) {
		// 	for (i = 7; i > l; i--)
		// 		cout << Kpath[k][i] << " -> ";
		// 	cout << endl;
		// }

		for (k = 0; k < K; k++) {
			hx = 0;
			// update the accumulated error
			for (j = 7; j > l; j--) 
				hx += R[l][j] * QAM16_normval[Kpath.at(k).at(j)];
			CurLevelErr[k] = Yexp[l] - hx;					// the rest y value
			// find the best first child
			kcand[k].second = 100000;
			for (i = 0; i < 4; i++) {
				err = CurLevelErr[k] - R[l][l] * QAM16_normval[i];
				accumErr = T[k] + err * err;
				if (kcand[k].second > accumErr) 
					kcand[k] = make_pair(i, accumErr);
			}
			Kpath.at(k).at(l) = kcand[k].first;
			if (kcand[k].first == 3)
				ZZstep[k] = -1;
			else if (kcand[k].first == 0)
				ZZstep[k] = 1;
			else if (err * R[l][l] > 0)
				ZZstep[k] = 1;
			else if (err * R[l][l] < 0)
				ZZstep[k] = -1;
		}
		
		// find the best path
		for (k = 0; k < K; k++) {
			minerr = 10000;
			bestpath = 0;
			for (j = 0; j < K; j++)					// Find the current best path
				if (kcand[j].second < minerr) {
					bestpath = j;
					minerr = kcand[j].second;
				}
			// store the selected path to temp
			temp.at(k).assign(Kpath[bestpath].begin(), Kpath[bestpath].end());
			tempT[k] = kcand[bestpath].second;

			// replace the selected path with its sibling, select the neighbor child
			// Zig-Zag searching
			kcand[bestpath].first += ZZstep[bestpath];
			if (kcand[bestpath].first < 0)
				kcand[bestpath].first = 2;
			else if (kcand[bestpath].first > 3)
				kcand[bestpath].first = 1;
			if (ZZstep[bestpath] > 0)
				ZZstep[bestpath] = -(ZZstep[bestpath] + 1);
			else if (ZZstep[bestpath] < 0)
				ZZstep[bestpath] = -(ZZstep[bestpath] - 1);

			err = CurLevelErr[bestpath] - R[l][l] * QAM16_normval[kcand[bestpath].first];
			accumErr = T[bestpath] + err * err;
			kcand[bestpath].second = accumErr;
			Kpath[bestpath][l] = kcand[bestpath].first;
		}

		Kpath.at(0).assign(temp.at(0).begin(), temp.at(0).end());	// the best path
		Kpath.at(1).assign(temp.at(1).begin(), temp.at(1).end());
		Kpath.at(2).assign(temp.at(2).begin(), temp.at(2).end());
		T[0] = tempT[0]; T[1] = tempT[1]; T[2] = tempT[2];
	}
	
	// for (k = 0; k < K; k++) {
	// 	for (i = 7; i >= 0; i--)
	// 		cout << Kpath[k][i] << " -> ";
	// 	cout << endl;
	// }

	for (int j = 0; j < 8; j++)
		x[j] = QAM16_val[Kpath[0][j]];
}

void MIMO_4x4::Kbest()
{
	const int K = 2;					// the number of candidates
	vector<int> path(8, 0);				// selected X
	vector<vector<int>> temp(K, path);	// store the next level candidate
	vector<vector<int>> Kpath(K, path);	// the K solutions
	double T[K] = {0}, tempT[K] = {0};	// the path error accumulator
	double CurLevelErr[K] = {0};
	pair<int, double> kcand[K];			// store the possible candidate
										// (index of QAM16_val[], error)

	double err = 0, minerr = 0, accumErr = 0,hx = 0;
	int i, j, k, bestpath = 0;

	for (i = 0; i < K; i++)
		kcand[i] = make_pair(0, 100);	// initialize the candidate container

	// initialize the first level of Kpath
	for (i = 0; i < 4; i++) {
		err = Yexp[7] - (R[7][7] * QAM16_normval[i]);
		err = err * err;								// square
		if (err < kcand[0].second) {
			kcand[1] = kcand[0];
			kcand[0] = make_pair(i, err);
		} else if (err < kcand[1].second) kcand[1] = make_pair(i, err);
	}

	Kpath.at(0).at(7) = kcand[0].first;			// lowest error path
	Kpath.at(1).at(7) = kcand[1].first;
	T[0] = kcand[0].second;
	T[1] = kcand[1].second;

	for (int l = 6; l >= 0; l--) {				// search each level
		// Find the first best child of each path
		for (k = 0; k < K; k++) {
			hx = 0;
			// update the accumulated error
			for (j = 7; j > l; j--) 
				hx += R[l][j] * QAM16_normval[Kpath.at(k).at(j)];
			CurLevelErr[k] = Yexp[l] - hx;					// the rest y value
			// find the best first child
			kcand[k].second = 100000;
			for (i = 0; i < 4; i++) {
				err = CurLevelErr[k] - R[l][l] * QAM16_normval[i];
				accumErr = T[k] + err * err;
				if (kcand[k].second > accumErr) 
					kcand[k] = make_pair(i, accumErr);
			}
			Kpath.at(k).at(l) = kcand[k].first;
		}
		
		// find the best path
		for (k = 0; k < K; k++) {
			minerr = 10000;
			bestpath = 0;
			for (j = 0; j < K; j++)					// Find the current best path
				if (kcand[j].second < minerr) {
					bestpath = j;
					minerr = kcand[j].second;
				}
			// store the selected path to temp
			temp.at(k).assign(Kpath[bestpath].begin(), Kpath[bestpath].end());
			tempT[k] = kcand[bestpath].second;

			// replace the selected path with its sibling, select the neighbor child
			// Zig-Zag searching
			if (kcand[bestpath].first == 3)
				kcand[bestpath].first = 2;
			else if (kcand[bestpath].first == 0)
				kcand[bestpath].first = 1;
			else {
				err = CurLevelErr[bestpath] - R[l][l] * QAM16_normval[kcand[bestpath].first];
				if (err * R[l][l] > 0)
					kcand[bestpath].first++;
				else 
					kcand[bestpath].first--;
			}
			// kcand[bestpath].first = (kcand[bestpath].first + 1) % 4;

			err = CurLevelErr[bestpath] - R[l][l] * QAM16_normval[kcand[bestpath].first];
			accumErr = T[bestpath] + err * err;
			kcand[bestpath].second = accumErr;
			Kpath[bestpath][l] = kcand[bestpath].first;
		}
		Kpath.at(0).assign(temp.at(0).begin(), temp.at(0).end());	// the best path
		Kpath.at(1).assign(temp.at(1).begin(), temp.at(1).end());
		T[0] = tempT[0]; T[1] = tempT[1];
	}

	for (int j = 0; j < 8; j++)
		x[j] = QAM16_val[Kpath[0][j]];
}

void MIMO_4x4::ML()
{
	int xminerr[8] = {0};
	int i0, i1, i2, i3, j;
	double minerr = 10000, err, errsum = 0;
	static bool is_firsttime = 1;

	if (is_firsttime) {
		is_firsttime = 0;
		for (int j = 0; j < 4; j++)
			xIndex[j] = j;		// initialize xIndex
	}

	for (int i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			Hr[2 * i][2 * j] 			= H[i][j].real();  
			Hr[2 * i][2 * j + 1] 		= -H[i][j].imag();	
			Hr[2 * i + 1][2 * j] 		= H[i][j].imag();
			Hr[2 * i + 1][2 * j + 1]	= H[i][j].real();
		}
	}

	Yr[0] = y[0].real();	Yr[4] = y[2].real();
	Yr[1] = y[0].imag();	Yr[5] = y[2].imag();
	Yr[2] = y[1].real();	Yr[6] = y[3].real();
	Yr[3] = y[1].imag();	Yr[7] = y[3].imag();

	for (i3 = 0; i3 < 16; i3++) {
		for (i2 = 0; i2 < 16; i2++) {
			for (i1 = 0; i1 < 16; i1++) {
				for (i0 = 0; i0 < 16; i0++) {
					errsum = 0;
					for (j = 0; j < 2 * N; j++) {
						err = Yr[j] - (Hr[j][0] * QAM16_val[i0 % 4] + Hr[j][1] * QAM16_val[i0 / 4] 
									+  Hr[j][2] * QAM16_val[i1 % 4] + Hr[j][3] * QAM16_val[i1 / 4]
									+  Hr[j][4] * QAM16_val[i2 % 4] + Hr[j][5] * QAM16_val[i2 / 4]
									+  Hr[j][6] * QAM16_val[i3 % 4] + Hr[j][7] * QAM16_val[i3 / 4]) / KMOD;
						errsum += err * err;
						if (errsum > minerr) 
							break;
					}
					if (errsum < minerr) {
						minerr = errsum;
						xminerr[0] = QAM16_val[i0 % 4];
						xminerr[1] = QAM16_val[i0 / 4];
						xminerr[2] = QAM16_val[i1 % 4];
						xminerr[3] = QAM16_val[i1 / 4];
						xminerr[4] = QAM16_val[i2 % 4];
						xminerr[5] = QAM16_val[i2 / 4];
						xminerr[6] = QAM16_val[i3 % 4];
						xminerr[7] = QAM16_val[i3 / 4];
					}
				}
			}
		}
	}

	xComplex[0].real(xminerr[0]);	xComplex[0].imag(xminerr[1]);
	xComplex[1].real(xminerr[2]);	xComplex[1].imag(xminerr[3]);	
	xComplex[2].real(xminerr[4]);	xComplex[2].imag(xminerr[5]);
	xComplex[3].real(xminerr[6]);	xComplex[3].imag(xminerr[7]);
}

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
	// double Hexp[2 * N][2 * M];
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
	}

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
		}
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

void MIMO_4x4::outputMatrix_CSV()
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

void MIMO_4x4::output_R_tmn()
{
	int i, j;
	for (i = 0; i < 2*N; i++) {
		for (j = 0; j < 2*N; j++) {
			cout << R[i][j] << ", ";
		}
		cout <<endl;
	}
}

