#include "MIMO_4x4.h"

void MIMO_4x4::detect()
{
	// decompositionHalf();
	decompositionFull();

	// output_Matrix_CSV();
	
	// Kbest();
	// Kbest3();
	Kbest4();
	// Kbest_optiK();
	// ML();

	xRearrange();
}

void MIMO_4x4::Kbest_optiK()
{
	const int K_MAX = 3;
	int K = 3;					// the number of candidates
	vector<int> path(8, 0);				// selected X
	vector<vector<int>> temp(K_MAX, path);	// store the next level candidate
	vector<vector<int>> Kpath(K_MAX, path);	// the K solutions
	double T[K_MAX] = {0}, tempT[K_MAX] = {0};	// the path error accumulator
	double CurLevelErr[K_MAX] = {0};
	pair<int, double> kcand[K_MAX];			// store the possible candidate
										// (index of QAM16_val[], error)
	// int ZZstep[K];
	pair<int, int> order[K_MAX];			// <type of enum path, current index>

	double err = 0, minerr = 0, accumErr = 0,hx = 0;
	double Err[3];
	int i, j, k, bestpath = 0;

	for (i = 0; i < K; i++)
		kcand[i] = make_pair(0, 10000);	// initialize the candidate container

	// run through all possible combination in the lowest two level
	// search the best path to the third lowest level
	for (i = 0; i < 4; i++) {		// L7
		err = Yexp[7] - (R[7][7] * QAM16_normval[i]);
		accumErr = err * err;
		for (j = 0; j < 4; j++) {		// L6
			err = Yexp[6] - (R[6][6] * QAM16_normval[j] + R[6][7] * QAM16_normval[i]);
			accumErr += Err[1] = err * err;

			for (int l = 0; l < 4; l++) {		// L5
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

	T[0] = kcand[0].second;
	T[1] = kcand[1].second;
	T[2] = kcand[2].second;

	for (int l = 4; l >= 0; l--) {				// search each level
		// Find the best child of each path
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
				if (kcand[k].second > accumErr) {
					kcand[k] = make_pair(i, accumErr);
					Kpath.at(k).at(l) = kcand[k].first;
					order[k].first = kcand[k].first * 2;
					if (err * R[l][l] > 0) order[k].first++;	// if err and R[l][l] are same sign, x > kcand[k].first so order + 1;
					order[k].second = 0;	// start from first

				}
			}
		}
		
		// find the best path
		for (k = 0; k < K; k++) {
			minerr = 10000;
			bestpath = 0;
			
			// Find the current best path
			for (j = 0; j < K; j++)
				if (kcand[j].second < minerr) {
					bestpath = j;
					minerr = kcand[j].second;
				}
			// store the selected path to temp
			temp.at(k).assign(Kpath[bestpath].begin(), Kpath[bestpath].end());
			tempT[k] = kcand[bestpath].second;

			// replace the selected path with its sibling, select the neighbor child
			// Zig-Zag searching
			if (k < 2) {
				kcand[bestpath].first =zigzag_path[order[bestpath].first][order[bestpath].second];
				order[bestpath].second++;

				// calculate the error of new path
				err = CurLevelErr[bestpath] - R[l][l] * QAM16_normval[kcand[bestpath].first];
				accumErr = T[bestpath] + err * err;
				kcand[bestpath].second = accumErr;
				Kpath[bestpath][l] = kcand[bestpath].first;
			}
		}

		Kpath.at(0).assign(temp.at(0).begin(), temp.at(0).end());	// the best path
		Kpath.at(1).assign(temp.at(1).begin(), temp.at(1).end());
		Kpath.at(2).assign(temp.at(2).begin(), temp.at(2).end());
		T[0] = tempT[0]; T[1] = tempT[1]; T[2] = tempT[2];
	}
	
	// mapping
	for (int j = 0; j < 8; j++)
		x[j] = QAM16_val[Kpath[0][j]];
}

void MIMO_4x4::Kbest4()
{
	const int K = 4;					// the number of candidates
	vector<int> path(8, 0);				// selected X
	vector<vector<int>> temp(K, path);	// store the next level candidate
	vector<vector<int>> Kpath(K, path);	// the K solutions
	double T[K] = {0}, tempT[K] = {0};	// the path error accumulator
	double CurLevelErr[K] = {0};
	pair<int, double> kcand[K];			// store the possible candidate
										// (index of QAM16_val[], error)
	pair<int, int> order[K];			// <type of enum path, current index>
	// int ZZstep[K];

	double err = 0, minerr = 0, accumErr = 0,hx = 0;
	int i, j, k, bestpath = 0;

	for (i = 0; i < K; i++)
		kcand[i] = make_pair(0, 100);	// initialize the candidate container

	// initialize the level 7 of the i-th Kpath
	for (i = 0; i < 4; i++) {
		err = Yexp[7] - (R[7][7] * QAM16_normval[i]);
		err = err * err;								// square
		T[i] = err;
		Kpath.at(i).at(7) = i;
	}

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
				if (kcand[k].second > accumErr) {
					kcand[k] = make_pair(i, accumErr);
					Kpath.at(k).at(l) = kcand[k].first;
					
					order[k].first = kcand[k].first * 2;
					if (err * R[l][l] > 0) order[k].first++;	// if err and R[l][l] are same sign, x > kcand[k].first so order + 1;
					order[k].second = 0;	// start from first
				}
			}
		}	// End finding best child
		
		// find the best path
		for (k = 0; k < K; k++) {
			minerr = 10000;
			bestpath = 0;

			// Find the current best path
			for (j = 0; j < K; j++)					
				if (kcand[j].second < minerr) {
					bestpath = j;
					minerr = kcand[j].second;
				}
			
			// store the selected path to temp
			temp.at(k).assign(Kpath[bestpath].begin(), Kpath[bestpath].end());
			tempT[k] = kcand[bestpath].second;

			// replace the selected path with its sibling, select the neighbor child
			// Zig-Zag searching
			if (k < 3) {
				kcand[bestpath].first = zigzag_path[order[bestpath].first][order[bestpath].second];
				order[bestpath].second++;

				// calculate the error of new path
				err = CurLevelErr[bestpath] - R[l][l] * QAM16_normval[kcand[bestpath].first];
				accumErr = T[bestpath] + err * err;
				kcand[bestpath].second = accumErr;
				Kpath[bestpath][l] = kcand[bestpath].first;
			}
		}

		Kpath.at(0).assign(temp.at(0).begin(), temp.at(0).end());	// the best path
		Kpath.at(1).assign(temp.at(1).begin(), temp.at(1).end());
		Kpath.at(2).assign(temp.at(2).begin(), temp.at(2).end());
		Kpath.at(3).assign(temp.at(3).begin(), temp.at(3).end());
		T[0] = tempT[0]; T[1] = tempT[1]; T[2] = tempT[2]; T[3] = tempT[3];
	}

	// mapping
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
				if (kcand[k].second > accumErr) {
					kcand[k] = make_pair(i, accumErr);
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
			}
		}
		
		// find the best path
		for (k = 0; k < K; k++) {
			minerr = 10000;
			bestpath = 0;

			// Find the current best path
			for (j = 0; j < K; j++)					
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
			
			// Update zigzag step
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
