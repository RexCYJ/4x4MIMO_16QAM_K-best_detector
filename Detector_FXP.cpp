// Fixed-point version of MIMO detector

#include "MIMO_4x4.h"

// input: complex NxN H, complex Nx1 Y
void MIMO_4x4::q_detect()
{
	// output_qHY_tmn();
	output_qHY_verilog();
	output_qHY_verilog_batch();

	q_decompositionFull();
	// q_decompositionFull_test();
	
	// output_qRY_tmn();
	// output_RY_tmn();
	// output_qRY_flp_tmn();
	// output_qRY_verilog();

	// setH(Hin);
	// setY(Yin);
	// decompositionFull();
	// output_RY_tmn();

	// Kbest_optiK();
	// q_Kbest_optiK();
	q_Kbest4();
	// Kbest4();
	// Kbest3();

	xRearrange();
}

// Quantilized optimized K-best algorithm
void MIMO_4x4::q_Kbest_optiK()
{
	const int K = 3;					// the number of candidates
	vector<int> path(8, 0);				// selected X
	vector<vector<int>> temp(K, path);	// store the next level candidate
	vector<vector<int>> Kpath(K, path);	// the K solutions
	int32_t T[K] = {0}, tempT[K] = {0};	// the path error accumulator
	int32_t CurLevelErr[K] = {0};
	pair<int, int32_t> kcand[K];		// store the possible candidate
										// (index of QAM16_val[], error)
	int ZZstep[K];

	int32_t err = 0, minerr = 0, accumErr = 0, hx = 0;
	int32_t Err[3];
	int i, j, k, bestpath = 0;
	int32_t MAX = (1 << 20) - 1;

	for (i = 0; i < K; i++)
		kcand[i] = make_pair(0, MAX);	// initialize the candidate container

	// run through all possible combination in the lowest two level
	// search the best path to the third lowest level
	for (i = 0; i < 4; i++) {
		err = (int32_t)qYtrans[7] - (((int32_t)qR[7][7] * QAM16_q_normval[i]) >> FRAC);
		accumErr = (err * err) >> FRAC;
		// printf("L1: %6d (%d)\n", accumErr, i);

		for (j = 0; j < 4; j++) {
			err = (int32_t)qYtrans[6] - (((int32_t)qR[6][6] * QAM16_q_normval[j]) >> FRAC)
						- (((int32_t)qR[6][7] * QAM16_q_normval[i]) >> FRAC);
			accumErr += Err[1] = (err * err) >> FRAC;
			// printf("L2: %6d (%d, %d)\n", accumErr, i, j);

			for (int l = 0; l < 4; l++) {
				err = (int32_t)qYtrans[5] - (((int32_t)qR[5][5] * QAM16_q_normval[l]) >> FRAC)
							- (((int32_t)qR[5][6] * QAM16_q_normval[j]) >> FRAC)
							- (((int32_t)qR[5][7] * QAM16_q_normval[i])	>> FRAC);
				accumErr += Err[2] = (err * err) >> FRAC;
				// printf("L3: %6d (%d, %d, %d)\n", accumErr, i, j, l);

				// update the best one
				if (accumErr < kcand[0].second) {
					kcand[2] = kcand[1];
					kcand[1] = kcand[0];
					kcand[0] = make_pair(l, accumErr);
					Kpath[2].assign(Kpath[1].begin(), Kpath[1].end());
					Kpath[1].assign(Kpath[0].begin(), Kpath[0].end());
					Kpath[0][7] = i; Kpath[0][6] = j; Kpath[0][5] = l;
					
					// printf("1. %6d", kcand[0].second);
					// printf("   %3d -> %3d -> %3d\n", Kpath[0][7], Kpath[0][6], Kpath[0][5]);
					
				} else if (accumErr < kcand[1].second) {
					kcand[2] = kcand[1];
					kcand[1] = make_pair(l, accumErr);
					Kpath[2].assign(Kpath[1].begin(), Kpath[1].end());
					Kpath[1][7] = i; Kpath[1][6] = j; Kpath[1][5] = l;

					// printf("2. %6d", kcand[1].second);
					// printf("   %3d -> %3d -> %3d\n", Kpath[1][7], Kpath[1][6], Kpath[1][5]);

				} else if (accumErr < kcand[2].second) {
					kcand[2] = make_pair(l, accumErr);
					Kpath[2][7] = i; Kpath[2][6] = j; Kpath[2][5] = l;

					// printf("3. %6d", kcand[2].second);
					// printf("   %3d -> %3d -> %3d\n", Kpath[2][7], Kpath[2][6], Kpath[2][5]);

				}
				accumErr -= Err[2];
			}
			accumErr -= Err[1];
		}
	}

	T[0] = kcand[0].second;
	T[1] = kcand[1].second;
	T[2] = kcand[2].second;

	// Finished Two Levels

	for (int l = 4; l >= 0; l--) {				// search each level
		// Find the first best child of each path

		for (k = 0; k < K; k++) {
			hx = 0;
			// update the accumulated error
			for (j = 7; j > l; j--) 
				hx += ((int32_t)qR[l][j] * QAM16_q_normval[Kpath.at(k).at(j)]) >> FRAC;
			CurLevelErr[k] = qYtrans[l] - hx;					// the rest y value
			// find the best first child
			kcand[k].second = MAX;
			for (i = 0; i < 4; i++) {
				err = CurLevelErr[k] - (((int32_t)qR[l][l] * QAM16_q_normval[i]) >> FRAC);
				accumErr = T[k] + ((err * err) >> FRAC);
				if (kcand[k].second > accumErr) 
					kcand[k] = make_pair(i, accumErr);
			}
			Kpath.at(k).at(l) = kcand[k].first;

			if (kcand[k].first == 3)
				ZZstep[k] = -1;
			else if (kcand[k].first == 0)
				ZZstep[k] = 1;
			else if (!((err & 0x8000000) ^ ((int32_t)qR[l][l] & 0x8000000)))
				ZZstep[k] = 1;
			else if ((err & 0x80000000) ^ ((int32_t)qR[l][l] & 0x8000000))
				ZZstep[k] = -1;
		}
		
		// find the best path
		for (k = 0; k < K; k++) {
			minerr = MAX;
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

			err = CurLevelErr[bestpath] - (((int32_t)qR[l][l] * QAM16_q_normval[kcand[bestpath].first]) >> FRAC);
			accumErr = T[bestpath] + ((err * err) >> FRAC);
			kcand[bestpath].second = accumErr;
			Kpath[bestpath][l] = kcand[bestpath].first;
		}

		Kpath.at(0).assign(temp.at(0).begin(), temp.at(0).end());	// the best path
		Kpath.at(1).assign(temp.at(1).begin(), temp.at(1).end());
		Kpath.at(2).assign(temp.at(2).begin(), temp.at(2).end());

		T[0] = tempT[0]; T[1] = tempT[1]; T[2] = tempT[2];
		// printf("%6d %6d %6d\n", T[0], T[1], T[2]);
	}
	
	for (int j = 0; j < 8; j++)
		x[j] = QAM16_val[Kpath[0][j]];

}

// Quantilized K-best algorithm with K = 4
void MIMO_4x4::q_Kbest4()
{
	const int K = 4;					// the number of candidates
	vector<int> path(8, 0);				// selected X
	vector<vector<int>> temp(K, path);	// store the next level candidate
	vector<vector<int>> Kpath(K, path);	// the K solutions
	int32_t T[K] = {0}, tempT[K] = {0};	// the path error accumulator
	int32_t CurLevelErr[K] = {0};
	pair<int, int32_t> kcand[K];			// store the possible candidate
										// (index of QAM16_val[], error)
	pair<int, int> order[K];			// <type of enum path, current index>
	// int ZZstep[K];

	int32_t err = 0, minerr = 0, accumErr = 0,hx = 0;
	int32_t MAX = 1 << 25;
	int i, j, k, bestpath = 0;

	// initialize the level 7 of the i-th Kpath
	for (i = 0; i < 4; i++) {
		err = (qYtrans[7] - (((int64_t)qR[7][7] * QAM16_q_normval[i]) >> FRAC)) >> (FRAC - ERR_FRAC);
		err = (err * err) >> ERR_FRAC;								// square
		T[i] = err;
		Kpath.at(i).at(7) = i;
	}

	for (int l = 6; l >= 0; l--) {				// search each level

		// Find the first best child of each path
		for (k = 0; k < K; k++) {
			hx = 0;
			// update the accumulated error
			for (j = 7; j > l; j--) 
				hx += ((int64_t)qR[l][j] * QAM16_q_normval[Kpath.at(k).at(j)]) >> FRAC;
			CurLevelErr[k] = (int64_t)qYtrans[l] - hx;					// the rest y value

			// find the best first child
			kcand[k].second = MAX;
			for (i = 0; i < 4; i++) {
				err = (CurLevelErr[k] - (((int64_t)qR[l][l] * QAM16_q_normval[i]) >> FRAC)) >> (FRAC - ERR_FRAC);
				accumErr = T[k] + ((err * (int64_t)err) >> ERR_FRAC);

				if (kcand[k].second > accumErr) {
					kcand[k] = make_pair(i, accumErr);
					Kpath.at(k).at(l) = kcand[k].first;
					
					order[k].first = kcand[k].first * 2;
					// if err and R[l][l] are same sign, x > kcand[k].first so order + 1;
					if (!((err & 0x80000000) ^ (qR[l][l] & 0x80000000))) order[k].first++;	
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
				kcand[bestpath].first =zigzag_path[order[bestpath].first][order[bestpath].second];
				order[bestpath].second++;

				// calculate the error of new path
				err = (CurLevelErr[bestpath] - (((int64_t)qR[l][l] * QAM16_q_normval[kcand[bestpath].first]) >> FRAC)) >> (FRAC - ERR_FRAC);
				accumErr = T[bestpath] + ((err * (int64_t)err) >> ERR_FRAC);
				kcand[bestpath].second = accumErr;
				Kpath[bestpath][l] = kcand[bestpath].first;
			}
		}	// End the level searching

		Kpath.at(0).assign(temp.at(0).begin(), temp.at(0).end());	// the best path
		Kpath.at(1).assign(temp.at(1).begin(), temp.at(1).end());
		Kpath.at(2).assign(temp.at(2).begin(), temp.at(2).end());
		Kpath.at(3).assign(temp.at(3).begin(), temp.at(3).end());
		T[0] = tempT[0]; T[1] = tempT[1]; T[2] = tempT[2]; T[3] = tempT[3];
	}	// End K-best

	// mapping
	for (int j = 0; j < 8; j++)
		x[j] = QAM16_val[Kpath[0][j]];

}