#include "Global_aligner.h"

Global_aligner::Global_aligner(Sequence a, Sequence b) {
	set_sequence1(a.toString());
	set_sequence2(b.toString());
	set_match_score(2);
	set_mismatch_score(-1);
	set_cap_penalty(-3); // -2
	initialize_matrix(a.toString(), b.toString());
}

Global_aligner::Global_aligner(string a, string b) {
	set_sequence1(a);
	set_sequence2(b);
	set_match_score(2);
	set_mismatch_score(-1);
	set_cap_penalty(-3);
	initialize_matrix(a, b);
}

/**
  * Aligns the sequence globally.
  * 
  * assume: matrix in class Sequence_aligner is initialised and filled with
  *  start values.
  * result: contains two sequences that are in alignment.
  * complexity: walks through the matrix. Time complexity is O(N*M).
  */
vector<Sequence>* Global_aligner::traceback() {
	string result_sequence1;
	string result_sequence2;

	// The algorithm starts from the last element --> both sequences globally aligned.
	int pos = x_koko * y_koko - 1;

	if (verbalize) {cout << "traceback started, please wait. Pos: " << pos << endl;}

	while (pos > 0) {
		// if match, then go there.jos match niin mene sinne
		// fetch the x and y-component:

		/*
			Path is checked from the traceback matrix -list.
			When going left diagonal (ldiag), both sequence letter is added
			left --> cap seq2, seq1 letter.
			top --> cap seq1, seq2 letter.

			Position must still be tracked.
		*/
		//cout << "Test, pos: " << pos << endl;

		if (tbmp[pos][1] == 1) {
			result_sequence1 += get_seq1_letter(pos);
			if (pos > x_koko) { // didn't go below seq2
				result_sequence2 += get_seq2_letter(pos);
			} else{
				result_sequence2 += "*";
			}
			pos = pos - x_koko - 1;
		} else if (tbmp[pos][0] == 1) { // horizontal move, add cap into seq2

			 //	check that won't go below the seq1:
			if (pos > 0) {
				result_sequence1 += get_seq1_letter(pos);
				result_sequence2 += "*";
				pos -= 1;
			} else{
				cout << "End1." << endl;
			}
		} else if (tbmp[pos][2] == 1){ // vertical move, add cap into seq1
			// Check that seq2 is read over or below the index.
			// Check that didn't go below the table bounds:
			if (pos >= x_koko) { // didn't go below the seq2 index
				result_sequence1 += "*";
				result_sequence2 += get_seq2_letter(pos);
			} else {
				cout << "End2. Pos = " << pos << ", x_koko = " << x_koko << endl;
			}
			pos -= x_koko;
		} else {
			cout << "WARNING, NO TRACEBACK POINTER FOUND!: " << "mp[pos]: " << mp[pos] << ", pos = " << pos << endl;
			print_matrix();
			pos--;
		}
		//mp[pos] = 9;
		//cout << "Printing the walked matrix: " << endl;
		//print_matrix();
	} // End While

	// cout << "Errors: " << error_flag << endl;

	result_sequence1 = reverse(result_sequence1);
	result_sequence2 = reverse(result_sequence2);

	//cout << "SEQ1: " << result_sequence1 << endl;
	//cout << "SEQ2: " << result_sequence2 << endl;

	Sequence seq1(result_sequence1);
	Sequence seq2(result_sequence2);

	vector<Sequence> *result = new vector<Sequence>;
	result->push_back(seq1);
	result->push_back(seq2);
	return result;

}

/* Override, because different initialization is wanted: */

void Global_aligner::initialize_matrix(string seq1, string seq2) {
	if (verbalize) {cout << "Alustan matriisin:" << endl;}
	y_koko = seq2.size() + 1;
	x_koko = seq1.size() + 1;

	if (verbalize) {cout << "\tx koko: " << x_koko << endl;}
	if (verbalize) {cout << "\ty koko: " << y_koko << endl;}

	mp = new int[x_koko * y_koko];
	if (verbalize) {cout << "\tLuon traceback matriisin." << endl;}
	tbmp = new int*[x_koko * y_koko];

    // The wanted initialization is as follows:
	//   ------
	//	| 0-1-2
	//	|-1....
	//	|-2....
	//  This means that the caps become more expensive.
	int add = 0;
	for (int x = 0; x < x_koko; x++) {
		mp[x] = 0 + add;
		add--;
	}
	add = -1;
	for (int y = 1; y < y_koko; y++) {
		mp[x_koko * y] = 0 + add;
		add--;
	}
	if (verbalize) {cout << "\tAlustan traceback matriisin..." << endl;}
	for (int i = 0; i < x_koko * y_koko; i++) {
		tbmp[i] = new int[2]; // [left,ldiag,top], possible values: 1, 0.
		for (int j = 0; j < 3; j++) {
			tbmp[i][j] = 0;	// values are initialized as zero
		}
	}
	if (verbalize) {cout << "Matrix initialized." << endl;}

}
