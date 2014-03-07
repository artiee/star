#include "Local_aligner.h"

Local_aligner::Local_aligner(Sequence a, Sequence b) {
  set_sequence1(a.toString());
  set_sequence2(b.toString());
  set_match_score(2);
  set_mismatch_score(-1);
  set_cap_penalty(-2);
  initialize_matrix(a.toString(), b.toString());
}

Local_aligner::Local_aligner(string a, string b) {
  set_sequence1(a);
  set_sequence2(b);
  set_match_score(2);
  set_mismatch_score(-1);
  set_cap_penalty(-2);
  initialize_matrix(a, b);
}

vector<Sequence>* Local_aligner::traceback() {
/**
  * assume: matrix in class Sequence_aligner is initialised and filled with
  * start values.
  * result: result contains two sequences that are in alignment.
  * time-complexity: walks through the matrix. Time complexity is O(N*M).
  */
  string result_sequence1;
  string result_sequence2;

  // start from the biggest element:
  int max = 0;
  int pos = 0;
  for (int i = 0; i < x_size * y_size; i++) {
  	if (mp[i] > max) {
      max = mp[i];
	  pos = i;
	}
  }

  if (verbalize) cout << "Traceback started, please wait. Pos: " << pos << endl;

  while (pos > 0) {
	// if match, go to that direction

	/*
		The path is checked from tbmp -list.
		When ldiag -> add letter from both sequences.
		left --> cap seq2, seq1 letter.
		top --> cap seq1, seq2 letter.

		Position is still tracked.
	*/

    if (tbmp[pos][1] == 1) {
	  result_sequence1 += get_seq1_letter(pos);
	  if (pos > x_size) { // don't go below seq2
	    result_sequence2 += get_seq2_letter(pos);
	  } else {
	    result_sequence2 += "*";
	  }
	  pos = pos - x_size - 1;
	} else if (tbmp[pos][0] == 1) { // horizontal move, add cap into seq2.
	  // check that didn't go below seq1 index:
	  if (pos > 0) {
        result_sequence1 += get_seq1_letter(pos);
        result_sequence2 += "*";
        pos -= 1;
	  } else {
	    if (verbalize) cout << "End horizontal reason." << endl;
	  }
	} else if (tbmp[pos][2] == 1) { // vertical move, add cap into seq1.
	  // Problem: seq2 is read over or below the index.
	  // Must check that didn't go below the table bounds.
      if (pos > x_size) { // don't go below seq2's index
        result_sequence1 += "*";
        result_sequence2 += get_seq2_letter(pos);
	  } else {
	    cout << "Debug: pos, x_size:" << pos << ", " << x_size << endl;
	  }
	  pos -= x_size;
	} else {
	  cout << "WARNING, NO TRACEBACK POINTER FOUND!" << endl;
	  pos--;
	}

  } // End While

  result_sequence1 = reverse(result_sequence1);
  result_sequence2 = reverse(result_sequence2);

  if (verbalize) cout << "SEQ1: " << result_sequence1 << endl;
  if (verbalize) cout << "SEQ2: " << result_sequence2 << endl;

  Sequence seq1(result_sequence1);
  Sequence seq2(result_sequence2);

  vector<Sequence> *result = new vector<Sequence>;
  result->push_back(seq1);
  result->push_back(seq2);
  return result;
}
