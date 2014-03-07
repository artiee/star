#ifndef GLOBAL_ALIGNER_H
#define GLOBAL_ALIGNER_H

#include <iostream>
#include <string>
#include <vector>
#include "Sequence_aligner.h"
#include "html_output.h"
#include "Sequence.h"

using namespace std;
/**

	Global_aligner is a global alignment algorithm.

	Aligns the two sequences globally as introduced in
	http://www.sbc.su.se/~per/molbioinfo2001/dynprog/adv_dynamic.html

	Scoring scheme:
		match = 2p
		mismatch = -1p
		cap	= -2p

  */
class Global_aligner : public Sequence_aligner {

public:
  Global_aligner(Sequence a, Sequence b);
  Global_aligner(string a, string b);
  vector<Sequence>* traceback();
  void initialize_matrix(string seq1, string seq2);

private:
  html_output * html_out;

};
#endif
