#ifndef MULTI_ALIGNER_H
#define MULTI_ALIGNER_H

#include <iostream>
#include <vector>
#include "Local_aligner.h"
#include "Global_aligner.h"
#include "Sequence.h"
#include "html_output.h"

/*

   Aligns two or more sequences.
   If only two sequences are given the dynamic programming is used,
   otherwise the STAR-method as introduced in
   http://staff.cs.utu.fi/kurssit/johdatus_bioinformatiikkaan_I/syksy_2004/rinnastukset2.pdf
   page 5.
*/
class Multi_aligner {
public:
	Multi_aligner(Sequence s[], int count);
	Multi_aligner(string   s[], int count);

	void print(); // for testing

	/*
		Uses the STAR-method to make the multi alignment
	*/
	vector<Sequence>* STAR();
	/*
	   ae: -
	   le: Sequences given to the class constructor are aligned.
	   ak: To build the table and get the center of star: O(k^2 * n^2), where
		   k is the number of sequences. To do the multiple sequence
		   alignment: O(kn^2 + k^2 * l), where l is length after gaps).
		   (Source:
			http://chortle.ccsu.ctstateu.edu/CS580_70/LectureNotes/Lecture5/Lecture5.PDF) page 76.)
	*/
private:
	vector<Sequence> seq;
	bool verbalize;
	int count_score(Sequence s1, Sequence s2, int match, int mismatch, int cap);
	/*
		Finds the sequence that provides the best match to other sequences
	*/
	int find_center_sequence();
};
#endif
