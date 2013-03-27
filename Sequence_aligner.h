#ifndef SEQUENCE_ALIGNER_H
#define SEQUENCE_ALIGNER_H

#include <iostream>
#include <string>
#include <vector>
#include "Sequence.h"

using namespace std;

/*

The abstract sequence aligner class.
Based on the tutorial at http://www.sbc.su.se/~per/molbioinfo2001/dynprog/dynamic.html

*/


class Sequence_aligner
{
public:
	Sequence_aligner();                  // Default constructor. If used the sequences must be set using set_sequence..-functions.
	Sequence_aligner(Sequence a, Sequence b);
	Sequence_aligner(string s1, string s2) ;
	~Sequence_aligner();

	void initialize_matrix(string seq1, string seq2); //Initialises the matrix. The first column and first row are filled with nulls by default. Different initialisation can be done by overwriting this method.
	void fill_matrix();                     // Fills the matrix by using the scoring set to class.
	void print_matrix();                    // Prints the values from matrice mp.
	virtual vector<Sequence>* traceback();  // Walks through the matrix and aligns the sequences in the way programmer wants.

	int get_cap_penalty() const;            // Tells how much cap costs in sequence.
	int get_match_score() const;
	int get_mismatch_score() const;
	void set_cap_penalty(const int cp);     //
	void set_match_score(const int ms);     //
	void set_mismatch_score(const int mms); //
	void set_sequence1(string a);           // The sequences that are aligned must be set if the default constructor is used.
	void set_sequence2(string b);
	void set_verbalize(bool v);             // If true program tells what it is doing.

	string reverse( const string& strReversible ); // Returns a copy that is reversed of string.

protected:
	int *mp;	// Pointer to valuematrix.

	/* Table sized of matrice mp. Each element is a pointer to a int-table.
	   Table is used in traback phase. Int-table holds values which tell
	   where to brach in traback.
	 		 ___________
	tbmp->	|p0|p1|p2|p3| p0--->[p0.1|p0.2|p0.3|...]
			|p4|p5|p6|p7|
			|p8|p9|pA|pB|
	 ¨¨¨¨¨¨¨¨¨¨¨
	*/
	int **tbmp; // Traceback matrices pointer.

	int y_koko;	   // Matrices y-size
	int x_koko;    // Matrices x-size
	string seq1;   // Sequence 1
	string seq2;
	int cap_penalty;	// How much costs to add a cap (empty space) to a sequence.
	int match_score;
	int mismatch_score;
	static const int MIN = -9999;
	bool verbalize; // Tells what is done while execution.

	/*
		Functions that involve walking in matrice:
	*/
	int left(int pos);	// Walks one left in matrice.
	int ldiag(int pos);	// ldiag = left diagonal (tässä ylävasen). Walks one top left.
	int top(int pos);	// Walks one up.

	/*
		Functions which check if element next to matrices position is the same.
	*/
	int left_match(int pos);  // Element in pos matches to a element one left from position.
	int ldiag_match(int pos); // Element in pos matches to a element one left and up (top left) from position.
	int top_match(int pos);   // Element in pos matches to a element one up from position.

	char get_seq1_letter(int pos) const; // Returns the letter from sequence1 (which can be thought to be on top column.)
	char get_seq2_letter(int pos) const; // Returns the letter from sequence2 (which can be thought to be on left row.)

    int max(int a, int b, int c) const;     // Returns max of these 3 integers.

	/*
	The following functions assume that i>0 and j>0.

    The purpose is to calculate from elements next to position
    (left, top left and top) which value they hold when the
    cap_penalty and match_score are taken into account.
    The left and top incorporate cap_penalty and top left gets
    match_score or mismatch score.
	*/

	int f1(int i, int j); // ldiag. Counts calue for left top.	(match/mismatch_score)
	int f2(int i, int j); // top. (cap_penalty)
	int f3(int i, int j); // left. (cap_penalty)

	int match(int i, int j) const; // Returns match_score if match, otherwise mismatch_score.
	int letter_match(int pos) const; // Returns 1 if letters match in sequence1 and sequence2 in given position in matrice mp. (The sequence1 is thought to be on top and sequence2 left to the matrice.)

};
#endif
