#include "Sequence_aligner.h"

int Sequence_aligner::left(int pos) {
    // Assuming that this can't go over the edge.
	if (pos - 1 < 0) {
	    return MIN;
	}
	return mp[pos - 1];
}
int Sequence_aligner::ldiag(int pos) {
	if (pos - x_koko - 1 < 0) {
	    return MIN;
	}
	return mp[pos - x_koko - 1];
}
int Sequence_aligner::top(int pos) {
	if (pos - x_koko < 0) {
	    return MIN;
	}
	return mp[pos - x_koko];
}
int Sequence_aligner::left_match(int pos) {
	if (get_seq1_letter(pos - 1) == get_seq2_letter(pos - 1)) {
	    return 1;
	}
	return 0;
}
int Sequence_aligner::ldiag_match(int pos) {
	if (get_seq1_letter(pos - x_koko-1) == get_seq2_letter(pos - x_koko - 1)) {
	    return 1;
	}
	return 0;
}
int Sequence_aligner::top_match(int pos) {
	if (get_seq1_letter(pos - x_koko) == get_seq2_letter(pos - x_koko)) {
	    return 1;
	}
	return 0;
}

char Sequence_aligner::get_seq1_letter(int pos) const {
	int z = pos;
	while (z > 0) {
		if (z - x_koko <= 0) {
		  break;
		}
		z = z - x_koko;
	}
	int x = z;
	//cout << "x: " << x << endl;
	return seq1[x - 1];
}

char Sequence_aligner::get_seq2_letter(int pos) const {
	/*
		1. Walk to the begin of row
		2. seq2:sta vastaava index = rivinro -1
			rivinro = index / x_koko


		  _ab
		 |012
		c|345

		x_koko    = 3
		c:n index = 3
		c:rivinro = 2
	*/
	int z = pos;
	while (z > x_koko) {
		z = z - x_koko;
	}
	if (z == x_koko) {
	    z = 0; // z not decremented in the loop, so we are in the start of the first row.
	}
	pos = pos - z; // we are in the begin of the row
	int rowNo = pos / x_koko;

	return seq2[rowNo - 1];
}


Sequence_aligner::Sequence_aligner() {
	//cout << "Creating sequence aligner without given sequences." << endl;
	verbalize = false;
}

Sequence_aligner::Sequence_aligner(Sequence a, Sequence b) {
	Sequence_aligner(a.toString(), b.toString());
}

Sequence_aligner::Sequence_aligner(string s1, string s2) {
	if (verbalize) {cout << "Creating sequence aligner." << endl;}
	seq1 = s1;
	seq2 = s2;
	// Scoring scheme:
	cap_penalty = 0;
	match_score = 1;
	mismatch_score = 0;
}

Sequence_aligner::~Sequence_aligner() {
	delete [] mp;
	// Destroy subtables:
	for (int i = 0; i < x_koko * y_koko; i++) {
		delete [] tbmp[i];
	}
	delete [] tbmp;
	if (verbalize) {
	    cout << "Sequence aligner destroyed." << endl;
	}
}

void Sequence_aligner::initialize_matrix(string seq1, string seq2) {
/*

Initialises the matrix. The first column and first row are filled with nulls by default.
Different initialisation can be done by overwriting this method.
*/
	if (verbalize) {
	    cout << "Initializing matrix:" << endl;
	}
	y_koko = seq2.size() + 1;
	x_koko = seq1.size() + 1;

	if (verbalize) {cout << "\tx koko: " << x_koko << endl;}
	if (verbalize) {cout << "\ty koko: " << y_koko << endl;}

	mp = new int[x_koko * y_koko];
	if (verbalize) {cout << "\tCreating traceback-matrix." << endl;}
	tbmp = new int*[x_koko * y_koko];


	/*   ------
		|000000
		|0.....
		|0.....
	*/
	for (int x = 0; x < x_koko; x++) {
		mp[x] = 0;
	}
	for (int y = 1; y < y_koko; y++) {
		mp[x_koko*y] = 0;
	}
	if (verbalize) {cout << "\tInitiatlizing traceback matrix..." << endl;}
	for (int i = 0; i < x_koko * y_koko; i++) {
		tbmp[i] = new int[2]; // [left,ldiag,top], possible values: 1, 0.
		for (int j = 0; j < 3; j++) {
			tbmp[i][j] = 0;	// initialize values with zero.
		}
	}
	if (verbalize) {cout << "Matrix initialized." << endl;}

}

void Sequence_aligner::fill_matrix() {
/*

Fills the matrix by using the scoring set to class.

*/

	if (verbalize) {cout << "Calculating values for matrix..." << endl;}
	if (mp == 0 || tbmp == 0) {
		cerr << "\tERROR: cannot access matrix." << endl;
		return;
	}
	int a, b, c;
	for (int x = 1; x < x_koko; x++) {
		for (int y = 1; y < y_koko; y++) {
			a = f1(x,y); // ldiag
			b = f2(x,y); // top
			c = f3(x,y); // left

            mp[x + x_koko * y] = max(a, b, c);
            /*
             * Check which f() returns the max value,
             * or if multiple f() return max-value.
             * Traceback pointers are saved for max values.
             */

			// Save the pointers:
			if (a > b) {
				if (a >= c) {
					tbmp[x + x_koko * y][1] = 1;
					if (a == c) {
						tbmp[x + x_koko * y][0] = 1;
					}
				} else {
				    tbmp[x + x_koko * y][0] = 1;
				}
			} else if (b >= c) {
				tbmp[x + x_koko * y][2] = 1;
				if (b == c) {
					tbmp[x + x_koko * y][0] = 1;
				}
			} else {
			    tbmp[x + x_koko * y][0] = 1;
			}

		}
	}
	/*
		First row and column are added the pointers too.
	*/
	for (int x = 1; x < x_koko; x++) {
		tbmp[x][0] = 1;
	}
	for (int y = 1; y < y_koko; y++) {
		tbmp[x_koko * y][2] = 1;
	}
	if (verbalize) {cout << "Values for matrix calculated. " << endl;}
}

/**
 *	Prints matrix mp. Used for testing.
 */
void Sequence_aligner::print_matrix() {
	cout << "    ";
	for (int i = 1; i < x_koko; i++) {
		cout << seq1[i - 1] << " ";
	}
	
	cout << endl;

	for (int y = 0; y < y_koko; y++) {
		if (y > 0) {
		    cout << seq2[y - 1] << " ";
		}
		else {
		    cout << "  ";
		}
		for (int x = 0; x < x_koko; x++) {
			cout << mp[x + x_koko * y] << " ";
		}
		cout << endl;
	}
}

/**
  *
  * Walks through the matrix and aligns the sequences in the way programmer wants.
  *
  * Example of cap (fin: aukko), match and mismatch:
  * Two sequences in alignment where "-" represents a cap,
  * match is when there is a same letter in both sequences
  * in the same position and mismatch when there is not:
  *
  * gctg-aa-cg
  * -ctataatc-
  */
vector<Sequence>* Sequence_aligner::traceback() {
	string result_sequence1;
	string result_sequence2;

	double matches = 0.0;
	int caps = 0;

    // Global alignment -algorithm starts from last position -> both sequences globally aligned.

	int pos = x_koko * y_koko;

	if (verbalize) {cout << "Traceback started, please wait. Pos: " << pos << endl;}


	while (pos > 0) {
		// If not match, go there,
		// fetch x- and y- components: (not needed here anymore?)
		int z = pos;
		int y = 0;
		while (z > 0) {
			if (z - x_koko <= 0) {
			    break;
			}
			y++;
			z = z - x_koko;
		}
		int x = z - 1;
		y = y - 1;

		if (verbalize) {
		    cout << "y = " << y << ", x = " << x << ". pos = " << pos << ". Arvo: " << mp[pos] << endl;
		}

		if (letter_match(pos)) {
			if (verbalize) {cout << "Match. Walking ldiag." << endl;}
			result_sequence1 += get_seq1_letter(pos);
			result_sequence2 += get_seq2_letter(pos);
			pos = pos - x_koko - 1;
			matches ++;
		}

		else if (left(pos) >= top(pos)) {
			if (left(pos) > ldiag(pos)) {
				if (verbalize) {cout << "Walking left. " << left(pos) << ">=" << ldiag(pos) << endl;}
				if (pos > 0) {
					result_sequence1 += get_seq1_letter(pos);
					result_sequence2 += "*";
					caps++;
				}
				pos -= 1;
			} else {
				if (verbalize) {cout << "Walking ldiag. " << endl;}

				if (pos > 0) {
					result_sequence1 += get_seq1_letter(pos);
					if (pos > x_koko) { // Don't go over from seq2
						result_sequence2 += get_seq2_letter(pos);
					} else {
						result_sequence2 += "*";
						caps++;
					}
				}
				pos = pos - x_koko - 1;
			}
		} else if (ldiag(pos) >= top(pos)) {
			if (verbalize) {cout << "Walking ldiag. " << endl;}

			if (pos > 0) {
				result_sequence1 += get_seq1_letter(pos);
				if (pos > x_koko) { // Don't go over from seq2
					result_sequence2 += get_seq2_letter(pos);
				} else {
					result_sequence2 += "*";
					caps++;
				}
			}
			pos = pos - x_koko - 1;
		} else {
			if (verbalize) {cout << "Walking up" << endl;}
			result_sequence1 += "*";
			result_sequence2 += get_seq2_letter(pos);
			pos -= x_koko;
			caps++;
		}

	} // End While

	if (verbalize) {
		cout << "Seq1 size: " << seq1.size() << " bytes." << endl;
		cout << "Seq2 size: " << seq2.size() << " bytes." << endl;
		cout << "Caps: " << caps << endl;
	}

	// Reverse result:
	result_sequence1 = reverse(result_sequence1);
	result_sequence2 = reverse(result_sequence2);

	Sequence seq1(result_sequence1);
	Sequence seq2(result_sequence2);

	vector<Sequence> *result = new vector<Sequence>;
	result->push_back(seq1);
	result->push_back(seq2);
	return result;
}

int  Sequence_aligner::get_cap_penalty() const {	
    return cap_penalty;
}

int  Sequence_aligner::get_match_score() const {
    return match_score;
}

int  Sequence_aligner::get_mismatch_score() const {
    return mismatch_score;
}

void Sequence_aligner::set_cap_penalty(const int cp) {
    cap_penalty = cp;
}

void Sequence_aligner::set_match_score(const int ms) {
    match_score = ms;
}

void Sequence_aligner::set_mismatch_score(const int mms) {
    mismatch_score = mms;
}

void Sequence_aligner::set_sequence1(const string a) {
	seq1 = a;
	if (verbalize) {cout << "Sequence 1 set." << endl;}
}

void Sequence_aligner::set_sequence2(string b) {
	seq2 = b;
	if (verbalize) {cout << "Sequence 2 set." << endl;}
}

void Sequence_aligner::set_verbalize(bool v) {
	verbalize = v;
}

/** returns max value of a, b, or c */
int Sequence_aligner::max(int a, int b, int c) const {
	if (a > b) {
		if (a > c) {
			return a;
		} else {
		    return c;
		}
	} else if (b > c) {
		return b;
	} else {
	    return c;
	}
}

/** Following functions assume that i > 0 and i > 0
  * ldiag
  */
int Sequence_aligner::f1(int i, int j) {
	int result = mp[(i - 1) + x_koko * (j - 1)] + match(i, j);
	return result;
}
/** top */
int Sequence_aligner::f2(int i, int j) {
	int result = mp[i + x_koko * (j - 1)] + cap_penalty;
	return result;
}
/** left */
int Sequence_aligner::f3(int i, int j) {
	int result = mp[(i - 1) + x_koko * j] + cap_penalty;
	return result;
}
/** Returns match_score if match, else mismatch_score. */
int Sequence_aligner::match(int i, int j) const {
	if(seq1[i - 1] == seq2[j - 1]) return match_score;
	return mismatch_score;
}

int Sequence_aligner::letter_match(int pos) const {
	if (get_seq1_letter(pos) == get_seq2_letter(pos) ) {
	    return 1;
	}
	return 0;
}
// Returns reversed copy from given string:
string Sequence_aligner::reverse(const string& strReversible) {
  string strCopy(strReversible);
   copy(strReversible.begin(), strReversible.end(), strCopy.rbegin());
   return strCopy;
}
