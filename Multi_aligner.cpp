#include "Multi_aligner.h"

Multi_aligner::Multi_aligner(Sequence s[], int count) {
	for (int i = 0; i < count; i++) {
		seq.push_back(s[i]);
	}
	verbalize = false;
}
Multi_aligner::Multi_aligner(string s[], int count) {

	// Transform strings into sequences:
	for (int i = 0; i < count; i++) {
		Sequence temp(s[i]);
		seq.push_back(temp);
	}
	verbalize = false;
}

// for testing purposes:
void Multi_aligner::print() {
	for (unsigned int i = 0; i < seq.size(); i++) {
		cout << seq[i].toString() << endl;
	}
}

/**
  *	Finds the sequence that provides the best match to other sequences
  *	ae: Sequence[1..n] != null
  * le: for each sequence in Sequence[] s, local/global alignment with s[n] and result provides most matches.
  *
  *	Miten summat lasketaan t�ss�:
  *	http://staff.cs.utu.fi/kurssit/johdatus_bioinformatiikkaan_I/syksy_2004/rinnastukset2.pdf
  *	Sivulla k�ytetty pisteytys:
  *		t�sm�ys = 1
  *		ep�t�sm�ys = -1
  *		aukko = -2
  *		aukko molemmissa = 0
  *
  *
  *		Kaikkia sekvenssej� verrataan muiden sekvenssien kanssa.
  *		Sekvenssin tulisi olla mahdollisimman samanlainen kaikkien muiden kanssa.
  *		T�m� tarkoittaa ett� valitaan sekvenssi jonka count_score(..) yhteenlaskettu
  *		summa muiden sekvenssien kanssa on suurin.
  */
int Multi_aligner::find_center_sequence() {
	int score = 0;
	int best_score = -9999;
	int best_sequence = 0;
	vector<Sequence> *temp;

	for (unsigned int i = 0; i < seq.size(); i++) {
		score = 0;
		for (unsigned int j = 0; j < seq.size(); j++) {
				if (i != j) {
					Global_aligner ga(seq[i], seq[j]);
					ga.fill_matrix();
					temp = ga.traceback();
					score += count_score(temp->at(0), temp->at(1), 1,-1,-2); /*last three params: match, mismatch, cap*/
				}
		}
		if (score > best_score) {
			best_score = score;
			best_sequence = i;
		}
		if (verbalize) {cout << "\tPisteet yhteensa: " << score << endl;}
	}
	if (verbalize) {cout << "Paras pistem. " << best_score << endl;}
	if (verbalize) {cout << "Sekvenssilla " << best_sequence << endl;}

	return best_sequence;

}

/**
 *	Uses the STAR-method to make the multi alignment	
 */
vector<Sequence>* Multi_aligner::STAR() {
	vector<Sequence> *temp;
	vector<Sequence> *result = new vector<Sequence>;

	if (verbalize) {cout << "Etsitaan keskimmainen sekvenssi:" << endl;}
	unsigned int center = find_center_sequence();

    // Rinnastetaan viel� kerran parhaimmin sopivan sekvenssin kanssa kaikki muut:

	if (verbalize) {cout << "Rinnastetaan keskimmainen sekvenssi muiden kanssa:" << endl;}

	bool center_added = false;

	for (unsigned int i = 0; i < seq.size(); i++) {
		if (i != center) {
			Global_aligner ga(seq[center], seq[i]);
			ga.fill_matrix();
			temp = ga.traceback();
			if (center_added == false) {
				result->push_back(temp->at(0)); // keskimm�inen sekvenssi vektorin alkuun.
				center_added = true;
			}
			result->push_back(temp->at(1));
		}
	}
/*
	Lis�t��n sekvenssien loppuun tarvittava m��r� aukkoja pisimm�n mukaan.
		Ensin etsit��n pisin sekvenssi:
*/
	if (verbalize) cout << "Etsitaan pisin:" << endl;
	int len = 0;
	int len_tmp = 0;
	for (unsigned int i = 0; i < result->size(); i++) {
		len_tmp = result->at(i).toString().size();
		if (len_tmp > len) {
			len = len_tmp;
		}
	}
/*
	Kaikille t�t� lyhyemmille sekvensseille lis�t��n aukkoja ('*')
	kunnes ne ovat yht� pitki� kuin pisin sekvenssi.
*/
	if (verbalize) cout << "Lisataan aukot:" << endl;
	for (unsigned int i = 0; i < result->size(); i++) {
		len_tmp = result->at(i).toString().size();
		if (len_tmp < len) {
			string temp;
			temp = result->at(i).toString();
			while (temp.size() < len) {
				temp = temp + '*';
			}
			result->at(i) = temp;
		}
	}

	return result;
}

int Multi_aligner::count_score(Sequence s1, Sequence s2, int match, int mismatch, int cap) {
/*
	Miten summat lasketaan t�ss�:
	http://staff.cs.utu.fi/kurssit/johdatus_bioinformatiikkaan_I/syksy_2004/rinnastukset2.pdf
*/
	// katsotaan kumpi on lyhyempi (vaikka niiden pit�isi olla t�ss� jo samanpituisia):
	string tmp_seq1;
	string tmp_seq2;
	if (s1.toString().size()<=s2.toString().size()) {
		tmp_seq1 = s1.toString();
		tmp_seq2 = s2.toString();
	} else {
		tmp_seq1 = s2.toString();
		tmp_seq2 = s1.toString();
	}
	// nyt tmp_seq1 sis�lt�� lyhyemm�n sekvenssin.
	int score = 0;

	for (unsigned int i = 0; i < tmp_seq1.size(); i++) {
		if (tmp_seq1[i] == tmp_seq2[i] && tmp_seq1[i] != '*') {
			score += match;
		} else if (tmp_seq1[i] == tmp_seq2[i]) {
			score += 0; // jotta ei lasketa aukkoja kahteen kertaan.
		} else if (tmp_seq1[i] != tmp_seq2[i] && tmp_seq1[i] != '*' && tmp_seq2[i] != '*') {
			score += mismatch;
		} else {
			score += cap;
		}
	}
	if (verbalize) cout << "\tPisteet: " << score << endl;
	return score;
}
