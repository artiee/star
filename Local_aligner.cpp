#include "Local_aligner.h"

Local_aligner::Local_aligner(Sequence a, Sequence b){
	set_sequence1(a.toString());
	set_sequence2(b.toString());
	set_match_score(2);
	set_mismatch_score(-1);
	set_cap_penalty(-2);
	initialize_matrix(a.toString(), b.toString());
}

Local_aligner::Local_aligner(string a, string b){
	set_sequence1(a);
	set_sequence2(b);
	set_match_score(2);
	set_mismatch_score(-1);
	set_cap_penalty(-2);
	initialize_matrix(a, b);
}

vector<Sequence>* Local_aligner::traceback(){
/*
ae: matrix in class Sequence_aligner is initialised and filled with
    start values.
le: result contains two sequences that are in alignment.
aikakompleksisuus: walks through the matrix. Time complexity is O(N*M).
*/
	string result_sequence1;
	string result_sequence2;

	/*
		aloitetaan suurimmasta alkiosta:
	*/
	int max = 0;
	int pos = 0;
	for(int i=0; i<x_koko*y_koko; i++){
		if(mp[i]>max){
			max = mp[i];
			pos = i;
		}
	}

	if(verbalize) cout << "Traceback started, please wait. Pos: " << pos << endl;

	while(pos>0){
		// jos match niin mene sinne

		/*
			Polku katsotaan tbmp listasta.
			Kun mennään ldiag --> lisää kummankin sekvenssin kirjain.
			left --> kappi seq2, seq1 kirjain.
			top --> kappi seq1, seq2 kirjain.

			Positionista tulee yhä pitää kirjaa.
		*/

		if(tbmp[pos][1] == 1){
			result_sequence1 += get_seq1_letter(pos);
			if(pos > x_koko){ // ei mennä ali seq2:sta
				result_sequence2 += get_seq2_letter(pos);
			} else{
				result_sequence2 += "*";
			}
			pos = pos - x_koko - 1;
		} else if(tbmp[pos][0] == 1){ // horisontaali liike, lisätään seq2:een cap.
			/*
				tarkistetaan ettei mennä seq1:en indeksin ali:
			*/
			if(pos>0){
				result_sequence1 += get_seq1_letter(pos);
				result_sequence2 += "*";
				pos -= 1;
			} else{
				if(verbalize) cout << "End horizontal reason." << endl;
			}
		} else if(tbmp[pos][2] == 1){ // vertikaali liike, lisätään seq1:een cap.
			/*
				Ongelma: seq2:sta luetaan yli tai ali indeksin.
				tarkistetaan ettei mennä taulukon rajojen ali:
			*/
			if(pos > x_koko){ // ei mennä ali seq2:sen indeksin.
				result_sequence1 += "*";
				result_sequence2 += get_seq2_letter(pos);
			} else{
				cout << "Debug: pos, x_koko:" << pos << ", " << x_koko << endl;
			}
			pos -= x_koko;
		} else{
			cout << "WARNING, NO TRACEBACK POINTER FOUND!" << endl;
			pos--;
		}

	} // End While

	// käännä result ympäri:
	result_sequence1 = reverse(result_sequence1);
	result_sequence2 = reverse(result_sequence2);

	if(verbalize) cout << "SEQ1: " << result_sequence1 << endl;
	if(verbalize) cout << "SEQ2: " << result_sequence2 << endl;

	Sequence seq1(result_sequence1);
	Sequence seq2(result_sequence2);

	vector<Sequence> *result = new vector<Sequence>;
	result->push_back(seq1);
	result->push_back(seq2);
	return result;
}
