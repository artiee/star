#include "Global_aligner.h"

Global_aligner::Global_aligner(Sequence a, Sequence b){
	set_sequence1(a.toString());
	set_sequence2(b.toString());
	set_match_score(2);
	set_mismatch_score(-1);
	set_cap_penalty(-3); // -2
	initialize_matrix(a.toString(), b.toString());
}

Global_aligner::Global_aligner(string a, string b){
	set_sequence1(a);
	set_sequence2(b);
	set_match_score(2);
	set_mismatch_score(-1);
	set_cap_penalty(-3);
	initialize_matrix(a, b);
}

vector<Sequence>* Global_aligner::traceback(){
/*
Aligns the sequence globally.

ae: matrix in class Sequence_aligner is initialised and filled with
    start values.
le: result contains two sequences that are in alignment.
aikakompleksisuus: walks through the matrix. Time complexity is O(N*M).

*/
	string result_sequence1;
	string result_sequence2;

	// Algoritmissa aloitetaan viimeisest‰ alkiosta --> kumpikin sekvenssi globaalisti rinnastettu.
	int pos = x_koko * y_koko - 1;

	if(verbalize) cout << "traceback started, please wait. Pos: " << pos << endl;

	while(pos>0){
		// jos match niin mene sinne
		// hae x ja y-komponentti:

		/*
			Polku katsotaan tbmp listasta.
			Kun menn‰‰n ldiag --> lis‰‰ kummankin sekvenssin kirjain.
			left --> kappi seq2, seq1 kirjain.
			top --> kappi seq1, seq2 kirjain.

			Positionista tulee yh‰ pit‰‰ kirjaa.
		*/
		//cout << "Test, pos: " << pos << endl;

		if(tbmp[pos][1] == 1){
			result_sequence1 += get_seq1_letter(pos);
			if(pos > x_koko){ // ei menn‰ ali seq2:sta
				result_sequence2 += get_seq2_letter(pos);
			} else{
				result_sequence2 += "*";
			}
			pos = pos - x_koko - 1;
		} else if(tbmp[pos][0] == 1){ // horisontaali liike, lis‰t‰‰n seq2:een cap.
			/*
				tarkistetaan ettei menn‰ seq1:en indeksin ali:
			*/
			if(pos>0){
				result_sequence1 += get_seq1_letter(pos);
				result_sequence2 += "*";
				pos -= 1;
			} else{
				cout << "End1." << endl;
			}
		} else if(tbmp[pos][2] == 1){ // vertikaali liike, lis‰t‰‰n seq1:een cap.
			/*
				Ongelma: seq2:sta luetaan yli tai ali indeksin.
				tarkistetaan ettei menn‰ taulukon rajojen ali:
			*/
			if(pos >= x_koko){ // ei menn‰ ali seq2:sen indeksin.
				result_sequence1 += "*";
				result_sequence2 += get_seq2_letter(pos);
			} else{
				cout << "End2. Pos = " << pos << ", x_koko = " << x_koko << endl;
			}
			pos -= x_koko;
		} else{
			cout << "WARNING, NO TRACEBACK POINTER FOUND!: " << "mp[pos]: " << mp[pos] << ", pos = " << pos << endl;
			print_matrix();
			pos--;
		}
		//mp[pos] = 9;
		//cout << "Printtaan kuljetun matriisin: " << endl;
		//print_matrix();
	} // End While

	// cout << "Errors: " << error_flag << endl;

	// k‰‰nn‰ result ymp‰ri:
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

/* Ylikuormitetaan koska halutaan erilainen alustus: */

void Global_aligner::initialize_matrix(string seq1, string seq2){
	if(verbalize) cout << "Alustan matriisin:" << endl;
	y_koko = seq2.size()+1;
	x_koko = seq1.size()+1;

	if(verbalize) cout << "\tx koko: " << x_koko << endl;
	if(verbalize) cout << "\ty koko: " << y_koko << endl;

	mp 		= new int[x_koko*y_koko];
	if(verbalize) cout << "\tLuon traceback matriisin." << endl;
	tbmp 	= new int*[x_koko*y_koko];

    // Haluttu alustus on seuraavanlainen:
	//   ------
	//	| 0-1-2
	//	|-1....
	//	|-2....
	//  T‰m‰ tarkoittaa ett‰ aukot tulevat kalliimmiksi.
	int add = 0;
	for(int x=0;x<x_koko;x++){
		mp[x] = 0 + add;
		add--;
	}
	add = -1;
	for(int y=1;y<y_koko;y++){
		mp[x_koko*y] = 0 + add;
		add--;
	}
	if(verbalize) cout << "\tAlustan traceback matriisin..." << endl;
	for(int i=0; i<x_koko*y_koko; i++){
		tbmp[i] = new int[2]; // [left,ldiag,top], mahdolliset arvot: 1, 0.
		for(int j=0; j<3; j++){
			tbmp[i][j] = 0;	// arvot alustetaan nolliksi.
		}
	}
	if(verbalize) cout << "Matriisin alustus valmis." << endl;

}
