#include <iostream>
#include <string>
//#include "sequence_aligner.h"
#include "html_output.h"

using namespace std;
/*

	adv_sequence_aligner on globaali rinnastusalgoritmi.
	Pisteytys:
		match = 2p
		mismatch = -1p
		cap	= -2p

*/
class adv_sequence_aligner : public sequence_aligner{
public:
	adv_sequence_aligner(string a, string b){
		set_sequence1(a);
		set_sequence2(b);
		set_match_score(2);
		set_mismatch_score(-1);
		set_cap_penalty(-2);
		initialize_matrix(a, b);
	}

/* ********************************************************************** */
	void traceback(){

		string result_sequence1;
		string result_sequence2;

		double matches = 0.0;
		int caps = 0;

		// adv. algoritmissa aloitetaan viimeisest‰ alkiosta --> kumpikin sekvenssi globaalisti rinnastettu.
		int pos = x_size * y_size - 1;
		// 2. aloita takaisinj‰ljitys.
		cout << "traceback started, please wait!! Pos: " << pos << endl;


		while(pos>0){
			// jos match niin mene sinne
			// hae x ja y-komponentti:
			int z = pos;
			int y=0;
			while(z>0){
				if(z-x_size<=0){break;}
				y++;
				z = z - x_size;
			}
			int x = z-1;
			y = y - 1;

			//cout << "y = " << y << ", x = " << x << ". pos = " << pos << ". Arvo: " << mp[pos] << endl;

			/*
				Polku katsotaan tbmp listasta.
				Kun menn‰‰n ldiag --> lis‰‰ kummankin sekvenssin kirjain.
				left --> kappi seq2, seq1 kirjain.
				top --> kappi seq1, seq2 kirjain.

				Positionista tulee yh‰ pit‰‰ kirjaa.
			*/

			//cout << "Test, pos: " << pos << endl;

			if(tbmp[pos][1] == 1){ //
				//cout << "got here1" << endl;
				result_sequence1 += get_seq1_letter(pos);
				if(pos > x_size){ // ei menn‰ ali seq2:sta
					result_sequence2 += get_seq2_letter(pos);
				} else{
					result_sequence2 += "*";
					caps++;
				}
				pos = pos - x_size - 1;
			} else if(tbmp[pos][0] == 1){ // horisontaali liike, lis‰t‰‰n seq2:een cap.
				//cout << "got here2" << endl;
				/*
					tarkistetaan ettei menn‰ seq1:en indeksin ali:
				*/
				if(pos>0){
					result_sequence1 += get_seq1_letter(pos);
					result_sequence2 += "*";
					pos -= 1;
					caps++;
				} else{
					cout << "End." << endl;
				}
			} else if(tbmp[pos][2] == 1){ // vertikaali liike, lis‰t‰‰n seq1:een cap.
				//cout << "got here3" << endl;
				/*
					Ongelma: seq2:sta luetaan yli tai ali indeksin.
					tarkistetaan ettei menn‰ taulukon rajojen ali:
				*/
				if(pos > x_size){ // ei menn‰ ali seq2:sen indeksin.
					result_sequence1 += "*";
					result_sequence2 += get_seq2_letter(pos);
					caps++;
				} else{
					cout << "End." << endl;
				}
				pos -= x_size;
			} else{
				cout << "WARNING, NO TRACEBACK POINTER FOUND!" << endl;
				pos--;
			}

			//cout << "SEQ1: " << result_sequence1 << endl;
			//cout << "SEQ2: " << result_sequence2 << endl;

			mp[pos] = 9;
			//print_matrix();
			cout << "Done." << endl;
		} // End While

		//cout << "Matches: " << matches << endl;
		//cout << "seq1 size: " << seq1.size() << " bytes." << endl;
		//cout << "seq2 size: " << seq2.size() << " bytes." << endl;
		//cout << "kappeja: " << caps << endl;

		//cout << "Match Percentage: " << (matches / (double)seq2.size())*100 << "%" << endl;

		// k‰‰nn‰ result ymp‰ri:
		cout << reverse(result_sequence1) << endl;
		cout << reverse(result_sequence2) << endl;

		html_out = new html_output();
		html_out->print_headers();
		html_out->output_html(result_sequence1, result_sequence2);

	}
private:
	html_output * html_out;

};
