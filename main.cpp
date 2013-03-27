#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include "Sequence_aligner.h"
#include "Global_aligner.h"
#include "Local_aligner.h"
#include "Multi_aligner.h"
#include "Sequence.h"

using namespace std;
/*

	STAR-menetelm‰ sekvenssien monirinnastusta varten.
	
	Jos ohjelmalle annetaan vain yksi parametri, 
    yritet‰‰n lukea sekvenssit parametrin nimisest‰ tiedostosta.
    T‰llˆin tiedoston jokainen rivi vastaa yht‰ sekvenssi‰.
    
    Kaksi tai useampi parametria k‰ynnist‰‰ rinnastuksen annetuille
    parametreille.
    
    Jos sekvenssej‰ on kaksi, k‰ytet‰‰n globaalia rinnastusta, muuten
    STAR-menetelm‰‰.

*/
int main(int argc, char * argv[]) {

	cout << "*** Sequence alignment ***" << endl;

	if(argc>2){
		string *sekv = new string[argc-1];

		for(int i=1; i<argc; i++){
			string tmp = argv[i];
			sekv[i-1] = tmp;
		}

		vector<Sequence> *se;

		if(argc > 3){
			cout << "\tUsing STAR-multialignment method.\n" << endl;
			Multi_aligner ma(sekv, argc-1);
			cout << "Input sequences: " << endl;
			ma.print();
			se = ma.STAR();
		} else{
			cout << "\tUsing global alignment method.\n" << endl;
			Global_aligner sa(sekv[0], sekv[1]);
			sa.fill_matrix();
			se = sa.traceback();
		}

		/*
			Tulostetaan rinnastuksen tulokset:
		*/

		cout << "\nResults:" << endl;

		for(unsigned int i=0;i<se->size();i++){
			cout << "\t" << se->at(i).toString() << endl;
		}

		html_output html;
		html.print_headers();
		html.output_html(se);

		delete []sekv;

	} else if(argc==2){
		cout << "Luen sekvenssit tiedostosta." << endl;
		ifstream seqfile (argv[1]);
		if(!seqfile.is_open()){
			cout << "Tiedoston avaaminen epaonnistui." << endl;
			return 1;
		}

		vector<string> v;
		string line;

		while(getline(seqfile, line)){
			v.push_back(line); // Add the line to the end
		}

		vector<Sequence> *se;
		string *sekv = new string[v.size()];

		for(int i = 0; i < v.size(); i++){
    		cout << i+1 << ": " << v[i] << endl;
			sekv[i] = v[i];
		}

		if(v.size() > 2){
			cout << "\tUsing STAR-multialignment method." << endl;
			Multi_aligner ma(sekv, v.size());
			se = ma.STAR();
		} else{
			cout << "\tUsing global alignment method.\n" << endl;
			Global_aligner sa(sekv[0], sekv[1]);
			sa.fill_matrix();
			se = sa.traceback();
		}

		/*
			Tulostetaan rinnastuksen tulokset:
		*/

		cout << "\nResults:" << endl;

		for(unsigned int i=0;i<se->size();i++){
			cout << "\t" << se->at(i).toString() << endl;
		}

		html_output html;
		html.print_headers();
		html.output_html(se);

		delete []sekv;
		

	} else {
		cout << "Anna parametrina rinnastettavat sekvenssit tai tiedoston nimi joka sisaltaa sekvenssit." << endl;
		return 1;
	}
    return EXIT_SUCCESS;
}



