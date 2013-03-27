#include "html_output.h"

html_output::html_output(){
	_font = "<font face=\"courier new\" size=+1 color=blue>";
	font_ = "</font>";
	_match = "<b><font color=\"lime\">";
	match_ = "</b></font>";
	_cap = "<font color=\"red\">";
	cap_ = "</b></font>";
	matches = 0;
	mismatches = 0;
	caps = 0;

	outhtml.open ("output.html", ofstream::out | ofstream::trunc);
}

html_output::~html_output(){
	outhtml.close();
}

void html_output::print_headers(){

	outhtml << "<html>\n<head></head>\n<body>\n" << endl;
	outhtml << _font << "\n<br/><br/><pre>" << endl;
}
void html_output::output_html(vector<Sequence> *s){
	/*
        Tulostaa html-dokumentin siten ett‰ sekvenssit alkiot
        joita rinnastetuissa sekvensseiss‰ on eniten maalataan
        vihre‰ll‰ ja aukot punaisella. Ep‰t‰sm‰ykset tulevat sinisell‰.

		K‰yd‰‰n kirjain kerrallaan l‰pi vektorissa olevia sekvenssej‰.
	*/

	/*
		Katsotaan mit‰ kirjainta on miss‰kin riviss‰ eniten.
		Tehd‰‰n taulukko jonka koko on sekvenssin koko.
		Taulukkoon ker‰t‰‰n eri kirjaimia.
		Taulukon arvot edustavat kirjainten m‰‰r‰‰.
		Taulukon indeksi edustaa kirjainta.
	*/

	int seq_len  = s->at(0).length();	// Sekvenssin pituus
	int vec_len  = s->size();			// Sekvenssien m‰‰r‰

	char *letters = new char[seq_len]; // tauluun tulee kirjaimet joita riveiss‰ eniten.
	char *let_tmp = new char[vec_len*2]; // tauluun tulee kirjaimet joita riviss‰ eniten. kirjain, m‰‰r‰
	bool found;
	int index = 0;

	for(int z=0;z<=2*vec_len;z++){
		let_tmp[z] = 0;
	}

	for(int i=0;i<seq_len;i++){
		for(int j=0;j<vec_len;j++){
			found = false;
			for(int z=0;z<vec_len;z++){
				if(let_tmp[z]==s->at(j).toString()[i]){ // Kirjain on let_tmp taulussa.
					found = true;
					let_tmp[z+vec_len]++; 					// -->Lis‰t‰‰n sen laskuria.
					break;
				}
			}
			if(found==false){ 							// Kirjain ei ole let_tmp taulukossa.
				if(s->at(j).toString()[i]!='*'){			 // Mik‰li kirjain ei ole aukko-merkki..
					let_tmp[index] = s->at(j).toString()[i]; // ..lis‰t‰‰n se let_tmp-taulukkoon.
					let_tmp[index+vec_len]++;
					index++;
				}
			}
		}
		/*
			Katsotaan mit‰ kirjainta oli eniten ja lis‰t‰‰n se letters-taulukkoon rivin kohdalle.
		*/
		int max = 0;	// indeksoi kirjaimen jota let_tmp-taulukossa eniten.
		for(int z=0;z<vec_len;z++){	// K‰yd‰‰n l‰pi let_tmp-taulukon arvo-osuus.
			if(let_tmp[z+vec_len] > let_tmp[max+vec_len]){
				max = z;		// Suurin arvo muistetaan. T‰ss‰ z vastaa arvoon liittyv‰‰ kirjainta.
			}
		}
		if(let_tmp[max+vec_len]>1){
			letters[i] = let_tmp[max];	// Talletetaan kirjain jota oli eniten riviss‰.
		} else{
			letters[i] = '*';
		}

		index=0;
		for(int z=0;z<=2*vec_len;z++){	// Nollataan let_tmp-taulukko.
			let_tmp[z] = 0;
		}
	}
    // Tulostetaan html-dokumentti:
	for(unsigned int i=0;i<s->size();i++){
		outhtml << "\n";
		for(int j=0;j<seq_len;j++){
			if(letters[j]==s->at(i).toString()[j] && letters[j]!='*'){
				outhtml << _match << s->at(i).toString()[j] << match_;
			} else if(s->at(i).toString()[j]=='*'){
				outhtml << _cap << s->at(i).toString()[j] << cap_;
			} else{
				outhtml << s->at(i).toString()[j];
			}
		}
	}

}
