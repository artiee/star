#ifndef HTML_OUTPUT_H
#define HTML_OUTPUT_H

#include <iostream>
#include <fstream>
#include <vector>
#include "Sequence.h"

using namespace std;

/*
	T‰m‰ on tarkoitettu _vain_ sekvenssien rinnustuksien tulostamiseen.
	Rinnustusalgoritmien pit‰isi osata luoda t‰st‰ olio ja kutsua output_html metodia
	tulos sekvensseill‰‰n.
*/

class html_output
{
public:
	html_output();
	~html_output(); // sulkee outhtml tiedoston.
	void print_headers();
	void output_html(vector<Sequence> *s);

protected:

private:
	  string _font;
	  string font_;
	  string _match;
	  string match_;
	  string _cap;
	  string cap_;

	  int matches;
	  int mismatches;
	  int caps;

	  ofstream outhtml;
};
#endif
