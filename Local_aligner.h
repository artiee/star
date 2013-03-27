#ifndef LOCAL_ALIGNER_H
#define LOCAL_ALIGNER_H

#include <iostream>
#include <vector>
#include "Sequence_aligner.h"
#include "html_output.h"
#include "Sequence.h"

using namespace std;
/*
	Lokaali rinnastus voi aloittaa rinnastamisen mistä tahansa kohtaa sekvenssiä.
	Sopiva algoritmi lyhyen sekvenssin löytämiseen.

    Aligns the two sequences locally as introduced in
    http://staff.cs.utu.fi/kurssit/johdatus_bioinformatiikkaan_I/syksy_2004/rinnastukset.pdf
    page 7, slide 42.

*/
class Local_aligner : public Sequence_aligner {

public:
	Local_aligner(Sequence a, Sequence b);
	Local_aligner(string a, string b);
	vector<Sequence>* traceback();

private:
	html_output * html_out;

};
#endif
