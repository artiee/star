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
/**
   STAR-method for multialigning sequences.
   
   If one parameter is given, the sequences are tried to be read from a file
   of that name. In this case each line in the file is considered as a sequence.

   Two or more parameters starts alignment for the given parameters.
    
   If there are two sequences, the global alignment is used, otherwise the
   STAR-method.
*/
int main(int argc, char * argv[]) {

	cout << "*** Sequence alignment ***" << endl;

	if (argc > 2) {
		string *sekv = new string[argc - 1];

		for (int i = 1; i < argc; i++) {
			string tmp = argv[i];
			sekv[i - 1] = tmp;
		}

		vector<Sequence> *se;

		if (argc > 3) {
			cout << "\tUsing STAR-multialignment method.\n" << endl;
			Multi_aligner ma(sekv, argc - 1);
			cout << "Input sequences: " << endl;
			ma.print();
			se = ma.STAR();
		} else {
			cout << "\tUsing global alignment method.\n" << endl;
			Global_aligner sa(sekv[0], sekv[1]);
			sa.fill_matrix();
			se = sa.traceback();
		}

		// Print the results of alignment:

		cout << "\nResults:" << endl;

		for (unsigned int i = 0; i < se->size(); i++) {
			cout << "\t" << se->at(i).toString() << endl;
		}

		html_output html;
		html.print_headers();
		html.output_html(se);

		delete []sekv;

	} else if (argc == 2) {
		cout << "Trying to read the sequences from a file." << endl;
		ifstream seqfile (argv[1]);
		if (!seqfile.is_open()) {
			cout << "Failed to open the file." << endl;
			return 1;
		}

		vector<string> v;
		string line;

		while (getline(seqfile, line)) {
			v.push_back(line); // Add the line to the end
		}

		vector<Sequence> *se;
		string *sekv = new string[v.size()];

		for (int i = 0; i < v.size(); i++) {
    		cout << i + 1 << ": " << v[i] << endl;
			sekv[i] = v[i];
		}

		if (v.size() > 2) {
			cout << "\tUsing STAR-multialignment method." << endl;
			Multi_aligner ma(sekv, v.size());
			se = ma.STAR();
		} else {
			cout << "\tUsing global alignment method.\n" << endl;
			Global_aligner sa(sekv[0], sekv[1]);
			sa.fill_matrix();
			se = sa.traceback();
		}

		// Print the alignment results:
		cout << "\nResults:" << endl;

		for (unsigned int i = 0; i < se->size(); i++) {
			cout << "\t" << se->at(i).toString() << endl;
		}

		html_output html;
		html.print_headers();
		html.output_html(se);

		delete []sekv;
		

	} else {
		cout << "Give the sequences to be alignment as parameter, or name of the file that contains the sequences." << endl;
		return 1;
	}
    return EXIT_SUCCESS;
}



