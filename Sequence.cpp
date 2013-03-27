#include "Sequence.h"

Sequence::Sequence(string s){
	seq = s;
}

int Sequence::length() const{
	return seq.size();
}

int Sequence::countAT() const{
	int count = 0;
	for(unsigned int i=0;i<seq.size();i++){
		if(seq[i]=='A' || seq[i]=='a' || seq[i]=='T' || seq[i]=='t'){
			count++;
		}
	}
	return count;
}

int Sequence::countCG() const{
	int count = 0;
	for(unsigned int i=0;i<seq.size();i++){
		if(seq[i]=='A' || seq[i]=='a' || seq[i]=='T' || seq[i]=='t'){
			count++;
		}
	}
	return count;
}

Sequence Sequence::reverse() const{
	char temp;
	int j = seq.size()-1;
	string result = seq;

	for(int i = 0; i <= j; i++){
		temp = seq[i];
		result[i] = seq[j];
		result[j] = temp;
		j--;
	}
	return result;
}

bool Sequence::equals(Sequence b) const {
	if(b.toString()==seq){
		return true;
	}
	return false;
}

bool Sequence::operator==(Sequence b) const {
	return equals(b);
}

int Sequence::number_of_matches(Sequence b) const {
	/*
		Ottaa huomioon muutkin kuin vain emäkset (a,t,c,g).
	*/
	// katsotaan kumpi on lyhyempi:
	string tmp_seq1;
	string tmp_seq2;
	if(seq.size()<=b.toString().size()){
		tmp_seq1 = seq;
		tmp_seq2 = b.toString();
	} else{
		tmp_seq1 = b.toString();
		tmp_seq2 = seq;
	}
	// nyt tmp_seq1 sisältää lyhyemmän sekvenssin.
	int count = 0;
	for(unsigned int i=0;i<tmp_seq1.size();i++){
		if(tmp_seq1[i]==tmp_seq2[i]){
			count++;
		}
	}
	return count;
}

bool Sequence::is_complementary(Sequence b) const{
	if(seq.size() != b.toString().size()){
		return false;
	}

	string temp = b.toString();
	for(unsigned int i=0;i<seq.size();i++){
		if(temp[i]=='a' || temp[i]=='A'){
			temp[i] = 't';
		} else if(temp[i]=='t' || temp[i]=='T'){
			temp[i] = 'a';
		} else if(temp[i]=='c' || temp[i]=='C'){
			temp[i] = 'g';
		} else if(temp[i]=='g' || temp[i]=='G'){
			temp[i] = 'c';
		} else{
			return false;
		}
	}
	if(this->reverse() == temp){
		return true;
	}
	return false;
}

int Sequence::number_of_complementary_matches(Sequence b) const{
	// katsotaan kumpi on lyhyempi:
	string tmp_seq1;
	string tmp_seq2;
	if(seq.size()<=b.toString().size()){
		tmp_seq1 = seq;
		tmp_seq2 = b.toString();
	} else{
		tmp_seq1 = b.toString();
		tmp_seq2 = seq;
	}
	// nyt tmp_seq1 sisältää lyhyemmän sekvenssin.
	int count = 0;
	for(unsigned int i=0;i<tmp_seq1.size();i++){
		if((tmp_seq1[i]=='a' || tmp_seq1[i]=='A') && (tmp_seq2[i]=='t' || tmp_seq2[i]=='T')){
			count++;
		} else if((tmp_seq1[i]=='t' || tmp_seq1[i]=='T') && (tmp_seq2[i]=='a' || tmp_seq2[i]=='A')){
			count++;
		} else if((tmp_seq1[i]=='c' || tmp_seq1[i]=='C') && (tmp_seq2[i]=='g' || tmp_seq2[i]=='G')){
			count++;
		} else if((tmp_seq1[i]=='g' || tmp_seq1[i]=='G') && (tmp_seq2[i]=='c' || tmp_seq2[i]=='C')){
			count++;
		}
	}
	return count;
}

string Sequence::toString() const{
	return seq;
}
