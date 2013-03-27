#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <iostream>
#include <string>

using namespace std;
/*
	Represents the DNA sequence.

*/

class Sequence {

public:
	/* Constuctor, string s is the DNA sequence that is saved. */
	Sequence(string s);

	/* Returns the lenght of the sequence. */
	int length() const;
	   //ae: -
	   //le: result == seq.size() as the sequence will be internally represented as a string.

	/* Returns the number of A and T bases */
	int countAT() const;
	   //ae: -
	   //le: result == sum of each b[n] == ‘A’ + sum of each b[n] == ‘T’ when n
	   //= 0..this.lenght()
	   //ak: O(n), depends only of the size of the sequence.

	/* Returns the number of C and G bases */
	int countCG() const;
	   //ae: -
	   //le: result == sum of each b[n] == ‘C’ + sum of each b[n] == ‘G’ when n
	   //= 0..this.lenght().
	   //ak: O(n), depends only of the size of the sequence.

	/* Reverses the sequence */
	Sequence reverse() const;
	   //ae: -
	   //le: elements of the sequence returned are in reverse order compared to
	   //the object.
	   //ak: read sequence: O(n) + copy it O(c). Thus O(n).

	/* Tests if two sequences are identical */
	bool equals(Sequence b) const;
	   //ae: b != null && this.lenght() > 0 && b.lenght() > 0.
	   //le: for each element in b holds: b[n] == this[n].
	   //ak: O(n)

	/* Same as equals(Sequence b) */
	bool operator==(Sequence b) const;

	/* Tests if the given sequence is complementary to this sequence */
	bool is_complementary(Sequence b) const;
	   //ae: b != null && this.lenght() > 0 && b.lenght() > 0.
	   //le: result == true if each for each element in b is done conversion: A
	   //T, T -> A, C -> G, G -> C and the result of the conversion equals
	   //to this.reverse().
	   //ak: conversion: O(n) + reverse: O(n). Thus O(2n).

	/* Returns how many bases are same in the sequence b */
	int number_of_matches(Sequence b) const;
	   //ae: b != null && this.lenght() == && b.lenght() > 0.
	   //le: result >= 0 && result <= b.lenght().
	   //ak: O(n)

	/* Returns how many bases match complementary to sequence */
	int number_of_complementary_matches(Sequence b) const;
	   //ae: b != null && this.lenght() == b.lenght() > 0.
	   //le: result >= 0 && result <= b.lenght().
	   //ak: O(n)
	string toString() const;

private:
	string seq;

};
#endif
