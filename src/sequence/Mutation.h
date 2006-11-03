// This may look like C code, but it's really -*- C++ -*-
#ifndef MUTATION_H_
#define MUTATION_H_

#include <iostream>
#include <set>

#include "Nucleotide.h"
#include "AminoAcid.h"

namespace seq {

/**
 * A mutation in a sequence of type Char.
 *
 * \sa AAMutation, NTMutation
 */
template<typename Char>
class Mutation
{
public:
  /**
   * Specify a mutation at position pos, changing "from" to "to".
   */
  Mutation(int pos, Char from, Char to)
    : pos_(pos),
      from_(from),
      to_(to)
  { }

  /**
   * Specify a mutation at position pos, changing the character to "to".
   */
  Mutation(int pos, Char to)
    : pos_(pos),
      to_(to)
  { }

  /**
   * An invalid mutation, for which isValid() == false.
   */
  Mutation()
    : pos_(-1)
  { }

  /**
   * Mutation position.
   */
  int pos()   const { return pos_; }

  /**
   * Mutation 'from' character.
   */
  Char from() const { return from_; }

  /**
   * Mutation 'to' character.
   */
  Char to()   const { return to_; }

  /**
   * Return a mutation that reverses this mutation.
   */
  Mutation<Char> reverse() const {
    return Mutation<Char>(pos_, to_, from_);
  }

  /**
   * Is the mutation valid, i.e. not constructed using the default constructor.
   */
  bool isValid() const { return pos_ >= 0; }

  /**
   * Two mutations are equal if they change the position to the same character.
   */
  bool operator== (const Mutation<Char>& other) const {
    return (pos_ == other.pos_) && (to_ == other.to_);
  }

  /**
   * Facilitate sorting of mutations: by position, and at the same position
   * by internal representation of 'to'.
   */
  bool operator< (const Mutation<Char>& other) const {
    return (pos_ < other.pos_)
      || ((pos_ == other.pos_)
	  && (to_.intRep() < other.to_.intRep()));
  }

 private:
  int  pos_;
  Char from_, to_;
};

/**
 * A typedef for nucleotide mutations.
 */
typedef class Mutation<Nucleotide> NTMutation;
/**
 * A typedef for amino acid mutations.
 */
typedef class Mutation<AminoAcid> AAMutation;

extern std::set<AAMutation> readMutations(std::istream& mutationFile,
					  std::string prefix)
  throw (ParseException);

};

#endif // MUTATION_H_
