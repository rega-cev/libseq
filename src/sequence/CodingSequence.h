// This may look like C code, but it's really -*- C++ -*-
#ifndef CODING_SEQUENCE_H_
#define CODING_SEQUENCE_H_

#include <iostream>

#include "NTSequence.h"
#include "AASequence.h"

namespace seq {

/**
 * A coding sequence represents a nucleotide sequence that codes for
 * an amino acid sequence (an oligo- or polypeptide).
 *
 * It is useful when one wants to track the effect of changes in the
 * nucleotide sequence for the amino acid sequence, and to investigate
 * properties of nucleotide mutations.
 */
class CodingSequence
{
 public:
  /**
   * Construct a coding sequence based on the given nucleotide
   * sequence. The sequence must be translatable as per
   * AASequence::translate(const NTSequence&).
   */
  CodingSequence(const NTSequence& aNtSequence);

  /**
   * Get the nucleotide sequence.
   */
  const NTSequence& ntSequence() const { return ntSequence_; }

  /**
   * Get the amino acid sequence.
   *
   * If needed, the amino acid sequence is updated to reflect changes
   * in the nucleotide sequence.
   */
  const AASequence& aaSequence() const;

  /**
   * Change a nucleotide at a given position in the nucleotide sequence to
   * a new value.
   */
  void changeNucleotide(int pos, Nucleotide value);

  /**
   * Investigate whether a give nucleotide mutation is synonymous or
   * non-synonymous with respect to the amino acid sequence.
   */
  bool isSynonymousMutation(int pos, Nucleotide value) const;

 protected:
  void updateAASequence() const;

 private:
  NTSequence         ntSequence_;
  mutable AASequence aaSequence_;

  bool               isDirty() const { return dirty_ != D_CLEAN; }

  mutable int        dirty_;

  static const int   D_CLEAN = -1;
  static const int   D_COMPLETE = -2;
};

};

#endif // CODING_SEQUENCE_H_
