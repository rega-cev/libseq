#include <ctype.h>

#include "Codon.h"
#include "CodingSequence.h"

namespace seq {

CodingSequence::CodingSequence(const NTSequence& aNtSequence)
  : ntSequence_(aNtSequence),
    aaSequence_(aNtSequence.size() / 3),
    dirty_(D_COMPLETE)
{ }

const AASequence& CodingSequence::aaSequence() const
{
  if (isDirty())
    updateAASequence();

  return aaSequence_;
}

void CodingSequence::changeNucleotide(int pos, Nucleotide value)
{
  // a small effort to avoid to avoid retranslation of the whole AA sequence.
  if (isDirty() && (dirty_ != D_COMPLETE))
    updateAASequence();

  ntSequence_[pos] = value;

  if (isDirty())
    dirty_ = D_COMPLETE;
  else
    dirty_ = pos;
}

bool CodingSequence::isSynonymousMutation(int pos, Nucleotide value) const
{
  if (isDirty())
    updateAASequence();

  const int aaPos = pos / 3;
  const int codonPos = pos % 3;

  NTSequence newcodon(ntSequence_.begin() + aaPos * 3,
		      ntSequence_.begin() + (aaPos * 3 + 3));
  newcodon[codonPos] = value;

  return (aaSequence_[aaPos] == Codon::translate(newcodon.begin()));
}

void CodingSequence::updateAASequence() const
{
  if (dirty_ == D_COMPLETE) {
    aaSequence_ = AASequence::translate(ntSequence_);
  } else {
    dirty_ /= 3;
    aaSequence_[dirty_]
      = Codon::translate(ntSequence_.begin() + (dirty_ * 3)); 
  }

  dirty_ = D_CLEAN;
}

};
