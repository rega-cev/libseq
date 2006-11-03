#include "Align.h"

namespace {

/*
 * A straight-forward implementation of Neeldeman-Wunsh algorithm
 * for a pairwise global alignment, with the difference that a
 * gapOpenScore is not added at the beginning or end of the sequence
 * (like ClustalW does).
 */
template <typename Symbol>
double NeedlemanWunsh(std::vector<Symbol>& seq1,
		      std::vector<Symbol>& seq2,
		      double gapOpenScore,
		      double gapExtensionScore,
		      double **weightMatrix)
{
  /*
   * Remove gaps, and warn that we did.
   */
  bool foundGaps = false;
  for (unsigned i = 0; i < seq1.size(); ++i) {
    if (seq1[i] == Symbol::GAP) {
      if (!foundGaps) {
	std::cerr << "Warning: NeedlemanWunsh: sequence contained gaps? "
	             "Removed them." << std::endl;
	foundGaps = true;
      }
      seq1.erase(seq1.begin() + i);
      --i;
    }
  }

  for (unsigned i = 0; i < seq2.size(); ++i) {
    if (seq2[i] == Symbol::GAP) {
      if (!foundGaps) {
	std::cerr << "Warning: NeedlemanWunsh: sequence contained gaps? "
	             "Removed them." << std::endl;
	foundGaps = true;
      }
      seq2.erase(seq2.begin() + i);
      --i;
    }
  }

  const int seq1Size = seq1.size();
  const int seq2Size = seq2.size();

  double **dnTable = new double* [seq1Size+1];
  for (unsigned i = 0; i < seq1Size+1; ++i)
    dnTable[i] = new double[seq2Size+1];
  int    **gapsLengthTable = new int *[seq1Size+1];
  for (unsigned i = 0; i < seq1Size+1; ++i)
    gapsLengthTable[i] = new int[seq2Size+1]; // >0: horiz, <0: vert

  /*
   * compute table
   */
  dnTable[0][0] = 0;
  gapsLengthTable[0][0] = 0;
  for (unsigned i = 1; i < seq1Size+1; ++i) {
    dnTable[i][0] = dnTable[i-1][0] + gapExtensionScore;
    gapsLengthTable[i][0] = gapsLengthTable[i-1][0] + 1;
  }
  for (unsigned j = 1; j < seq2Size+1; ++j) {
    dnTable[0][j] = dnTable[0][j-1] + gapExtensionScore;
    gapsLengthTable[0][j] = gapsLengthTable[0][j-1] - 1;
  }

  for (unsigned i = 1; i < seq1Size+1; ++i) {
    for (unsigned j = 1; j < seq2Size+1; ++j) {

      double sextend
	= dnTable[i-1][j-1]
	+ weightMatrix[seq1[i-1].intRep()][seq2[j-1].intRep()];

      double horizGapScore = ((gapsLengthTable[i-1][j] > 0) || (j == seq2Size)
			      ? gapExtensionScore
			      : gapOpenScore + gapExtensionScore);
      double sgaphoriz
	= dnTable[i-1][j] + horizGapScore;

      double vertGapScore = (gapsLengthTable[i][j-1] < 0 || (i == seq1Size)
			     ? gapExtensionScore
			     : gapOpenScore + gapExtensionScore);
      double sgapvert
	= dnTable[i][j-1] + vertGapScore;

      if ((sextend >= sgaphoriz) && (sextend >= sgapvert)) {
	dnTable[i][j] = sextend;
	gapsLengthTable[i][j] = 0;
      } else {
	if (sgaphoriz > sgapvert) {
	  dnTable[i][j] = sgaphoriz;
	  gapsLengthTable[i][j] = std::max(0, gapsLengthTable[i-1][j]) + 1;
	} else {
	  dnTable[i][j] = sgapvert;
	  gapsLengthTable[i][j] = std::min(0, gapsLengthTable[i][j-1]) - 1;
	}
      }
    }
  }

  /*
   * reconstruct best solution alignment.
   */
  int i = seq1Size+1, j = seq2Size+1;
  do {
    if (gapsLengthTable[i-1][j-1] == 0) {
      --i; --j;
    } else if (gapsLengthTable[i-1][j-1] > 0) {
      --i;
      seq2.insert(seq2.begin() + j-1, Symbol::GAP);
    } else {
      --j;
      seq1.insert(seq1.begin() + i-1, Symbol::GAP);
    }
  } while (i > 1 || j > 1);

  double score = dnTable[seq1Size][seq2Size];

  for (unsigned i = 0; i < seq1Size+1; ++i) {
    delete[] dnTable[i];
    delete[] gapsLengthTable[i];
  }

  return score;
}

};

namespace seq {
  

double Align(NTSequence& seq1, NTSequence& seq2,
	     double gapOpenScore, double gapExtensionScore,
	     double **weightMatrix)
{
  return NeedlemanWunsh(seq1, seq2, gapOpenScore, gapExtensionScore,
			weightMatrix);
}

double Align(AASequence& seq1, AASequence& seq2,
	     double gapOpenScore, double gapExtensionScore,
	     double **weightMatrix)
{
  return NeedlemanWunsh(seq1, seq2, gapOpenScore, gapExtensionScore,
			weightMatrix);
}

double ComputeAlignScore(const NTSequence& seq1, const NTSequence& seq2,
			 double gapOpenScore, double gapExtensionScore,
			 double **ntWeightMatrix)
{
  double score = 0;
  bool seq1Gap = true;
  bool seq2Gap = true;

  for (unsigned i = 0; i < seq1.size(); ++i) {
    if (seq1[i] == Nucleotide::GAP) {
      if (!seq1Gap) {
	seq1Gap = true;
	score += gapOpenScore;
      }
      score += gapExtensionScore;
    } else {
      seq1Gap = false;

      if (seq2[i] == Nucleotide::GAP) {
	if (!seq2Gap) {
	  seq2Gap = true;
	  score += gapOpenScore;
	}
	score += gapExtensionScore;
      } else {
	seq2Gap = false;

	score += ntWeightMatrix[seq1[i].intRep()][seq2[i].intRep()];
      }
    }
  }

  return score;
}

double AlignLikeAA(NTSequence& seq1, NTSequence& seq2,
		   int ORF, double gapOpenScore, double gapExtensionScore,
		   double **ntWeightMatrix, const AASequence& seqAA1,
		   const AASequence& seqAA2)
{
  NTSequence seq2ORFLead(seq2.begin(), seq2.begin() + ORF);
  seq2.erase(seq2.begin(), seq2.begin() + ORF);
  int aaLength = (seq2.size() - ORF) / 3;
  NTSequence seq2ORFEnd(seq2.begin() + aaLength*3, seq2.end());
  seq2.erase(seq2.begin() + aaLength*3, seq2.end());

  int firstNonGap = -1;
  int lastNonGap = -1;

  for (unsigned i = 0; i < seqAA1.size(); ++i) {
    if (seqAA1[i] == AminoAcid::GAP)
      seq1.insert(seq1.begin() + (i*3), 3, Nucleotide::GAP);
    if (seqAA2[i] == AminoAcid::GAP)
      seq2.insert(seq2.begin() + (i*3), 3, Nucleotide::GAP);
    else {
      if (firstNonGap == -1)
	firstNonGap = i*3;
      lastNonGap = i*3 + 3;
    }
  }

  for (unsigned i = 0; i < seq2ORFLead.size(); ++i)
    seq2[firstNonGap - seq2ORFLead.size() + i] = seq2ORFLead[i];
  for (unsigned i = 0; i < seq2ORFEnd.size(); ++i)
    seq2[lastNonGap + i] = seq2ORFEnd[i];

  return ComputeAlignScore(seq1, seq2, gapOpenScore, gapExtensionScore,
			   ntWeightMatrix);
}

bool haveGaps(const NTSequence& seq, int from, int to)
{
  for (unsigned i = std::max(from, 0); i < std::min((int)seq.size(), to); ++i)
    if (seq[i] == Nucleotide::GAP)
      return true;

  return false;
}

std::pair<double, int>
CodonAlign(NTSequence& ref, NTSequence& target, int maxFrameShifts,
	   double gapOpenScore, double gapExtensionScore,
	   double **ntWeightMatrix, double **aaWeightMatrix)
{
  /*
   * 1. translate the reference sequence
   * 2. for every open reading frame:
   *   - translate the target sequence
   *   - perform the alignment
   * 3. take the alignment with best score and align nucleotide
   *    sequence like amino acid sequence
   * 4. compute nucleotide alignment score
   * 5. make nucleotide sequence alignment, compare score, if difference
   *    too big then correct the frame shift and repeat.
   */
  AASequence refAA = AASequence::translate(ref);

  int bestFrameShift = -1;
  double bestScore = -1E10;
  AASequence bestRefAA;
  AASequence bestTargetAA;

  for (unsigned i = 0; i < 3; ++i) {
    int last = i + ((target.size() - i) / 3) * 3;
    AASequence targetAA
      = AASequence::translate(target.begin() + i, target.begin() + last);

    AASequence refCopyAA = refAA;
    double score = Align(refCopyAA, targetAA, gapOpenScore, gapExtensionScore,
			 aaWeightMatrix);

    if (score > bestScore) {
      bestFrameShift = i;
      bestScore = score;
      bestRefAA = refCopyAA;
      bestTargetAA = targetAA;
    }
  }

  NTSequence refCodonAligned = ref;
  NTSequence targetCodonAligned = target;

  double ntCodonScore = AlignLikeAA(refCodonAligned,
				    targetCodonAligned,
				    bestFrameShift,
				    gapOpenScore, gapExtensionScore,
				    ntWeightMatrix, bestRefAA,
				    bestTargetAA);

  NTSequence refNTAligned = ref;
  NTSequence targetNTAligned = target;
  double ntScore = Align(refNTAligned, targetNTAligned, gapOpenScore,
			 gapExtensionScore, ntWeightMatrix);

  if (ntScore - ntCodonScore > 100) {
    /*
     * a possible frameshift
     */
    if (maxFrameShifts) {
      /*
       * try to fix: walk through the nucleotide alignment, and find
       * an "isolated" gap that is not of size multiple of 3.
       */
      const int BOUNDARY=10;
      int seq2pos = 0;
      int refGapStart = 0;
      int targetGapStart = 0;
      bool fixed = false;

      for (unsigned i = 0; i < refNTAligned.size(); ++i) {
	if (refNTAligned[i] == Nucleotide::GAP) {
	  if (refGapStart == -1)
	    refGapStart = i;
	} else { 
	  if (refGapStart > 0) {
	    int refGapStop = i;

	    if ((refGapStop - refGapStart) % 3) {
	      /*
	       * check it is isolated: no gaps in either sequence around
	       * this gap
	       */
	      if (haveGaps(refNTAligned,
			   refGapStart - BOUNDARY, refGapStart)
		  || haveGaps(refNTAligned,
			      refGapStop, refGapStop + BOUNDARY)
		  || haveGaps(targetNTAligned,
			      refGapStart - BOUNDARY, refGapStart)
		  || haveGaps(targetNTAligned,
			      refGapStop, refGapStop + BOUNDARY)) {
		/*
		 * not isolated: skip this gap.
		 */
	      } else {
		/*
		 * fix it !
		 */
		target.insert(target.begin() + seq2pos,
			      3 - (refGapStop - refGapStart) % 3,
			      Nucleotide::N);
		fixed = true;
		break;		
	      }
	    }
	  }

	  refGapStart = -1;
	}

	if (targetNTAligned[i] == Nucleotide::GAP) {
	  if (targetGapStart == -1)
	    targetGapStart = i;
	} else {
	  if (targetGapStart > 0) {
	    int targetGapStop = i;

	    if ((targetGapStop - targetGapStart) % 3) {
	      /*
	       * check it is isolated: no gaps in either sequence around
	       * this gap
	       */
	      if (haveGaps(refNTAligned,
			   targetGapStart - BOUNDARY, targetGapStart)
		  || haveGaps(refNTAligned, targetGapStop,
			      targetGapStop + BOUNDARY)
		  || haveGaps(targetNTAligned,
			      targetGapStart - BOUNDARY, targetGapStart)
		  || haveGaps(targetNTAligned,
			      targetGapStop, targetGapStop + BOUNDARY)) {
		/*
		 * not isolated: skip this gap.
		 */
	      } else {
		/*
		 * fix it !
		 */
		target.insert(target.begin() + seq2pos,
			      (targetGapStop - targetGapStart) % 3,
			      Nucleotide::N);
		fixed = true;
		break;
	      }
	    }
	  }

	  targetGapStart = -1;
	  ++seq2pos;
	}
      }

      if (!fixed)
	throw FrameShiftError(ntScore, ntCodonScore,
			      refNTAligned, targetNTAligned);
      else {
	std::pair<double, int> result
	  = CodonAlign(ref, target, maxFrameShifts - 1, gapOpenScore,
		       gapExtensionScore, ntWeightMatrix, aaWeightMatrix);
	++result.second;
	return result;
      }
    } else {
      throw FrameShiftError(ntScore, ntCodonScore,
			    refNTAligned, targetNTAligned);
    }
  } else {
    ref = refCodonAligned;
    target = targetCodonAligned;

    return std::make_pair(ntCodonScore, 0);
  }
}

FrameShiftError::FrameShiftError(double ntScore, double codonScore,
				 const NTSequence& ntRef,
				 const NTSequence& ntTarget)
{ }

FrameShiftError::~FrameShiftError() throw()
{ }

namespace weights {

extern double **IUB()
{
  static double rowA[] = { 5,-4,-4,-4,1,1,1,-4,-4,-4,-1,-1,-1,-4,-2 };
  static double rowC[] = { -4,5,-4,-4,1,-4,-4,1,1,-4,-1,-1,-4,-1,-2 };
  static double rowG[] = { -4,-4,5,-4,-4,1,-4,1,-4,1,-1,-4,-1,-1,-2 };
  static double rowT[] = { -4,-4,-4,5,-4,-4,1,-4,1,1,-4,-1,-1,-1,-2 };
  static double rowM[] = { 1,1,-4,-4,-1,-2,-2,-2,-2,-4,-1,-1,-3,-3,-1 };
  static double rowR[] = { 1,-4,1,-4,-2,-1,-2,-2,-4,-2,-1,-3,-1,-3,-1 };
  static double rowW[] = { 1,-4,-4,1,-2,-2,-1,-4,-2,-2,-3,-1,-1,-3,-1 };
  static double rowS[] = { -4,1,1,-4,-2,-2,-4,-1,-2,-2,-1,-3,-3,-1,-1 };
  static double rowY[] = { -4,1,-4,1,-2,-4,-2,-2,-1,-2,-3,-1,-3,-1,-1 };
  static double rowK[] = { -4,-4,1,1,-4,-2,-2,-2,-2,-1,-3,-3,-1,-1,-1 };
  static double rowV[] = { -1,-1,-1,-4,-1,-1,-3,-1,-3,-3,-1,-2,-2,-2,-1 };
  static double rowH[] = { -1,-1,-4,-1,-1,-3,-1,-3,-1,-3,-2,-1,-2,-2,-1 };
  static double rowD[] = { -1,-4,-1,-1,-3,-1,-1,-3,-3,-1,-2,-2,-1,-2,-1 };
  static double rowB[] = { -4,-1,-1,-1,-3,-3,-3,-1,-1,-1,-2,-2,-2,-1,-1 };
  static double rowN[] = { -2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 };

  static double *iub[] = { rowA, rowC, rowG, rowT, rowM, rowR, rowW, rowS,
			   rowY, rowK, rowV, rowH, rowD, rowB, rowN };

  return iub;
}

extern double **BLOSUM30()
{
  static double rowA[] = { 4,-3,0,0,-2,0,-2,0,0,-1,1,0,-1,1,-1,1,1,1,-5,-4,-7,0,0,0,0,0 };
  static double rowC[] = { -3,17,-3,1,-3,-4,-5,-2,-3,0,-2,-1,-3,-2,-2,-2,-2,-2,-2,-6,-7,0,0,0,-2,-2 };
  static double rowD[] = { 0,-3,9,1,-5,-1,-2,-4,0,-1,-3,1,-1,-1,-1,0,-1,-2,-4,-1,-7,0,0,0,5,-1 };
  static double rowE[] = { 0,1,1,6,-4,-2,0,-3,2,-1,-1,-1,1,2,-1,0,-2,-3,-1,-2,-7,0,5,0,0,-1 };
  static double rowF[] = { -2,-3,-5,-4,10,-3,-3,0,-1,2,-2,-1,-4,-3,-1,-1,-2,1,1,3,-7,0,-4,0,-3,-1 };
  static double rowG[] = { 0,-4,-1,-2,-3,8,-3,-1,-1,-2,-2,0,-1,-2,-2,0,-2,-3,1,-3,-7,0,-2,0,0,-1 };
  static double rowH[] = { -2,-5,-2,0,-3,-3,14,-2,-2,-1,2,-1,1,0,-1,-1,-2,-3,-5,0,-7,0,0,0,-2,-1 };
  static double rowI[] = { 0,-2,-4,-3,0,-1,-2,6,-2,2,1,0,-3,-2,-3,-1,0,4,-3,-1,-7,0,-3,0,-2,0 };
  static double rowK[] = { 0,-3,0,2,-1,-1,-2,-2,4,-2,2,0,1,0,1,0,-1,-2,-2,-1,-7,0,1,0,0,0 };
  static double rowL[] = { -1,0,-1,-1,2,-2,-1,2,-2,4,2,-2,-3,-2,-2,-2,0,1,-2,3,-7,0,-1,0,-1,0 };
  static double rowM[] = { 1,-2,-3,-1,-2,-2,2,1,2,2,6,0,-4,-1,0,-2,0,0,-3,-1,-7,0,-1,0,-2,0 };
  static double rowN[] = { 0,-1,1,-1,-1,0,-1,0,0,-2,0,8,-3,-1,-2,0,1,-2,-7,-4,-7,0,-1,0,4,0 };
  static double rowP[] = { -1,-3,-1,1,-4,-1,1,-3,1,-3,-4,-3,11,0,-1,-1,0,-4,-3,-2,-7,0,0,0,-2,-1 };
  static double rowQ[] = { 1,-2,-1,2,-3,-2,0,-2,0,-2,-1,-1,0,8,3,-1,0,-3,-1,-1,-7,0,4,0,-1,0 };
  static double rowR[] = { -1,-2,-1,-1,-1,-2,-1,-3,1,-2,0,-2,-1,3,8,-1,-3,-1,0,0,-7,0,0,0,-2,-1 };
  static double rowS[] = { 1,-2,0,0,-1,0,-1,-1,0,-2,-2,0,-1,-1,-1,4,2,-1,-3,-2,-7,0,-1,0,0,0 };
  static double rowT[] = { 1,-2,-1,-2,-2,-2,-2,0,-1,0,0,1,0,0,-3,2,5,1,-5,-1,-7,0,-1,0,0,0 };
  static double rowV[] = { 1,-2,-2,-3,1,-3,-3,4,-2,1,0,-2,-4,-3,-1,-1,1,5,-3,1,-7,0,-3,0,-2,0 };
  static double rowW[] = { -5,-2,-4,-1,1,1,-5,-3,-2,-2,-3,-7,-3,-1,0,-3,-5,-3,20,5,-7,0,-1,0,-5,-2 };
  static double rowY[] = { -4,-6,-1,-2,3,-3,0,-1,-1,3,-1,-4,-2,-1,0,-2,-1,1,5,9,-7,0,-2,0,-3,-1 };
  static double rowSTP[] = { -7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,-7,1,0,-7,0,-7,-7 };
  static double rowGAP[] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  static double rowZ[] = { 0,0,0,5,-4,-2,0,-3,1,-1,-1,-1,0,4,0,-1,-1,-3,-1,-2,-7,0,4,0,0,0 };
  static double rowU[] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  static double rowB[] = { 0,-2,5,0,-3,0,-2,-2,0,-1,-2,4,-2,-1,-2,0,0,-2,-5,-3,-7,0,0,0,5,-1 };
  static double rowX[] = { 0,-2,-1,-1,-1,-1,-1,0,0,0,0,0,-1,0,-1,0,0,0,-2,-1,-7,0,0,0,-1,-1 };

  static double *mat[] = { rowA, rowC, rowD, rowE, rowF, rowG, rowH, rowI,
			   rowK, rowL, rowM, rowN, rowP, rowQ, rowR, rowS,
			   rowT, rowV, rowW, rowY, rowSTP, rowGAP,
			   rowZ, rowU, rowB, rowX };

  return mat;
}

};

};
