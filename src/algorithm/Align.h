// This may look like C code, but it's really -*- C++ -*-
#ifndef ALIGN_H_
#define ALIGN_H_

#include <NTSequence.h>
#include <AASequence.h>

/**
 * libseq namespace
 */
namespace seq {

/**
 * namespace for similarity weights matrices.
 */
namespace weights {
  /**
   * Similarity weights matrix for nucleotides.
   *
   * Compares also IUB ambiuguity codes, and is the matrix used by BLAST.
   *
   * Taken from: ftp://ftp.ncbi.nih.gov/blast/matrices/NUC.4.4
   */
  extern double **IUB();

  /**
   * Similarity weights matrix for amino acids.
   *
   * This is from the famous BLOSUM series of weight matrices, the one
   * that is the default use by ClustalX.
   *
   * From: ftp://ftp.ncbi.nih.gov/blast/matrices/BLOSUM30
   */
  extern double **BLOSUM30();
};

/**
 * Pair-wise align two nucleotide sequences, using a modified
 * NeedleMan-Wunsh algorithm.
 *
 * The two sequences seq1 and seq2 are aligned in-place: gaps are inserted
 * according to a global alignment, and they will have equal length.
 *
 * The algorithm is NeedleMan-Wunsh, with two popular modifications:
 *  - there is a different cost for opening a gap or for extending a gap.
 *  - there is no gap open cost for a gap at the beginning or the end.
 */
extern double Align(NTSequence& seq1, NTSequence& seq2,
		    double gapOpenScore = -10, double gapExtensionScore = -3.3,
		    double **weightMatrix = weights::IUB());

/**
 * Pair-wise align two amino acid sequences, using a modified
 * NeedleMan-Wunsh algorithm.
 *
 * The two sequences seq1 and seq2 are aligned in-place: gaps are inserted
 * according to a global alignment, and they will have equal length.
 *
 * The algorithm is NeedleMan-Wunsh, with two popular modifications:
 *  - there is a different cost for opening a gap or for extending a gap.
 *  - there is no gap open cost for a gap at the beginning or the end.
 */
extern double Align(AASequence& seq1, AASequence& seq2,
		    double gapOpenScore = -10, double gapExtensionScore = -3.3,
		    double **weightMatrix = weights::BLOSUM30());

/**
 * Error thrown by CodonAlign when apparent frame shifts cannot be corrected.
 *
 * Details in CodonAlign.
 */
class FrameShiftError : public std::exception
{
public:
  FrameShiftError(double ntScore, double codonScore,
		  const NTSequence& ntRef, const NTSequence& ntTarget);
  ~FrameShiftError() throw();

  const char *what() const throw() { return "Frameshift error"; }

  /** %Nucleotide alignment score.
   */
  double nucleotideAlignmentScore() const { return ntScore_; }

  /** Codon-based alignemnt score.
   */
  double codonAlignmentScore() const { return codonScore_; }

  /** %Nucleotide aligned reference sequence
   */
  const NTSequence& nucleotideAlignedRef() const { return ntRef_; }

  /** %Nucleotide aligned target sequence
   */
  const NTSequence& nucleotideAlignedTarget() const { return ntTarget_; }

private:
  double ntScore_, codonScore_;
  NTSequence ntRef_, ntTarget_;
};

/**
 * Perform codon-based alignment of nucleotide sequences.
 *
 * Two nucleotide sequences are pair-wise aligned, but so that gaps are
 * at codon boundaries. Optionally, frameshifts may be detected and corrected.
 *
 * The reference sequence must be of length a multiple of 3, and is assumed
 * to represent an Open Reading Frame (ORF).
 *
 * The procedure translates the target sequence in the 3 ORFs,
 * and for each ORF performs an amino-acid alignment against the translated
 * reference sequence. The best alignment is used to create the nucleotide
 * alignment.
 *
 * Then, the score of the codon aligned nucleotide alignment is computed, and
 * compared with a direct nucleotide alignment of both nucleotide sequences.
 * The codon alignment is accepted only if the difference is smaller than 100.
 * Otherwise, if maxFrameShifts > 0, the frameshift is searched, corrected
 * by inserting 1 or 2 'N' symbols in the target sequence, and repeating the
 * codon alignment. This is repeated for up to maxFrameShifts of times.
 *
 * The result is the nucleotide alignment score of the codon alignment, and
 * the number of frameshifts that have been corrected.
 *
 * @throws FrameShiftError when frameshifts could not be corrected, or
 *         the number of detected frameshifts exceeds maxFrameShifts.
 */
extern std::pair<double, int>
CodonAlign(NTSequence& ref, NTSequence& target,
	   int maxFrameShifts = 1,
	   double gapOpenScore = -10,
	   double gapExtensionScore = -3.3,
	   double **ntWeightMatrix = weights::IUB(),
	   double **aaWeightMatrix = weights::BLOSUM30());

};

#endif // ALIGN_H_
