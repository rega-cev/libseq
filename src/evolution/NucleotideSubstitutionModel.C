#include "NucleotideSubstitutionModel.h"

namespace seq {

NucleotideSubstitutionModel::NucleotideSubstitutionModel
   (double piA, double piC, double piG, double piT,
    double rAC, double rAG, double rAT,
    double rCG, double rCT, double rGT,
    double errorRate)
{ 
  /*
   * mu = - piA . QAA - piC . QCC - piG . QGG - piT . QTT
   *    = - C . (piA . .... ).
   *
   * thus we can simply rescale all matrix entries with a same
   * factor to obtain the expected mu.
   */

  double rCA = rAC, rGA = rAG, rTA = rAT, rGC = rCG, rTC = rCT, rTG = rGT;

  matrix_[Nucleotide::NT_A][Nucleotide::NT_C] = rAC * piC;
  matrix_[Nucleotide::NT_A][Nucleotide::NT_G] = rAG * piG;
  matrix_[Nucleotide::NT_A][Nucleotide::NT_T] = rAT * piT;
  matrix_[Nucleotide::NT_A][Nucleotide::NT_A] =
    - (matrix_[Nucleotide::NT_A][Nucleotide::NT_C]
       + matrix_[Nucleotide::NT_A][Nucleotide::NT_G]
       + matrix_[Nucleotide::NT_A][Nucleotide::NT_T]);

  matrix_[Nucleotide::NT_C][Nucleotide::NT_A] = rCA * piA;
  matrix_[Nucleotide::NT_C][Nucleotide::NT_G] = rCG * piG;
  matrix_[Nucleotide::NT_C][Nucleotide::NT_T] = rCT * piT;
  matrix_[Nucleotide::NT_C][Nucleotide::NT_C] =
    - (matrix_[Nucleotide::NT_C][Nucleotide::NT_A]
       + matrix_[Nucleotide::NT_C][Nucleotide::NT_G]
       + matrix_[Nucleotide::NT_C][Nucleotide::NT_T]);
  
  matrix_[Nucleotide::NT_G][Nucleotide::NT_A] = rGA * piA;
  matrix_[Nucleotide::NT_G][Nucleotide::NT_C] = rGC * piC;
  matrix_[Nucleotide::NT_G][Nucleotide::NT_T] = rGT * piT;
  matrix_[Nucleotide::NT_G][Nucleotide::NT_G] =
    - (matrix_[Nucleotide::NT_G][Nucleotide::NT_A]
       + matrix_[Nucleotide::NT_G][Nucleotide::NT_C]
       + matrix_[Nucleotide::NT_G][Nucleotide::NT_T]);

  matrix_[Nucleotide::NT_T][Nucleotide::NT_A] = rTA * piA;
  matrix_[Nucleotide::NT_T][Nucleotide::NT_C] = rTC * piC;
  matrix_[Nucleotide::NT_T][Nucleotide::NT_G] = rTG * piG;
  matrix_[Nucleotide::NT_T][Nucleotide::NT_T] =
    - (matrix_[Nucleotide::NT_T][Nucleotide::NT_A]
       + matrix_[Nucleotide::NT_T][Nucleotide::NT_C]
       + matrix_[Nucleotide::NT_T][Nucleotide::NT_G]);

  double mu = - piA * matrix_[Nucleotide::NT_A][Nucleotide::NT_A]
    - piC * matrix_[Nucleotide::NT_C][Nucleotide::NT_C]
    - piG * matrix_[Nucleotide::NT_G][Nucleotide::NT_G]
    - piT * matrix_[Nucleotide::NT_T][Nucleotide::NT_T];

  double C = errorRate / mu;

  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      matrix_[i][j] *= C;
}

double NucleotideSubstitutionModel::getMu(Nucleotide fromNT,
					  Nucleotide toNT) const
{
  return matrix_[fromNT.intRep()][toNT.intRep()];
}

};
