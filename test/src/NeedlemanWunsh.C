#include "NeedlemanWunsh.h"

#include <fstream>

using namespace seq;

int main(int argc, char **argv)
{
  std::ifstream s1(argv[1]);
  std::ifstream s2(argv[2]);

  AASequence seq1, seq2;
  s1 >> seq1;
  s2 >> seq2;
  std::cerr << seq1 << std::endl;
  std::cerr << seq2 << std::endl;
  NeedlemanWunsh needlemanWunsh(-10, -3.3);
  double result = needlemanWunsh.align(seq1, seq2);
  std::cerr << std::endl;
  std::cerr << seq1 << std::endl;
  std::cerr << seq2 << std::endl;
  return 0;
}
