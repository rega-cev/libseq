#include <CodonAlign.h>
#include <NeedlemanWunsh.h>

#include <fstream>

using namespace seq;

int main(int argc, char **argv)
{
  std::ifstream s1(argv[1]);
  std::ifstream s2(argv[2]);

  NTSequence seq1, seq2;
  s1 >> seq1;
  s2 >> seq2;
  AlignmentAlgorithm* needlemanWunsh;
  try {
	needlemanWunsh = new NeedlemanWunsh(-10,-3.3);
	CodonAlign ca(needlemanWunsh);
    std::pair<double, int> result = ca.align(seq1, seq2, 2);
    std::cerr << "Aligned successfully, score: " << result.first
	      << ", frameshifts: " << result.second << std::endl;
    std::cerr << seq1 << seq2 << std::endl;
  } catch (std::exception& e) {
    std::cerr << "Alignment problem: " << e.what() << std::endl;
  }
  delete needlemanWunsh;
  return 0;
}
