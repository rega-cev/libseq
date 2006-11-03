#include "Align.h"

#include <fstream>

using namespace seq;

int main(int argc, char **argv)
{
  std::ifstream s1(argv[1]);
  std::ifstream s2(argv[2]);

  NTSequence seq1, seq2;
  s1 >> seq1;
  s2 >> seq2;

  try {
    std::pair<double, int> result = CodonAlign(seq1, seq2, 2);
    std::cerr << "Aligned successfully, score: " << result.first
	      << ", frameshifts: " << result.second << std::endl;

    std::cerr << seq1 << seq2 << std::endl;
  } catch (std::exception& e) {
    std::cerr << "Alignment problem: " << e.what() << std::endl;
  }

  return 0;
}
