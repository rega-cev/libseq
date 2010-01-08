#include "NeedlemanWunsh.h"
#include <fstream>

using namespace seq;

int main(int argc, char **argv) {
	std::ifstream s1(argv[1]);
	std::ifstream s2(argv[2]);

	NTSequence seq1, seq2;
	s1 >> seq1;
	s2 >> seq2;

	NeedlemanWunsh nmw;
	AASequence aaseq1 = AASequence::translate(seq1);
	AASequence aaseq2 = AASequence::translate(seq2);
	double result = nmw.align(aaseq1,aaseq2);

	std::cerr << std::endl;
	std::cerr << aaseq1 << std::endl;
	std::cerr << aaseq2 << std::endl;
	return 0;
}

