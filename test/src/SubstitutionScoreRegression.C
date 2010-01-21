#include <iostream>

#include <AminoAcid.h>
#include <Nucleotide.h>
#include <NTSequence.h>
#include <AlignmentAlgorithm.h>
#include <Codon.h>

using namespace seq;

int main(int argc, char **argv)
{
	NTSequence cdn1(3);
	NTSequence cdn2(3);
	double ** submat = AlignmentAlgorithm::AmbiguousSubMat();

	std::cout << "n1,n2,n3,n4,n5,n6,score" << std::endl;
	for(int a = 1; a < 16; ++a){
		cdn1[0] = Nucleotide::fromRep(a);
		for(int b = 1; b < 16; ++b){
			cdn1[1] = Nucleotide::fromRep(b);
			for(int c = 1; c < 16; ++c){
				cdn1[2] = Nucleotide::fromRep(c);
				for(int d = 1; d < 16; ++d){
					cdn2[0] = Nucleotide::fromRep(d);
					for(int e = 1; e < 16; ++e){
						cdn2[1] = Nucleotide::fromRep(e);
						for(int f = 1; f < 16; ++f){
							cdn2[2] = Nucleotide::fromRep(f);

							AminoAcid aa1 = Codon::translate(cdn1.begin());
							AminoAcid aa2 = Codon::translate(cdn2.begin());
							if(aa1.intRep() > AminoAcid::AA_GAP || aa2.intRep() > AminoAcid::AA_GAP){
								double score = submat[aa1.intRep()][aa2.intRep()];
								std::cout << a << "," << b << "," << c << "," << d << "," << e << "," << f << "," << score << std::endl;
							}
						}
					}
				}
			}
		}
	}
	return 0;
}


