#include "NTSequence.h"
#include "AASequence.h"

#include <iterator>
#include <fstream>
#include <math.h>
#include <limits>
#include <map>
#include <boost/lexical_cast.hpp>

double MU = 2.5E-5;

using namespace seq;

void readSequences(const char *fName, std::vector<NTSequence>& naive_seqs, std::vector<NTSequence>& treated_seqs)
{
  std::ifstream seqs(fName);

  /*
   * Iterate over all nucleotide sequences in the file.
   */
  try {
    std::cerr << "start reading sequences ..." << std::endl;

      int naive = 0;
      int treated = 0; 
    for (std::istream_iterator<NTSequence> i(seqs);
	 i != std::istream_iterator<NTSequence>();
	 ++i) {
      const NTSequence& s = *i;
      if (s.name()[0] == 'N') 
      {
      naive++;
      	naive_seqs.push_back(s);
      }
      else
      {
      	treated++;
      	treated_seqs.push_back(s);
      }
    }
    std::cerr << "naive sequences: " << naive << std::endl;
    std::cerr << "treated sequences: " << treated << std::endl;
  } catch (ParseException& e) {
    std::cerr << "Error reading " << fName << ": "
	      << e.message() << std::endl;
  }

}

double calculateGeneticDiversity(std::vector<NTSequence>& naive_seqs, std::vector<NTSequence>& treated_seqs)
{
    std::cerr << "calculating genetic diversity ..." << std::endl;
	
	long long includedCount = 0;
	long long diffCount = 0;
	std::cerr << "progress ";
	int seqLength = treated_seqs[0].size();
	for(int i = 0; i<naive_seqs.size(); i++)
	{
    		std::cerr << ".";
		for(int j = 0; j<treated_seqs.size(); j++)
		{
			for(int k = 0; k<seqLength; k++)
			{
				if(naive_seqs[i][k]!=Nucleotide::GAP && treated_seqs[j][k]!=Nucleotide::GAP)
				{
					includedCount++;
					if(naive_seqs[i][k]!=treated_seqs[j][k])
					{
						diffCount++;
					}
				}
			}
		}
	}
    		
	std::cerr << std::endl; 
	
    std::cerr << "included nucleotides: " << includedCount << std::endl;
    std::cerr << "different nucleotides: " << diffCount << std::endl;

    return (double)diffCount/includedCount;
}

int main(int argc, char **argv)
{
  if (argc < 1) {
    std::cerr << "Usage: " << std::endl
	      << argv[0] << " sequences.fasta" << std::endl;
    exit(1);
  }

  std::vector<NTSequence> naive_seqs;
  std::vector<NTSequence> treated_seqs;
  readSequences(argv[1], naive_seqs, treated_seqs);
  double geneticDiversity = calculateGeneticDiversity(naive_seqs, treated_seqs);
  std::cerr << "genetic diversity: " << geneticDiversity << std::endl;
}
