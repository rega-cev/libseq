#include <NTSequence.h>
#include <set>
#include <iostream>
#include <fstream>
#include <stdlib.h>

int main(int argc, char **argv)
{
  if (argc < 3) {
    std::cerr << "Usage: stockholm sequences.fasta sequences.stockholm [line-length]"
	      << std::endl;
    return 1;
  }
  int length = 10000;
  if(argc == 4){
    length = atoi(argv[3]);
  }

  std::vector<seq::NTSequence> sequences;

  std::ifstream f_seqs(argv[1]);
  while (f_seqs) {
    seq::NTSequence s;

    try {
      f_seqs >> s;

      if (f_seqs) {
	sequences.push_back(s);
      }
    } catch (seq::ParseException& e) {
      if (e.recovered())
	std::cout << e.name() << ",\"" << e.message() << "\"" << std::endl;
      else
	std::cerr << "Fatal error: " << e.message() << std::endl;

      if (!e.recovered())
	return 1;
    }
  }

  std::ofstream f_out(argv[2]);

  seq::writeStockholm(f_out,sequences,length);

  return 0;
}

