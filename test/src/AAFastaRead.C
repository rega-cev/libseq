#include "NTSequence.h"
#include "AASequence.h"

#include <iterator>
#include <fstream>

using namespace seq;

int main(int argc, char **argv)
{
  std::ifstream s(argv[1]);

  /*
   * Iterate over all amino acid sequences in the file.
   */
  try {
    for (std::istream_iterator<AASequence> i(s);
	 i != std::istream_iterator<AASequence>();
	 ++i) {
      const AASequence& seq = *i;

      // write out a FASTA file entry
      std::cout << seq;
    }
  } catch (ParseException& e) {
    std::cerr << "Error reading " << argv[1] << ": "
	      << e.message() << std::endl;
  }
}
