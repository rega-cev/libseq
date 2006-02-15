#include <boost/tokenizer.hpp>

#include "Mutation.h"

namespace seq {

std::set<AAMutation> readMutations(std::istream& mutationFile,
				   std::string prefix)
{
  std::set<AAMutation> result;

  std::string line;
  typedef boost::tokenizer<boost::escaped_list_separator<char> > csv_tok;
  getline(mutationFile, line);
  csv_tok tok(line);

  for (csv_tok::iterator i = tok.begin(); i != tok.end(); ++i) {
    std::string mutation = *i;
    AminoAcid aa(mutation[mutation.length() - 1]);
    int pos
      = atoi(mutation.substr(prefix.length(),
                             mutation.length() - prefix.length() - 1).c_str());

    result.insert(AAMutation(pos, aa));
  }

  return result;
}

};
