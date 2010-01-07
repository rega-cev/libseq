#ifndef AMBIGUOUSTABLE_H_
#define AMBIGUOUSTABLE_H_

#include "AminoAcid.h"
#include <map>
#include <set>

namespace seq {

class AmbiguousTable {
public:

	AmbiguousTable();
	int getIndex(std::set<AminoAcid> ambiguities) {
		return lookup_[ambiguities];
	}

private:

	std::map< std::set<AminoAcid> , int > lookup_;
};

}
#endif /* AMBIGUOUSTABLE_H_ */
