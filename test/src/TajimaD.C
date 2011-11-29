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

void readSequences(const char *fName, std::vector<NTSequence>& result,
		   int from, int to, std::string group)
{
  std::ifstream seqs(fName);

  /*
   * Iterate over all nucleotide sequences in the file.
   */
  try {
    for (std::istream_iterator<NTSequence> i(seqs);
	 i != std::istream_iterator<NTSequence>();
	 ++i) {
      const NTSequence& s = *i;
      if (s.name().find(group) != std::string::npos) {

	if (from >= to || (to > s.size())) {
	  std::cerr << "Error: need from < to < "
		    << s.size() << std::endl;
	  exit(1);
	}

	NTSequence part(s.begin() + from - 1, s.begin() + to);
	part.setName(s.name());
	part.setDescription(s.description());

	result.push_back(part);
      }
    }
  } catch (ParseException& e) {
    std::cerr << "Error reading " << fName << ": "
	      << e.message() << std::endl;
  }

}

const char* conclusions[] = { "negative or selective sweeps",
			      "diversifying selection",
			      "neutral",
			      "unknown" };

void computeTajimaD(const std::vector<NTSequence>& sequences,
		    int& n, double& S, double& khat, double& theta,
		    double& theta2, double& D, int& conclusion)
{
  n = sequences.size();

  if (n < 4) {
    std::cerr << "Error: need at least 4 sequences" << std::endl;
    exit(1);
  }

  /*
   * compute S: number of segregating sites
   */

  S = 0;
  for (int i = 0; i < sequences[0].size(); ++i) {
    for (int j = 1; j < sequences.size(); ++j) {
      if (sequences[j][i] != sequences[0][i]) {
	S += 1.0;
	break;
      }
    }
  }

  /*
   * compute khat: number of pairwise nucleotide differences
   */

  khat = 0;
  for (int i = 0; i < sequences[0].size(); ++i) {
    for (int j = 0; j < sequences.size(); ++j) {
      for (int k = j+1; k < sequences.size(); ++k) {
	if (sequences[j][i] != sequences[k][i])
	  ++khat;
      }
    }
  }

  khat /= (n * (n-1) / 2);

  theta = khat / sequences[0].size();

  /*
   * compute D
   */
  double a1 = 0;
  double a2 = 0;

  for (int i = 1; i <= n - 1; ++i) {
    a1 += 1.0 / i;
    a2 += 1.0 / (i * i);
  }

  theta2 = S / a1 / sequences[0].size();

  double b1 = (n + 1.0)/(3.0 * (n - 1.0));
  double b2 = 2*(n*n + n + 3.0)/(9.0 * n * (n - 1.0));

  double c1 = b1 - 1.0 / a1;
  double c2 = b2 - (n + 2.0)/(a1 * n) + a2 / (a1 * a1);

  double e1 = c1 / a1;
  double e2 = c2 / (a1*a1 + a2);

  double d = (khat - S / a1);
  double Vd = e1 * S + e2 * S * (S - 1);

  if (Vd <= 0)
    D = 0;
  else
    D = d / sqrt(Vd);

  double betalim[][2] = { { -0.876, 2.232 },
			  { -1.269, 1.834 },
			  { -1.478, 1.999 },
			  { -1.608, 1.932 },
			  { -1.663, 1.975 },
			  { -1.713, 1.954 },
			  { -1.733, 1.975 },
			  { -1.757, 1.966 },
			  { -1.765, 1.979 },
			  { -1.779, 1.976 },
			  { -1.783, 1.985 },
			  { -1.791, 1.984 },
			  { -1.793, 1.990 },
			  { -1.798, 1.990 },
			  { -1.799, 1.996 },
			  { -1.802, 1.996 },
			  { -1.803, 2.001 },
			  { -1.805, 2.001 },
			  { -1.804, 2.005 } 
  };

  if (n >= 4 && n <= 22) {
    if (D < betalim[n - 4][0])
      conclusion = 0;
    else if (D > betalim[n - 4][1])
      conclusion = 1;
    else
      conclusion = 2;
  } else
    conclusion = 3;
}

void exportInfile(const std::vector<NTSequence>& sequences)
{
  std::ofstream f("infile");

  f << "1" << std::endl
    << sequences.size() << " " << sequences[0].size() << std::endl;

  for (int i = 0; i < sequences.size(); ++i) {
    f << "seq" << i << "        " << std::endl;
    f << sequences[i].asString() << std::endl;
  }
}

std::string summary(std::vector<double> &v, int n)
{
  if (n == 1)
    return boost::lexical_cast<std::string>(v[0]);

  double sum = 0;
  double min = std::numeric_limits<double>::max();
  double max = -std::numeric_limits<double>::max();
  for (int i = 0; i < n; ++i) {
    if (v[i] < min)
      min = v[i];
    if (v[i] > max)
      max = v[i];
    sum += v[i];
  }

  double average = sum / n;
  double var = 0;

  for (int i = 0; i < n; ++i)
    var += (v[i] - average)*(v[i] - average);

  double stdev = sqrt(var / n);

  return boost::lexical_cast<std::string>(average)
    + "," + boost::lexical_cast<std::string>(stdev)
    + "," + boost::lexical_cast<std::string>(min)
    + "," + boost::lexical_cast<std::string>(max);
}


std::string summary(std::vector<int> &v, int n)
{
  if (n == 1)
    return conclusions[v[0]];

  std::map<int, int> counts;
  for (int i = 0; i < n; ++i)
    ++counts[v[i]];

  int max_v;
  int max_count = std::numeric_limits<int>::min();
  for (std::map<int, int>::const_iterator i = counts.begin();
       i != counts.end(); ++i) {
    if (i->second > max_count) {
      max_count = i->second;
      max_v = i->first;
    }
  }

  return std::string(conclusions[max_v])
    + "," + boost::lexical_cast<std::string>(max_count);
}

int main(int argc, char **argv)
{
  if (argc < 5) {
    std::cerr << "Usage: " << std::endl
	      << argv[0] << " sequences.fasta from.nt to.nt "
	      << "group [replicates]" << std::endl;
    exit(1);
  }

  int replicates = (argc == 6 ? atoi(argv[5]) : 0);

  std::vector<NTSequence> allsequences;
  readSequences(argv[1], allsequences, atoi(argv[2]), atoi(argv[3]),
		argv[4]);

  if (replicates) {
    if (allsequences.size() % replicates != 0) {
      std::cerr << "Error: replicates needs to divide total sample."
		<< std::endl;
      exit(1);
    }

    int n = allsequences.size() / replicates;

    std::vector<int> partn(n);
	std::vector<int> conclusion(n);
    std::vector<double> S(n);
	std::vector<double> khat(n);
	std::vector<double> theta(n);
	std::vector<double> theta2(n);
	std::vector<double> D(n);

    for (int i = 0; i < replicates; ++i) {
      std::vector<NTSequence> sequences(allsequences.begin() + i*n,
					allsequences.begin() + (i+1)*n);

      computeTajimaD(sequences,
		     partn[i], S[i], khat[i], theta[i], theta2[i], D[i],
		     conclusion[i]);

      if (i == 0)
	exportInfile(sequences);
    }

    std::cout << n
	      << "," << replicates
	      << "," << allsequences[0].size()
	      << "," << summary(S, replicates)
	      << "," << summary(khat, replicates)
	      << "," << summary(theta, replicates)
	      << "," << summary(theta2, replicates)
	      << "," << summary(D, replicates)
	      << "," << summary(conclusion, replicates) << std::endl;
  } else {
    int n, conclusion;
    double S, khat, theta, theta2, D;

    computeTajimaD(allsequences,
		   n, S, khat, theta, theta2, D, conclusion);

    std::cout << "sample size: " << n << std::endl
	      << "sites: " << allsequences[0].size() << std::endl
	      << "number of segregating sites: " << S << std::endl
	      << "average nucleotide diversity: " << khat << std::endl
	      << "Theta (based on diversity): " << theta << std::endl
	      << " Ne (mu = " << MU << "): " << (theta / 2 / MU) << std::endl
	      << "Theta (based on segregating sites: " << theta2 << std::endl
	      << " Ne (mu = " << MU << "): " << (theta2 / 2 / MU) << std::endl
	      << "D = " << D << " " << conclusions[conclusion] << std::endl;

    exportInfile(allsequences);
  }
}
