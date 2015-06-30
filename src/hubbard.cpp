#include <cstdlib>
#include <unistd.h>
#include "Vector.h"
#include "GrandPotential.h"

void printUsage(PsimagLite::String prog)
{
	std::cerr<<"USAGE: "<<prog<<" -b begin -t total -s step -m mu -U U\n";

}

typedef double RealType;

int main(int argc, char** argv)
{
	int opt = 0;
	SizeType meshTotal = 1000;
	SizeType tt = 0;
	RealType tb;
	RealType ts;
	RealType mu;
	RealType U;
	while ((opt = getopt(argc, argv,"b:t:s:m:U:M:")) != -1) {
		switch (opt) {
		case 'b':
			tb = atof(optarg);
			break;
		case 't':
			tt = atoi(optarg);
			break;
		case 's':
			ts = atof(optarg);
			break;
		case 'm':
			mu = atof(optarg);
			break;
		case 'U':
			U = atof(optarg);
			break;
		case 'M':
			meshTotal = atoi(optarg);
			break;
		default:
			printUsage(argv[0]);
			return 1;
		}
	}

	if (tt == 0 or ts <= 0) {
		printUsage(argv[0]);
		return 2;
	}

	BetheAnsatz::GrandPotential<RealType> grandPotential(mu,U,meshTotal);
	for (SizeType i = 0; i < tt; ++i) {
		RealType t = tb + i*ts;
		RealType omegaValue = grandPotential.at(t);
		std::cout<<t<<" "<<omegaValue<<"\n";
	}
}

