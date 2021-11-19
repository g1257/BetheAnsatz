#define  USE_PTHREADS_OR_NOT_NG
#include "Models/Heisenberg/TwoSpinonHeisenberg.h"
#include "Models/Heisenberg/FourSpinonHeisenberg.h"

void usageMain(const char* name)
{
	std::cout<<"USAGE: "<<name<<" [-w stepEnergy] [-W energyPoints] ";
	std::cout<<"[-k stepK] [-K kPoints] what\n";
}

int main(int argc, char** argv)
{
	SizeType total = 31;
	double step = 0.1;
	SizeType totalE = 10;
	double stepE = 0.1;
	int opt = 0;
	double cutoff = 0;
	while ((opt = getopt(argc, argv,"w:W:k:K:C:")) != -1) {
		switch (opt) {
		case 'w':
			stepE = atof(optarg);
			break;
		case 'W':
			totalE = atoi(optarg);
			break;
		case 'k':
			step = atof(optarg);
			break;
		case 'K':
			total = atoi(optarg);
			break;
		case 'C':
			cutoff = atof(optarg);
			break;
		default:
			usageMain(argv[0]);
			return 1;
		}
	}

	BetheAnsatz::Mesh<double> kmesh(total,0.0,step);
	BetheAnsatz::Mesh<double> emesh(totalE,0.0,stepE);

	PsimagLite::String what = (optind < argc) ? argv[optind] : "";

	if (what == "") return 1;

	if (what.find("two") != PsimagLite::String::npos) {
		BetheAnsatz::TwoSpinonHeisenberg<double> twoSpinon(kmesh,emesh);
		twoSpinon.toGnuplot(std::cout,cutoff);
	}

	if (what.find("four") != PsimagLite::String::npos) {
		BetheAnsatz::FourSpinonHeisenberg<double> fourSpinon(kmesh,emesh);
		fourSpinon.toGnuplot(std::cout,cutoff);
	}
}

