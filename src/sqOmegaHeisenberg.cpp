#include "Models/Heisenberg/TwoSpinonHeisenberg.h"

void usageMain(const char* name)
{
	std::cout<<"USAGE: "<<name<<" [-w stepEnergy] [-W energyPoints] ";
	std::cout<<"[-k stepK] [-K kPoints]\n";
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
	TwoSpinonHeisenberg<double> twoSpinon(kmesh,emesh);
	twoSpinon.toGnuplot(std::cout,cutoff);
}

