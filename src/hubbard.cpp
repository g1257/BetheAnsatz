#include <cstdlib>
#include <unistd.h>
#include "Vector.h"
#include "Grounded.h"
#include "GrandPotential.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "Parameters.h"

void printUsage(PsimagLite::String prog)
{
	std::cerr<<"USAGE: "<<prog<<" -f filename\n";

}

typedef double RealType;
typedef PsimagLite::InputNg<BetheAnsatz::InputCheck> InputNgType;

int main(int argc, char** argv)
{
	int opt = 0;
	PsimagLite::String filename;
	BetheAnsatz::InputCheck inputCheck;

	while ((opt = getopt(argc, argv,"f:")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		default:
			printUsage(argv[0]);
			return 1;
		}
	}

	if (filename == "") {
		inputCheck.usageMain(PsimagLite::String(argv[0]) + " -f filename\n");
		return 2;
	}

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	typedef BetheAnsatz::Parameters<RealType, InputNgType::Readable> ParametersType;
	ParametersType params(io);
	std::cerr<<"Echo of Parameters read from "<<filename<<"\n";
	std::cerr<<params;
	BetheAnsatz::Grounded<ParametersType> grounded(params);

	RealType ts = (params.te - params.tb)/params.tt;
	RealType ms = (params.me - params.mb)/params.mt;
	for (SizeType i = 0; i < params.tt; ++i) {
		RealType t = params.tb + i*ts;
		for (SizeType j = 0; j < params.mt; ++j) {
			RealType mu = params.mb + i*ms;
			BetheAnsatz::GrandPotential<ParametersType> grandPotential(params,
			                                                           grounded,
			                                                           mu,
			                                                           t);
			RealType omegaValue = grandPotential();
			std::cout<<omegaValue<<" ";
		}

		std::cout<<"\n";
	}
}

