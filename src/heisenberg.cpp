/*
Copyright (c) 2015, UT-Battelle, LLC

BetheAnsatz, Version 0.1

This file is part of BetheAnsatz.
BetheAnsatz is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
BetheAnsatz is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with BetheAnsatz. If not, see <http://www.gnu.org/licenses/>.
*/
#include <cstdlib>
#include <unistd.h>
#define USE_PTHREADS_OR_NOT_NG
#include "Vector.h"
#include "InputNg.h"
#include "Engine/InputCheck.h"
#include "Models/Heisenberg/ParametersHeisenberg.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include "TypeToString.h"
#include "Models/Heisenberg/Heisenberg.h"

template<typename ParametersType>
class ParallelTemperature {

	typedef typename ParametersType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef BetheAnsatz::Heisenberg<ParametersType> HeisenbergType;
	typedef std::pair<RealType, RealType> PairRealType;
	typedef typename PsimagLite::Vector<PairRealType>::Type VectorPairRealType;

public:

	ParallelTemperature(const ParametersType& params)
	    : params_(params), omegaValue_(params.tt)
	{}

	SizeType tasks() const { return params_.tt; }

	void doTask(SizeType taskNumber ,SizeType threadNum)
	{
		RealType ts = (params_.te - params_.tb)/params_.tt;
		RealType t = params_.tb + taskNumber*ts;
		HeisenbergType heisenberg(params_,t);
		omegaValue_[taskNumber] = PairRealType(heisenberg.energy(),
		                              heisenberg.sz());

	}

	void print(std::ostream& os) const
	{
		os<<omegaValue_;
	}

private:

	const ParametersType& params_;
	VectorPairRealType omegaValue_;
}; // class ParallelTemperature

typedef double RealType;
typedef PsimagLite::InputNg<BetheAnsatz::InputCheck> InputNgType;

int main(int argc, char** argv)
{
	typedef BetheAnsatz::ParametersHeisenberg<RealType, InputNgType::Readable>
	        ParametersType;
	typedef ParallelTemperature<ParametersType> ParallelTemperatureType;
	typedef PsimagLite::Parallelizer<ParallelTemperatureType> ParallelizerType;

	int opt = 0;
	SizeType nthreads = 1;
	PsimagLite::String filename;
	BetheAnsatz::InputCheck inputCheck;
	PsimagLite::String usage(argv[0]);
	usage += " -f filename [-t threads]\n";

	while ((opt = getopt(argc, argv,"f:t:")) != -1) {
		switch (opt) {
		case 'f':
			filename = optarg;
			break;
		case 't':
			nthreads = atoi(optarg);
			break;
		default:
			inputCheck.usageMain(usage);
			return 1;
		}
	}

	if (filename == "") {
		inputCheck.usageMain(usage);
		return 2;
	}

	PsimagLite::Concurrency concurrency(&argc,&argv,nthreads);

	InputNgType::Writeable ioWriteable(filename,inputCheck);
	InputNgType::Readable io(ioWriteable);

	ParametersType params(io);
	std::cerr<<"Echo of Parameters read from "<<filename<<"\n";
	std::cerr<<params;
	std::cerr<<"Threads="<<PsimagLite::Concurrency::codeSectionParams.npthreads<<"\n";
	inputCheck.checkForThreads(PsimagLite::Concurrency::codeSectionParams.npthreads);

	ParallelizerType threadObject(PsimagLite::Concurrency::codeSectionParams);
	ParallelTemperatureType parallelTemperature(params);

	threadObject.loopCreate(parallelTemperature);

	parallelTemperature.print(std::cout);
}


