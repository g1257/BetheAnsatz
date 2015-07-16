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
#include "Vector.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "Parameters.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include "TypeToString.h"
#include "Heisenberg.h"

template<typename ParametersType>
class ParallelTemperature {

	typedef typename ParametersType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

public:

	ParallelTemperature(const ParametersType& params,
	                    SizeType threads)
	    : params_(params), omegaValue_(params.tt)
	{
		for (SizeType i = 0; i < threads; ++i) {
			PsimagLite::String name = filenameForThread(i);
			unlink(name.c_str());
		}
	}

	void thread_function_(SizeType threadNum,
	                      SizeType blockSize,
	                      SizeType total,
	                      PsimagLite::Concurrency::MutexType*)
	{
		RealType ts = (params_.te - params_.tb)/params_.tt;

		PsimagLite::String name = filenameForThread(threadNum);
		std::ofstream fout(name.c_str(),std::ios::app);
		for (SizeType p = 0; p < blockSize; p++) {
			SizeType i = threadNum*blockSize + p;
			if (i >= total) break;
			RealType t = params_.tb + i*ts;
			omegaValue_[i] = 0.0;
		}

		fout.close();
	}

	void print(std::ostream& os) const
	{
		os<<omegaValue_;
	}

private:

	PsimagLite::String filenameForThread(SizeType thread) const
	{
		return params_.logroot + ttos(thread) + ".txt";
	}

	const ParametersType& params_;
	VectorRealType omegaValue_;
}; // class ParallelTemperature

typedef double RealType;
typedef PsimagLite::InputNg<BetheAnsatz::InputCheck> InputNgType;

int main(int argc, char** argv)
{
	typedef BetheAnsatz::Parameters<RealType, InputNgType::Readable> ParametersType;
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
	std::cerr<<"Threads="<<PsimagLite::Concurrency::npthreads<<"\n";
	inputCheck.checkForThreads(PsimagLite::Concurrency::npthreads);

	ParallelizerType threadObject(PsimagLite::Concurrency::npthreads,
	                              PsimagLite::MPI::COMM_WORLD);
	ParallelTemperatureType parallelTemperature(params,nthreads);

	threadObject.loopCreate(params.tt,parallelTemperature);

	parallelTemperature.print(std::cout);
}


