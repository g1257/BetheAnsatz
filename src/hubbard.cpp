#include <cstdlib>
#include <unistd.h>
#include "Vector.h"
#include "Grounded.h"
#include "GrandPotential.h"
#include "InputNg.h"
#include "InputCheck.h"
#include "Parameters.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include "TypeToString.h"

template<typename ParametersType, typename GroundedType>
class ParallelTemperature {

	typedef typename ParametersType::RealType RealType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;

public:

	ParallelTemperature(const ParametersType& params,
	                    const GroundedType& grounded,
	                    SizeType threads)
	    : params_(params), grounded_(grounded), omegaValue_(params.tt, params.mt)
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
		RealType ms = (params_.me - params_.mb)/params_.mt;

		PsimagLite::String name = filenameForThread(threadNum);
		std::ofstream fout(name.c_str(),std::ios::app);
		for (SizeType p = 0; p < blockSize; p++) {
			SizeType i = threadNum*blockSize + p;
			if (i >= total) break;
			RealType t = params_.tb + i*ts;
			for (SizeType j = 0; j < params_.mt; ++j) {
				RealType mu = params_.mb + i*ms;
				BetheAnsatz::GrandPotential<ParametersType> grandPotential(params_,
				                                                           grounded_,
				                                                           mu,
				                                                           t,
				                                                           fout);
				omegaValue_(i,j) = grandPotential();
			}

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
	const GroundedType& grounded_;
	MatrixRealType omegaValue_;
}; // class ParallelTemperature

typedef double RealType;
typedef PsimagLite::InputNg<BetheAnsatz::InputCheck> InputNgType;

int main(int argc, char** argv)
{
	typedef BetheAnsatz::Parameters<RealType, InputNgType::Readable> ParametersType;
	typedef BetheAnsatz::Grounded<ParametersType> GroundedType;
	typedef ParallelTemperature<ParametersType, GroundedType> ParallelTemperatureType;
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

	GroundedType grounded(params);

	ParallelizerType threadObject(PsimagLite::Concurrency::npthreads,
	                              PsimagLite::MPI::COMM_WORLD);
	ParallelTemperatureType parallelTemperature(params,grounded,nthreads);

	threadObject.loopCreate(params.tt,parallelTemperature);

	parallelTemperature.print(std::cout);
}


