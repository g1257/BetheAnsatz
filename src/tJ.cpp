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
#include "Engine/InputCheck.h"
#include "Concurrency.h"
#include "Parallelizer.h"
#include "TypeToString.h"
#include "Models/Tj/ParametersTj.h"

typedef double RealType;
typedef PsimagLite::InputNg<BetheAnsatz::InputCheck> InputNgType;

int main(int argc, char** argv)
{
	typedef BetheAnsatz::ParametersTj<RealType, InputNgType::Readable> ParametersType;

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

}


