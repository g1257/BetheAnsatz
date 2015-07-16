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
#ifndef BETHE_PARAMETERS_BASE_H
#define BETHE_PARAMETERS_BASE_H
#include "Vector.h"

namespace BetheAnsatz {

template<typename RealType_, typename InputType>
struct ParametersBase {

	typedef RealType_ RealType;

	ParametersBase(InputType& io)
	    : meshLambdaTotal(10000),
	      nMax(50),
	      tt(0.0),
	      iterations(1),
	      infty(1e6),
	      tb(0.0),
	      te(0.0),
	      logroot("ba")
	{
		io.readline(nMax,"Nmax=");
		io.readline(tt,"TemperatureTotal=");
		io.readline(iterations,"Iterations=");
		io.readline(tb,"TemperatureBegin=");
		io.readline(te,"TemperatureEnd=");

		try {
			io.readline(meshLambdaTotal,"MeshLambdaTotal=");
		} catch (std::exception&) {}

		try {
			io.readline(infty,"Infty=");
		} catch (std::exception&) {}

		try {
			io.readline(logroot,"LogRoot=");
		} catch (std::exception&) {}
	}

	SizeType meshLambdaTotal;
	SizeType nMax;
	SizeType tt;
	SizeType iterations;
	RealType infty;
	RealType tb;
	RealType te;
	PsimagLite::String logroot;
}; // struct ParametersBase

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const ParametersBase<T1,T2>& params)
{
	os<<"MeshLambdaTotal="<<params.meshLambdaTotal<<"\n";
	os<<"Nmax="<<params.nMax<<"\n";
	os<<"TemperatureTotal="<<params.tt<<"\n";
	os<<"Iterations="<<params.iterations<<"\n";
	os<<"Infty="<<params.infty<<"\n";
	os<<"TemperatureBegin="<<params.tb<<"\n";
	os<<"TemperatureEnd="<<params.te<<"\n";
	os<<"LogRoot="<<params.logroot<<"\n";

	return os;
}
} // namespace BetheAnsatz
#endif // BETHE_PARAMETERS_BASE_H

