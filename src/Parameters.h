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
#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "Vector.h"

namespace BetheAnsatz {

template<typename RealType_, typename InputType>
struct Parameters {

	typedef RealType_ RealType;

	Parameters(InputType& io)
	    : meshKtotal(1000),
	      meshLambdaTotal(10000),
	      nMax(50),
	      tt(0.0),
	      iterations(1),
	      infty(1e6),
	      tb(0.0),
	      te(0.0),
	      U(0.0),
	      logroot("ba")
	{
		io.readline(nMax,"Nmax=");
		io.readline(tt,"TemperatureTotal=");
		io.readline(mt,"MuTotal=");
		io.readline(iterations,"Iterations=");
		io.readline(tb,"TemperatureBegin=");
		io.readline(te,"TemperatureEnd=");
		io.readline(mb,"MuBegin=");
		io.readline(me,"MuEnd=");
		io.readline(U,"U=");

		try {
			io.readline(meshKtotal,"MeshKTotal=");
		} catch (std::exception&) {}

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

	SizeType meshKtotal;
	SizeType meshLambdaTotal;
	SizeType nMax;
	SizeType tt;
	SizeType mt;
	SizeType iterations;
	RealType infty;
	RealType tb;
	RealType te;
	RealType mb;
	RealType me;
	RealType U;
	PsimagLite::String logroot;
}; // struct Parameters

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const Parameters<T1,T2>& params)
{
	os<<"MeshKTotal="<<params.meshKtotal<<"\n";
	os<<"MeshLambdaTotal="<<params.meshLambdaTotal<<"\n";
	os<<"Nmax="<<params.nMax<<"\n";
	os<<"TemperatureTotal="<<params.tt<<"\n";
	os<<"MuTotal="<<params.mt<<"\n";
	os<<"Iterations="<<params.iterations<<"\n";
	os<<"Infty="<<params.infty<<"\n";
	os<<"TemperatureBegin="<<params.tb<<"\n";
	os<<"TemperatureEnd="<<params.te<<"\n";
	os<<"MuBegin="<<params.mb<<"\n";
	os<<"MuEnd="<<params.me<<"\n";
	os<<"U="<<params.U<<"\n";
	os<<"LogRoot="<<params.logroot<<"\n";

	return os;
}
} // namespace BetheAnsatz
#endif // PARAMETERS_H

