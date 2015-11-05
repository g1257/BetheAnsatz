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
#ifndef BETHE_PARAMETERS_TJ_H
#define BETHE_PARAMETERS_TJ_H
#include "../../Engine/ParametersBase.h"
#include "Vector.h"

namespace BetheAnsatz {

template<typename RealType_, typename InputType>
struct ParametersTj : public ParametersBase<RealType_, InputType> {

	typedef RealType_ RealType;

	ParametersTj(InputType& io)
	    : ParametersBase<RealType_, InputType>(io),
	      densityTotal(0),
	      densityBegin(0.0),
	      densityStep(0.0)
	{
		RealType densityEnd = 0;
		io.readline(densityBegin,"DensityBegin=");
		bool mode = false;
		try {
			io.readline(densityEnd,"DensityEnd=");
			mode = true;
		} catch (std::exception&) {}

		io.readline(densityTotal,"DensityTotal=");

		if (mode) {
			densityStep = (densityEnd-densityBegin)/densityTotal;
		} else {
			io.readline(densityStep,"DensityStep=");
		}
	}

	SizeType densityTotal;
	RealType densityBegin;
	RealType densityStep;
}; // struct ParametersTj

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const ParametersTj<T1,T2>& params)
{
	const ParametersBase<T1,T2>& parent = params;
	os<<parent;
	os<<"DensityTotal="<<params.densityTotal<<"\n";
	os<<"DensityBegin="<<params.densityBegin<<"\n";
	os<<"DensityStep="<<params.densityStep<<"\n";
	return os;
}
} // namespace BetheAnsatz
#endif // BETHE_PARAMETERS_TJ_H

