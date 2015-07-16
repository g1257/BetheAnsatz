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
#ifndef BETHE_PARAMETERS_HEISENBERG_H
#define BETHE_PARAMETERS_HEISENBERG_H
#include "ParametersBase.h"
#include "Vector.h"

namespace BetheAnsatz {

template<typename RealType_, typename InputType>
struct ParametersHeisenberg : public ParametersBase<RealType_, InputType> {

	typedef RealType_ RealType;

	ParametersHeisenberg(InputType& io)
	    : ParametersBase<RealType_, InputType>(io)
	{
		io.readline(J,"J=");
	}

	RealType J;
}; // struct ParametersHeisenberg

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const ParametersHeisenberg<T1,T2>& params)
{
	const ParametersBase<T1,T2>& parent = params;
	os<<parent;
	os<<"J="<<params.J<<"\n";

	return os;
}
} // namespace BetheAnsatz
#endif // BETHE_PARAMETERS_HEISENBERG_H

