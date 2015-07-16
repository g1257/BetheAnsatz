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
#ifndef BETHE_PARAMETERS_HUBBARD_H
#define BETHE_PARAMETERS_HUBBARD_H
#include "ParametersBase.h"
#include "Vector.h"

namespace BetheAnsatz {

template<typename RealType_, typename InputType>
struct ParametersHubbard : public ParametersBase<RealType_, InputType> {

	typedef RealType_ RealType;

	ParametersHubbard(InputType& io)
	    : ParametersBase<RealType_, InputType>(io),
	      meshKtotal(1000),
	      U(0.0)
	{
		io.readline(mt,"MuTotal=");
		io.readline(mb,"MuBegin=");
		io.readline(me,"MuEnd=");
		io.readline(U,"U=");

		try {
			io.readline(meshKtotal,"MeshKTotal=");
		} catch (std::exception&) {}
	}

	SizeType meshKtotal;
	SizeType mt;
	RealType mb;
	RealType me;
	RealType U;
}; // struct ParametersHubbard

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const ParametersHubbard<T1,T2>& params)
{
	const ParametersBase<T1,T2>& parent = params;
	os<<parent;
	os<<"MeshKTotal="<<params.meshKtotal<<"\n";
	os<<"MuTotal="<<params.mt<<"\n";
	os<<"MuBegin="<<params.mb<<"\n";
	os<<"MuEnd="<<params.me<<"\n";
	os<<"U="<<params.U<<"\n";

	return os;
}
} // namespace BetheAnsatz
#endif // BETHE_PARAMETERS_HUBBARD_H

