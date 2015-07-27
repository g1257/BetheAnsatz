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
#ifndef BETHE_HEISENBERG_H
#define BETHE_HEISENBERG_H
#include "LogEta.h"
#include "Rho.h"

namespace BetheAnsatz {

template<typename ParametersType>
class Heisenberg {

	typedef LogEta<ParametersType> LogEtaType;
	typedef typename LogEtaType::MeshType MeshType;
	typedef Rho<LogEtaType> RhoType;
	typedef typename ParametersType::RealType RealType;

public:

	Heisenberg(const ParametersType& params,
	           RealType temperature,
	           std::ostream& clog)
	    : J_(params.J),
	      logEta_(params,temperature,clog),
	      mesh_(logEta_.mesh()),
	      rho_(params,temperature,clog,logEta_),
	      energy_(0.25*J_),
	      sz_(0.5)
	{
		//clog<<logEta_;
		//clog<<rho_;
		for (SizeType n = 0; n < params.nMax; ++n) {
			energy_ += integralEnergy(n);
			sz_ -= (n+1)*RhoType::integralSz(n,rho_.matrix(),mesh_);
		}
	}

	const RealType& energy() const { return energy_; }

	const RealType& sz() const { return sz_; }

private:

	RealType integralEnergy(SizeType n) const
	{
		RealType sum = 0;
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			RealType k = mesh_.x(i);
			sum += g(n+1,k)*rho_(n,i);
		}

		return sum*mesh_.step();
	}

	RealType g(SizeType n, RealType k) const
	{
		return -2.0*n*J_/(k*k + n*n);
	}

	RealType J_;
	LogEtaType logEta_;
	const MeshType& mesh_;
	RhoType rho_;
	RealType energy_;
	RealType sz_;
}; // class Heisenberg
} // namespace BetheAnsatz

#endif // BETHE_HEISENBERG_H

