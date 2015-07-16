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
	typedef Rho<LogEtaType> RhoType;
	typedef typename ParametersType::RealType RealType;

public:

	Heisenberg(const ParametersType& params,
	           RealType temperature,
	           std::ostream& clog)
	    : logEta_(params,temperature,clog),rho_(logEta_)
	{}

private:

	LogEtaType logEta_;
	RhoType rho_;
}; // class Heisenberg
} // namespace BetheAnsatz

#endif // BETHE_HEISENBERG_H

