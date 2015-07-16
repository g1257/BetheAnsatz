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
#include "Mesh.h"

namespace BetheAnsatz {

template<typename ParametersType>
class Heisenberg {

	typedef typename ParametersType::RealType RealType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef Mesh<RealType> MeshType;

	Heisenberg(const ParametersType& params)
	    : mesh_(2*params.infty,-params.infty,2.0*params.infty/params.meshLambdaTotal),
	      eta_(params.nMax,mesh_.total()),
	      minusTwoJOverT_(-2.0*params.J/params.T)
	{}

private:

	void calcEta1()
	{
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			RealType k = mesh_.x(k);
			RealType logn1k = minusTwoJOverT_/(1.0 + k*k) + conv1(i) - conv2(i);
			eta_(0,i) = exp(logn1k);
		}
	}

	void calcEtaN(SizeType n)
	{
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			RealType k = mesh_.x(k);
		}
	}

	RealType conv1(SizeType i) const
	{
		return 0.0;
	}

	RealType conv2(SizeType i) const
	{
		return 0.0;
	}

	MeshType mesh_;
	MatrixRealType eta_;
	RealType minusTwoJOverT_;
}; // class Heisenberg
} // namespace BetheAnsatz

#endif // BETHE_HEISENBERG_H

