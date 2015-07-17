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

#ifndef BETHE_CONVOLUTION_H
#define BETHE_CONVOLUTION_H
#include "Vector.h"

namespace BetheAnsatz {

template<typename BracketableFunctionType,
         typename AuxiliaryType,
         typename MeshType>
class Convolution {

	typedef typename MeshType::RealType RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

public:

	Convolution(SizeType n,
	            BracketableFunctionType& f,
	            const AuxiliaryType& a)
	    : result_(a.mesh.total()),mesh_(a.mesh)
	{

		for (SizeType i = 0; i < mesh_.total(); ++i) {
			RealType k = mesh_.x(i);
			RealType sum = 0.0;
			for (SizeType j = 0; j < mesh_.total(); ++j) {
				RealType kdiff = k - mesh_.x(j);
				RealType tmp = n*n + kdiff*kdiff;
				sum += n*f(j,a)/tmp;
			}

			assert(i < result_.size());
			result_[i] = sum*mesh_.step()/M_PI;
		}
	}

	const RealType& operator()(SizeType i) const
	{
		assert(i < result_.size());
		return result_[i];
	}

private:

	VectorRealType result_;
	const MeshType& mesh_;
}; // class Convolution

} // namespace BetheAnsatz

#endif // BETHE_CONVOLUTION_H

