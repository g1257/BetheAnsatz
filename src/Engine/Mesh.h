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
#ifndef MESH_H
#define MESH_H
#include "Vector.h"
#include <cassert>

namespace BetheAnsatz {

template<typename RealType>
class Mesh {

public:

	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	Mesh(SizeType total, RealType start, RealType step)
	    : step_(step), data_(total)
	{
		for (SizeType i = 0; i < total; ++i) {
			data_[i] = start + i*step;
		}
	}

	SizeType total() const { return data_.size(); }

	RealType step() const { return step_; }

	const RealType& x(SizeType i) const
	{
		assert(i < data_.size());
		return data_[i];
	}

private:

	RealType step_;
	VectorRealType data_;
}; // class GrandPotential

} // namespace BetheAnsatz
#endif // MESH_H

