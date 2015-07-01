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

