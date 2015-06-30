#ifndef MESH_H
#define MESH_H
#include "Vector.h"
#include <cassert>

namespace BetheAnsatz {

template<typename RealType>
class Mesh {

public:

	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	Mesh(SizeType total, RealType start, RealType step) : x_(total), y_(total,0)
	{
		for (SizeType i = 0; i < total; ++i) {
			x_[i] = start + i*step;
		}
	}

	SizeType total() const { return x_.size(); }

	const RealType& x(SizeType i) const
	{
		assert(i < x_.size());
		return x_[i];
	}

	const RealType& operator[](SizeType i) const
	{
		assert(i < x_.size());
		return y_[i];
	}

	RealType& operator[](SizeType i)
	{
		assert(i < x_.size());
		return y_[i];
	}

private:

	VectorRealType x_;
	VectorRealType y_;
}; // class GrandPotential

} // namespace BetheAnsatz
#endif // MESH_H

