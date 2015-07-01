#ifndef GRANDPOTENTIAL_H
#define GRANDPOTENTIAL_H
#include "Grounded.h"

namespace BetheAnsatz {

template<typename RealType>
class GrandPotential {

public:

	typedef Grounded<RealType> GroundedType;

	GrandPotential(RealType, RealType U, SizeType meshTotal)
	    : grounded_(U, meshTotal)
	{}

	RealType at(RealType)
	{
		return 0.0;
	}

	// for testing
	RealType rho0(SizeType i) const
	{
		return grounded_.rho0(i);
	}

private:

	GroundedType grounded_;
}; // class GrandPotential

} // namespace BetheAnsatz
#endif // GRANDPOTENTIAL_H

