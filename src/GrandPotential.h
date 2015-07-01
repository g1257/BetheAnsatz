#ifndef GRANDPOTENTIAL_H
#define GRANDPOTENTIAL_H
#include "Grounded.h"

namespace BetheAnsatz {

template<typename RealType>
class GrandPotential {

public:

	typedef Grounded<RealType> GroundedType;

	GrandPotential(const GroundedType& grounded,
	               RealType mu,
	               RealType T,
	               SizeType meshLambdaTotal)
	    : grounded_(grounded),mu_(mu),T_(T),meshLambdaTotal_(meshLambdaTotal)
	{}

	RealType operator()() const
	{
		return 0.0;
	}

private:

	const GroundedType& grounded_;
	RealType mu_;
    RealType T_;
    SizeType meshLambdaTotal_;
}; // class GrandPotential

} // namespace BetheAnsatz
#endif // GRANDPOTENTIAL_H

