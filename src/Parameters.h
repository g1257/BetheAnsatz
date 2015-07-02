#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "Vector.h"

namespace BetheAnsatz {

template<typename RealType_, typename InputType>
struct Parameters {

	typedef RealType_ RealType;

	Parameters(InputType&)
	    : meshKtotal(1000),
	      meshLambdaTotal(10000),
	      nMax(50),
	      tt(0.0),
	      iterations(1),
	      infty(1e6),
	      tb(0.0),
	      ts(0.0),
	      U(0.0)
	{}

	SizeType meshKtotal;
	SizeType meshLambdaTotal;
	SizeType nMax;
	SizeType tt;
	SizeType iterations;
	RealType infty;
	RealType tb;
	RealType ts;
	RealType U;
}; // struct Parameters

} // namespace BetheAnsatz
#endif // PARAMETERS_H

