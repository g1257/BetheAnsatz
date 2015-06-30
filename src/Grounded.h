#ifndef GROUNDED_H
#define GROUNDED_H
#include "Mesh.h"
#include "SigmaZero.h"

namespace BetheAnsatz {

template<typename RealType>
class Grounded {

public:

	typedef Mesh<RealType> MeshType;
	typedef SigmaZero<RealType> SigmaZeroType;

	Grounded(RealType U, SizeType meshTotal)
	    : sigma0_(U),rho0_(meshTotal,-acos(-1),2*acos(-1)/meshTotal)
	{

		for (SizeType i = 0; i < rho0_.total(); ++i) {
			RealType k = rho0_.x(i);
			//IntegrableFunctionType f(sigma0_,k);
			rho0_[i] = 0.5 + U*cos(k); //*integralOf(f,-infty,infty);
		}
	}

private:

	SigmaZeroType sigma0_;
	MeshType rho0_;
}; // class Grounded
} // namespace BetheAnsatz

#endif // GROUNDED_H

