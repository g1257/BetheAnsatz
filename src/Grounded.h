#ifndef GROUNDED_H
#define GROUNDED_H
#include "Mesh.h"
#include "SigmaZero.h"

namespace BetheAnsatz {

template<typename RealType_>
class Grounded {

	typedef SigmaZero<RealType_> SigmaZeroType;
	typedef typename SigmaZeroType::VectorRealType VectorRealType;

	class SomethingIntegrand {

		struct Params {
			Params(const SigmaZeroType& sigma0_,
			       RealType_ k_)
			    : sigma0(sigma0_),k(k_),Usquared(sigma0.U()*sigma0.U())
			{}

			const SigmaZeroType& sigma0;
			RealType_ k;
			RealType_ Usquared;
		};

	public:

		typedef RealType_ RealType;

		SomethingIntegrand(const SigmaZeroType& sigma0,
		                   RealType k)
		    : params_(sigma0,k)
		{}

		static RealType function(RealType lambda, void* vp)
		{
			Params* p = static_cast<Params*>(vp);
			RealType tmp = (lambda - sin(p->k));
			RealType den = tmp*tmp + p->Usquared;
			return p->sigma0(lambda)/den;
		}


		Params& params() { return params_; }

	private:

		Params params_;
	};

public:

	typedef RealType_ RealType;
	typedef Mesh<RealType> MeshType;

	Grounded(RealType U, SizeType meshTotal)
	    : oneOverPi_(1.0/acos(-1)),
	      sigma0_(U),
	      rho0_(meshTotal,-acos(-1),2*acos(-1)/meshTotal)
	{

		for (SizeType i = 0; i < rho0_.total(); ++i) {
			RealType k = rho0_.x(i);
			//IntegrableFunctionType f(sigma0_,k);
			rho0_[i] = oneOverPi_*(0.5 + U*cos(k)*integralOfSomething(k));
		}
	}

	RealType rho0(SizeType i) const
	{
		return rho0_[i];
	}

private:

	RealType integralOfSomething(RealType k) const
	{
		SomethingIntegrand somethingIntegrand(sigma0_,k);
		PsimagLite::Integrator<SomethingIntegrand> integrator(somethingIntegrand);
		return integrator();
	}

	RealType oneOverPi_;
	SigmaZeroType sigma0_;
	MeshType rho0_;
}; // class Grounded
} // namespace BetheAnsatz

#endif // GROUNDED_H

