#ifndef SIGMAZERO_H
#define SIGMAZERO_H
#include "Integrator.h"

namespace BetheAnsatz {

template<typename RealType_>
class SigmaZero {

	class SigmaZeroIntegrand {

		struct Params {
			Params(RealType_ U, RealType_ lambda_)
			 :   lambda(lambda_), factor1(0.25/U), factor2(0.5*acos(-1)/U)
			{}

			RealType_ lambda;
			RealType_ factor1;
			RealType_ factor2;
		};

	public:

		typedef RealType_ RealType;

		SigmaZeroIntegrand(RealType U, RealType lambda)
		    : params_(U,lambda)
		{}

		static RealType function(RealType k, void *vp)
		{
			Params* p = static_cast<Params*>(vp);
			return oneOver2Pi_*s(p->lambda-sin(k),p->factor1,p->factor2);
		}

		Params& params() { return params_; }

	private:

		static RealType s(RealType lambda, RealType factor1, RealType factor2)
		{
			return factor1/cosh(factor2*lambda);
		}

		static RealType oneOver2Pi_;
		Params params_;

	}; // class SigmaZeroIntegrand

public:

	typedef RealType_ RealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;

	SigmaZero(RealType U) : U_(U)
	{}

	RealType operator()(RealType lambda) const
	{
		SigmaZeroIntegrand sigmaZeroIntegrand(U_,lambda);
		PsimagLite::Integrator<SigmaZeroIntegrand> integrator(sigmaZeroIntegrand);
		VectorRealType pts(2,0);
		pts[0] = -acos(-1);
		pts[1] = -pts[0];
		return integrator(pts);
	}

	const RealType& U() const { return U_; }

private:

	RealType U_;
}; // class SigmaZero

template<typename RealType>
RealType SigmaZero<RealType>::SigmaZeroIntegrand::oneOver2Pi_ = 0.5/acos(-1);
} // namespace BetheAnsatz
#endif // SIGMAZERO_H

