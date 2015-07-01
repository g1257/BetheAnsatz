#ifndef SIGMAZERO_H
#define SIGMAZERO_H
#include "Integrator.h"

namespace BetheAnsatz {

template<typename RealType_>
class SigmaZero {

public:

	enum TypeEnum {
		TYPE_SIGMAZERO, TYPE_KAPPA0
	};

private:

	class SigmaZeroIntegrand {

		struct Params {
			Params(RealType_ U_, RealType_ lambda_, TypeEnum type_)
			 : U(U_),
			   lambda(lambda_),
			   factor1(0.25/U),
			   factor2(0.5*acos(-1)/U),
			   type(type_)
			{}

			RealType_ U;
			RealType_ lambda;
			RealType_ factor1;
			RealType_ factor2;
			TypeEnum type;
		};

	public:

		typedef RealType_ RealType;

		SigmaZeroIntegrand(RealType U, RealType lambda, TypeEnum type)
		    : params_(U,lambda, type)
		{}

		static RealType function(RealType k, void *vp)
		{
			Params* p = static_cast<Params*>(vp);
			return (p->type == TYPE_SIGMAZERO) ?
			            oneOver2Pi_*s(p->lambda-sin(k),p->factor1,p->factor2) :
			            kappa0Integrand(k,p->lambda,p->U,p->factor1,p->factor2);
		}

		Params& params() { return params_; }

	private:

		static RealType s(RealType lambda, RealType factor1, RealType factor2)
		{
			return factor1/cosh(factor2*lambda);
		}

		static RealType kappa0Integrand(RealType lambda,
		                                RealType k,
		                                RealType Uvalue,
		                                RealType factor1,
		                                RealType factor2)
		{

			RealType a = 1 + Uvalue*Uvalue - lambda*lambda;
			RealType b = 2*Uvalue*lambda;
			RealType tmp = a + sqrt(a*a + b*b);
			return s(lambda - sin(k), factor1, factor2) * sqrt(0.5*tmp);
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
		SigmaZeroIntegrand sigmaZeroIntegrand(U_, lambda, TYPE_SIGMAZERO);
		PsimagLite::Integrator<SigmaZeroIntegrand> integrator(sigmaZeroIntegrand);
		VectorRealType pts(2,0);
		pts[0] = -acos(-1);
		pts[1] = -pts[0];
		return integrator(pts);
	}

	RealType kappa0Part(RealType k)
	{
		SigmaZeroIntegrand kappa0Integrand(U_, k, TYPE_KAPPA0);
		PsimagLite::Integrator<SigmaZeroIntegrand> integrator(kappa0Integrand);
		return integrator();
	}

	const RealType& U() const { return U_; }

private:

	RealType U_;
}; // class SigmaZero

template<typename RealType>
RealType SigmaZero<RealType>::SigmaZeroIntegrand::oneOver2Pi_ = 0.5/acos(-1);
} // namespace BetheAnsatz
#endif // SIGMAZERO_H

