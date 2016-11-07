#ifndef TWOSPINONHEISENBERG_H
#define TWOSPINONHEISENBERG_H

#include <iostream>
#include "../../Engine/Mesh.h"
#include "Matrix.h"
#include "Integrator.h"
#include "SpecialFunctions.h"
#include <math.h>

namespace BetheAnsatz {

// Reference [1] M. Karbach et al., PRB 55, 12510 (1997)
template<typename RealType_>
class TwoSpinonHeisenberg {

	typedef typename PsimagLite::Vector<RealType_>::Type VectorRealType;
	typedef RealType_ RealType;

	static const RealType gamma_;

	class Fintegrand {

	public:

		enum {F1, F2};

		struct Params {
			Params(RealType_ rho_,
			       SizeType type_)
			    : rho(rho_),type(type_)
			{}

			RealType_ rho;
			SizeType type;
		};

		typedef RealType_ RealType;

		Fintegrand(RealType rho,SizeType type)
		    : params_(rho,type)
		{}

		Params& params() { return params_; }

		static RealType function(RealType x, void* ptr)
		{
			Params* paramsPtr = static_cast<Params*>(ptr);

			bool f1 = (paramsPtr->type == F1);

			RealType tmp = (f1) ? cosh(x) : coth(x);
			RealType num = cos(paramsPtr->rho*x);
			RealType den = x*tmp*tmp;
			return num/den;
		}

	private:

		static RealType coth(RealType x)
		{
			return cosh(x)/sinh(x);
		}

		Params params_;
	}; // class Fintegrand

public:

	TwoSpinonHeisenberg(const BetheAnsatz::Mesh<RealType>& kmesh,
	                    const BetheAnsatz::Mesh<RealType>& emesh)
	    : kmesh_(kmesh),
	      emesh_(emesh),
	      m_(kmesh.total(),emesh.total()),
	      i0Over2_(0.5*(gamma_ + f1Function(0.0) - f2Function(0.0)))
	{
		std::cerr<<"f1(0)= "<<f1Function(0.0)<<" f2(0)= "<<f2Function(0.0);
		std::cerr<<"I0="<<(2.0*i0Over2_)<<"\n";
		for (SizeType i = 0; i < kmesh.total(); ++i) {
			RealType k = kmesh.x(i);
			RealType wu = M_PI*sin(0.5*k);
			RealType wl = 0.5*M_PI*fabs(sin(k));
			for (SizeType j = 0; j < emesh.total(); ++j) {
				RealType w = emesh.x(j);
				// Eq. (1.7) of Ref [1]
				m_(i,j) = mFunction(wu,wl,w)*dFunction(wu,wl,w);
			}
		}
	}

	void toGnuplot(std::ostream& os, RealType cutoff)
	{
		for (SizeType j = 0; j < emesh_.total(); ++j) {
			RealType w = emesh_.x(j);
			for (SizeType i = 0; i < kmesh_.total(); ++i) {
				RealType k = kmesh_.x(i);
				RealType value = m_(i,j);
				if (cutoff > 0)
					value = (value > cutoff) ? cutoff : value;
				os<<k<<" "<<w<<" "<<value<<"\n";
			}
		}
	}

	RealType eMinusIfunction(RealType t) const
	{
		RealType tmp = 0.5*hFunction(t);
		return exp(tmp)*sqrt(t)*sinh(M_PI*t*0.25)*exp(i0Over2_);
	}

private:

	// Eq.(1.5) of Ref [1]
	RealType dFunction(RealType wu, RealType wl, RealType w)
	{
		if (w >= wu) return 0.0;
		if (w <= wl) return 0.0;
		RealType den = wu*wu-w*w;
		if (isTooSmall(den)) return 0.0;
		return 1.0/sqrt(den);
	}

	// Eq. (2.3) of Ref [1]
	RealType mFunction(RealType wu, RealType wl, RealType w)
	{
		if (w >= wu) return 0.0;
		if (w <= wl) return 0.0;
		RealType den = (w*w-wl*wl);
		if (isTooSmall(den)) return 0.0;
		RealType tmp2 = (wu*wu-wl*wl)/den;
		RealType t = 4.0*log(sqrt(tmp2) + sqrt(tmp2-1))/M_PI;
		return eMinusIfunction(t);
	}

	RealType hFunction(RealType t) const
	{
		return PsimagLite::Ci(t) + f1Function(t) - f2Function(t);
	}

	RealType f1Function(RealType t) const
	{
		Fintegrand f1Integrand(t,Fintegrand::F1);
		PsimagLite::Integrator<Fintegrand> integrator(f1Integrand);
		return integrator.toInfinity(1.0);
	}

	RealType f2Function(RealType t) const
	{
		Fintegrand f2Integrand(t,Fintegrand::F2);
		PsimagLite::Integrator<Fintegrand> integrator(f2Integrand);
		VectorRealType pts(2,0.0);
		pts[1] = 1.0;
		return integrator(pts);
	}

	bool isTooSmall(RealType x) const
	{
		return (fabs(x)<1e-10);
	}

	const BetheAnsatz::Mesh<RealType>& kmesh_;
	const BetheAnsatz::Mesh<RealType>& emesh_;
	PsimagLite::Matrix<RealType> m_;
	RealType i0Over2_;
};

template<typename RealType>
const RealType TwoSpinonHeisenberg<RealType>::gamma_ = 0.577215664901532860606512090082;
} // namespace BetheAnsatz

#endif // TWOSPINONHEISENBERG_H

