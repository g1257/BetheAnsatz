#ifndef TWOSPINONHEISENBERG_H
#define TWOSPINONHEISENBERG_H

#include <iostream>
#include "../../Engine/Mesh.h"
#include "Matrix.h"
#include "Integrator.h"

// Reference [1] M. Karbach et al., PRB 55, 12510 (1997)
template<typename RealType_>
class TwoSpinonHeisenberg {

	typedef typename PsimagLite::Vector<RealType_>::Type VectorRealType;
	typedef RealType_ RealType;

	static const RealType gamma_ = 0.571572664901532860606512090082;

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

	class CosIntegral {

	public:

		CosIntegral(SizeType pre = 100, RealType tol = 1e-6)
		    : pre_(pre),toleranceSeries_(tol)
		{
			RealType prev = 0.0;
			for (SizeType i = 0; i < pre_.size(); ++i) {
				pre_[i] = preCompute(prev,i+1);
			}
		}

		RealType operator()(RealType x) const
		{
			return (x > 20) ? 0.0 : gamma_ + log(x) + 0.5*series(x);
		}

	private:

		RealType series(RealType x) const
		{
			RealType sum = 0.0;
			RealType prev = 0.0;
			SizeType total = pre_.size();
			for (SizeType i = 0; i < total; ++i) {
				SizeType j = i + 1;
				prev = sum;
				RealType tmp = pre_[i]*pow(x,2*j);
				RealType sign = (i&1) ? 1 : -1;
				sum += tmp*sign;
			}

			if (fabs(sum-prev) > toleranceSeries_) {
				PsimagLite::String str("CosIntegral: series ");
				str += "tolerance not achieved\n";
				throw PsimagLite::RuntimeError(str);
			}

			return sum;
		}

		RealType preCompute(RealType& prev,SizeType i) const
		{
			return 1.0/(i*factorialTwo(prev,i));
		}

		RealType factorialTwo(RealType& prev, SizeType x) const
		{
			if (x == 1) {
				prev = 2;
			}

			SizeType start = 2*x-1;
			SizeType end = 2*x + 1;
			for (SizeType i = start; i < end; ++i) prev *= i;

			return prev;
		}

		VectorRealType pre_;
		RealType toleranceSeries_;
	}; // class CosIntegral

public:

	TwoSpinonHeisenberg(const BetheAnsatz::Mesh<RealType>& kmesh,
	                    const BetheAnsatz::Mesh<RealType>& emesh)
	    : kmesh_(kmesh),
	      emesh_(emesh),
	      m_(kmesh.total(),emesh.total()),
	      cosIntegral_()
	{
		RealType i0Over2 = 0.5*(gamma_ + f1Function(0.0) - f2Function(0.0));
		for (SizeType i = 0; i < kmesh.total(); ++i) {
			RealType k = kmesh.x(i);
			RealType wu = M_PI*sin(0.5*k);
			RealType wl = 0.5*M_PI*fabs(sin(k));
			for (SizeType j = 0; j < emesh.total(); ++j) {
				RealType w = emesh.x(j);
				// Eq. (1.7) of Ref [1]
				m_(i,j) = mFunction(wu,wl,w)*dFunction(wu,wl,w)*exp(i0Over2);
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
		RealType t = 4.0*(sqrt(tmp2) + sqrt(tmp2+1))/M_PI;
		RealType tmp = 0.5*hFunction(t);
		return 0.5*exp(tmp)*sqrt(t)*sinh(M_PI*t*0.25);
	}

	RealType hFunction(RealType t)
	{
		return cosIntegral_(t) + f1Function(t) - f2Function(t);
	}

	RealType f1Function(RealType t)
	{
		Fintegrand f1Integrand(t,Fintegrand::F1);
		PsimagLite::Integrator<Fintegrand> integrator(f1Integrand);
		return integrator.toInfinity(1.0);
	}

	RealType f2Function(RealType t)
	{
		return 0.0;
		Fintegrand f2Integrand(t,Fintegrand::F2);
		PsimagLite::Integrator<Fintegrand> integrator(f2Integrand);
		VectorRealType pts(2,0.0);
		pts[1] = 1.0;
		return integrator(pts);
	}

	bool isTooSmall(RealType x) const
	{
		return (fabs(x)<1e-8);
	}

	const BetheAnsatz::Mesh<RealType>& kmesh_;
	const BetheAnsatz::Mesh<RealType>& emesh_;
	PsimagLite::Matrix<RealType> m_;
	CosIntegral cosIntegral_;
};

#endif // TWOSPINONHEISENBERG_H

