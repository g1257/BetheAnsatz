#ifndef GROUNDED_H
#define GROUNDED_H
#include "Mesh.h"
#include "SigmaZero.h"

namespace BetheAnsatz {

template<typename ParametersType>
class Grounded {

	typedef typename ParametersType::RealType RealType_;
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

	Grounded(const ParametersType& params)
	    : oneOverPi_(1.0/acos(-1)),
	      sigma0_(params.U),
	      kIndex_(params.meshKtotal,-acos(-1),2*acos(-1)/params.meshKtotal),
	      lambdaIndex_(params.meshLambdaTotal,
	                   -params.infty,
	                   2.0*params.infty/params.meshLambdaTotal),
	      rho0_(params.meshKtotal,0.0),
	      kappa0_(params.meshKtotal,0.0),
	      sigma0vector_(params.meshLambdaTotal,0.0)
	{
		for (SizeType i = 0; i < kIndex_.total(); ++i) {
			RealType k = kIndex_.x(i);
			rho0_[i] = oneOverPi_*(0.5 + params.U*cos(k)*integralOfSomething(k));
		}

		for (SizeType i = 0; i < kIndex_.total(); ++i) {
			RealType k = kIndex_.x(i);
			kappa0_[i] = 2.0*cos(k)-4*sigma0_.kappa0Part(k);
		}

		initSigmaZeroVector();
		std::cerr<<"Grounded::ctor: done\n";
	}

	RealType rho0(SizeType i) const
	{
		assert(i < rho0_.size());
		return rho0_[i];
	}

	RealType kappa0(SizeType i) const
	{
		assert(i < kappa0_.size());
		return kappa0_[i];
	}

	RealType sigma0(SizeType index) const
	{
		assert(index < sigma0vector_.size());
		return sigma0vector_[index];
	}

	const MeshType& kIndex() const { return kIndex_; }

	const MeshType& lambdaIndex() const { return lambdaIndex_; }

	const RealType& U() const { return sigma0_.U(); }

private:

	RealType integralOfSomething(RealType k) const
	{
		SomethingIntegrand somethingIntegrand(sigma0_,k);
		PsimagLite::Integrator<SomethingIntegrand> integrator(somethingIntegrand);
		return integrator();
	}

	void initSigmaZeroVector()
	{
		SizeType meshLambdaTotal = lambdaIndex_.total();
		for (SizeType j = 0; j < meshLambdaTotal; ++j) {
			sigma0vector_[j] = sigma0_(lambdaIndex_.x(j));
		}
	}

	RealType oneOverPi_;
	SigmaZeroType sigma0_;
	MeshType kIndex_;
	MeshType lambdaIndex_;
	VectorRealType rho0_;
	VectorRealType kappa0_;
	VectorRealType sigma0vector_;
}; // class Grounded
} // namespace BetheAnsatz

#endif // GROUNDED_H

