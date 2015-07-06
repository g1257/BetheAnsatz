/*
Copyright (c) 2015, UT-Battelle, LLC

BetheAnsatz, Version 0.1

This file is part of BetheAnsatz.
BetheAnsatz is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
BetheAnsatz is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with BetheAnsatz. If not, see <http://www.gnu.org/licenses/>.
*/
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
	}; // class SomethingIntegrand

	class Energy0Integrand {

	public:

		typedef RealType_ RealType;

		Energy0Integrand(RealType U)
		    : oneOverPi_(1.0/acos(-1)),sigma0_(U)
		{}

		static RealType function(RealType k, void* vp)
		{
			SigmaZeroType* p = static_cast<SigmaZeroType*>(vp);
			const SigmaZeroType& sigma0 = *p;
			RealType U = sigma0.U();
			return cos(k)*(0.5 + U*cos(k)*integralOfSomething(k,sigma0))/M_PI;
		}

		SigmaZeroType& params() { return sigma0_; }

		RealType rho0(RealType k) const
		{
			return oneOverPi_*(0.5 + sigma0_.U()*cos(k)*integralOfSomething(k,sigma0_));
		}

		RealType kappa0Part(RealType k) const
		{
			return sigma0_.kappa0Part(k);
		}

		RealType sigma0(RealType lambda) const
		{
			return sigma0_(lambda);
		}

		const RealType& U() const
		{
			return sigma0_.U();
		}

	private:

		static RealType integralOfSomething(RealType k,
		                                    const SigmaZeroType& sigma0)
		{
			SomethingIntegrand somethingIntegrand(sigma0,k);
			PsimagLite::Integrator<SomethingIntegrand> integrator(somethingIntegrand);
			return integrator();
		}

		RealType oneOverPi_;
		SigmaZeroType sigma0_;
	}; // class Energy0Integrand

public:

	typedef RealType_ RealType;
	typedef Mesh<RealType> MeshType;

	Grounded(const ParametersType& params)
	    : e0_(0.0),
	      energy0Integrand_(params.U),
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
			rho0_[i] = energy0Integrand_.rho0(k);
		}

		for (SizeType i = 0; i < kIndex_.total(); ++i) {
			RealType k = kIndex_.x(i);
			kappa0_[i] = 2.0*cos(k)-4*energy0Integrand_.kappa0Part(k);
		}

		initSigmaZeroVector();

		e0_ = calculateE0();
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

	const RealType& e0() const { return e0_; }

	const MeshType& kIndex() const { return kIndex_; }

	const MeshType& lambdaIndex() const { return lambdaIndex_; }

	const RealType& U() const { return energy0Integrand_.U(); }

private:

	void initSigmaZeroVector()
	{
		SizeType meshLambdaTotal = lambdaIndex_.total();
		for (SizeType j = 0; j < meshLambdaTotal; ++j) {
			sigma0vector_[j] = energy0Integrand_.sigma0(lambdaIndex_.x(j));
		}
	}

	RealType calculateE0()
	{
		PsimagLite::Integrator<Energy0Integrand> integrator(energy0Integrand_);
		VectorRealType pts(2,0.0);
		pts[0] = -acos(-1);
		pts[1] = -pts[0];
		return integrator(pts)*(-2.0);
	}

	RealType e0_;
	Energy0Integrand energy0Integrand_;
	MeshType kIndex_;
	MeshType lambdaIndex_;
	VectorRealType rho0_;
	VectorRealType kappa0_;
	VectorRealType sigma0vector_;
}; // class Grounded
} // namespace BetheAnsatz

#endif // GROUNDED_H

