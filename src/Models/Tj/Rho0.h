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
#ifndef BETHE_RHO0_H
#define BETHE_RHO0_H
#include "../../Engine/Mesh.h"
#include "Integrator.h"
#include "RootFindingBisection.h"

// See PRB 44, 130 (1991)
namespace BetheAnsatz {

template<typename ParametersType_>
class Rho0 {

	typedef ParametersType_ ParametersType;
	typedef typename ParametersType::RealType RealType;
	typedef Mesh<RealType> MeshType;
	typedef typename MeshType::VectorRealType VectorRealType;

	class ShibaIntegrand {

	public:

		typedef typename ParametersType::RealType RealType;

		static RealType function(RealType lambda, void* vp)
		{
			RealType* p = static_cast<RealType*>(vp);
			return cos((*p)*lambda)/(1.0+exp(fabs(lambda)));
		}

		RealType& params() { return xOver2_; }

		ShibaIntegrand(RealType x) : xOver2_(0.5*x){}

	private:

		RealType xOver2_;
	}; // class ShibaIntegrand

	class EqThreePointFiveB {

	public:

		typedef typename ParametersType::RealType RealType;

		EqThreePointFiveB(const MeshType& mesh,
		                  const VectorRealType& rho0,
		                  RealType density)
		    : mesh_(mesh),rho0_(rho0),densityOver2_(0.5*density)
		{}

		RealType operator()(RealType q0) const
		{
			SizeType q0Index = static_cast<SizeType>((q0-mesh_.x(0))/mesh_.step());
			RealType sum = -densityOver2_;
			assert(q0Index <= rho0_.size());
			for (SizeType i = 0; i < q0Index; ++i) sum += rho0_[i];

			RealType start = mesh_.total() - q0Index;
			for (SizeType i = start; i < mesh_.total(); ++i) sum += rho0_[i];

			return sum;
		}

	private:

		const MeshType& mesh_;
		const VectorRealType& rho0_;
		RealType densityOver2_;
	}; // class EqThreePointFiveB

public:

	Rho0(const ParametersType& params,
	     RealType density,
	     std::ostream&)
	    : mesh_(params.meshLambdaTotal,-params.infty,
	            2.0*params.infty/params.meshLambdaTotal),
	      density_(density),
	      rho0_(mesh_.total())
	{
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			RealType nu = mesh_.x(i);
			rho0_[i] = 1.0/(M_PI*(1.0 + nu*nu));
		}

		VectorRealType rho0(rho0_.size());
		for (SizeType it = 0; it < params.iterations; ++it) {
			SizeType q0Index = computeQ0Index();
			computeNewRho0(rho0,q0Index);
			rho0_ = rho0;
		}
	}

	const MeshType& mesh() const { return mesh_; }

	const RealType& operator()(SizeType i) const
	{
		assert(i < rho0_.size());
		return rho0_[i];
	}

	const RealType energy() const
	{
		SizeType zeroIndex = static_cast<SizeType>(-1.0*mesh_.x(0)/mesh_.step());
		return 2.0*(1.0-density_-M_PI*rho0_[zeroIndex]);
	}

	template<typename T>
	friend std::ostream& operator<<(std::ostream& os, const Rho0<T>& r);

private:

	void computeNewRho0(VectorRealType& rho0, SizeType q0Index) const
	{
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			RealType nu = mesh_.x(i);
			rho0[i] = 2*shibaR(2*nu) + integralOfSomething(nu,q0Index);
		}
	}

	SizeType computeQ0Index() const
	{
		SizeType end = mesh_.total() - 1;
		EqThreePointFiveB eqThreePointFive(mesh_,rho0_,density_);
		PsimagLite::RootFindingBisection<EqThreePointFiveB> rfb(eqThreePointFive,
		                                                        mesh_.x(0),
		                                                        mesh_.x(end));
		RealType q0 = rfb();
		return static_cast<SizeType>((q0-mesh_.x(0))/mesh_.step());
	}

	RealType shibaR(RealType x) const
	{
		ShibaIntegrand shibaIntegrand(x);
		PsimagLite::Integrator<ShibaIntegrand> integrator(shibaIntegrand);
		return integrator()/(4.0*M_PI);
	}

	RealType integralOfSomething(RealType nu,SizeType q0Index) const
	{
		SizeType end = mesh_.total() - q0Index;
		RealType sum = 0.0;
		for (SizeType i = q0Index; i < end; ++i) {
			RealType nuPrime = mesh_.x(i);
			sum += shibaR(2*(nu-nuPrime))*rho0_[i];
		}

		return 2.0*sum*mesh_.step();
	}

	MeshType mesh_;
	RealType density_;
	VectorRealType rho0_;
}; // class Rho0

template<typename T>
std::ostream& operator<<(std::ostream& os, const Rho0<T>& r)
{
	os<<"Rho0\n";
	os<<r.mesh();
	os<<r.rho_;
	return os;
}

} // namespace BetheAnsatz

#endif // BETHE_RHO0_H

