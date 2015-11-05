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

namespace BetheAnsatz {

template<typename ParametersType_>
class Rho0 {

public:

	typedef ParametersType_ ParametersType;
	typedef typename ParametersType::RealType RealType;
	typedef Mesh<RealType> MeshType;
	typedef typename MeshType::VectorRealType VectorRealType;

	Rho0(const ParametersType& params,
	     std::ostream& clog)
	    : mesh_(params.meshLambdaTotal,-params.infty,
	            2.0*params.infty/params.meshLambdaTotal),
	      rho0_(mesh_.total())
	{
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			RealType nu = mesh_(i);
			rho0_[i] = 1.0/(M_PI*(1.0 + nu*nu));
		}

		VectorRealType rho0(rho0_.size());
		for (SizeType it = 0; it < params.iterations; ++it) {
			RealType q0 = computeQ0();
			computeNewRho0(rho0,q0);
			rho0_ = rho0;
		}
	}

	const MeshType& mesh() const { return mesh_; }

	const RealType& operator()(SizeType i) const
	{
		assert(i < rho0_.size());
		return rho0_[i];
	}

	template<typename T>
	friend std::ostream& operator<<(std::ostream& os, const Rho0<T>& r);

private:

	void computeNewRho0(VectorRealType& rho0, RealType q0) const
	{
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			RealType nu = mesh_(i);
			rho0[i] = 2*shibaR(2*nu) + integralOfSomething(nu,q0);

		}
	}

	RealType computeQ0() const
	{
		throw PsimagLite::RuntimeError("computeQ0 unimplemented");
	}

	RealType shibaR(RealType) const
	{
		throw PsimagLite::RuntimeError("shibaR unimplemented");
	}

	RealType integralOfSomething(RealType nu,RealType q0) const
	{
		throw PsimagLite::RuntimeError("integralOfSomething unimplemented");
	}

	MeshType mesh_;
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

