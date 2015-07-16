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
#ifndef BETHE_HEISENBERG_H
#define BETHE_HEISENBERG_H
#include "Mesh.h"

namespace BetheAnsatz {

template<typename ParametersType>
class Heisenberg {

	typedef typename ParametersType::RealType RealType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef typename PsimagLite::Vector<RealType>::Type VectorRealType;
	typedef Mesh<RealType> MeshType;

	struct Auxiliary {
		Auxiliary(SizeType n_,const MatrixRealType& logEta_)
		    : n(n_),logEta(logEta_)
		{}

		SizeType n;
		const MatrixRealType& logEta;
	};

	typedef Auxiliary AuxiliaryType;
	typedef RealType (BracketableFunctionType)(SizeType, const AuxiliaryType&);

	class Convolution {

	public:

		Convolution(SizeType n,
		            BracketableFunctionType& f,
		            const AuxiliaryType& a,
		            const MeshType& mesh)
		    : result_(mesh.total()),mesh_(mesh)
		{

			for (SizeType i = 0; i < mesh_.total(); ++i) {
				RealType k = mesh_.x(i);
				RealType sum = 0.0;
				for (SizeType j = 0; j < mesh_.total(); ++j) {
					RealType kdiff = k - mesh_.x(j);
					RealType tmp = n*n + kdiff*kdiff;
					sum += n*f(j,a)/tmp;
				}

				assert(i < result_.size());
				result_[i] = sum*mesh_.step()/M_PI;
			}
		}

		const RealType& operator()(SizeType i) const
		{
			assert(i < result_.size());
			return result_[i];
		}

	private:

		VectorRealType result_;
		const MeshType& mesh_;
	}; // class Convolution

public:

	Heisenberg(const ParametersType& params,
	           RealType temperature,
	           std::ostream& clog)
	    : mesh_(2*params.infty,-params.infty,2.0*params.infty/params.meshLambdaTotal),
	      logEta_(params.nMax,mesh_.total()),
	      minusTwoJOverT_(-2.0*params.J/temperature)
	{
		clog<<"T="<<temperature<<" J="<<params.J<<"\n";
		RealType controlOld = 0;
		for (SizeType it = 0; it < params.iterations; ++it) {
			RealType control = calcEta1();
			for (SizeType n = 1; n < params.nMax; ++n) {
				control += calcEtaN(n);
			}

			clog<<it<<" "<<control<<"\n";
			if (fabs(1.0 - controlOld/control) < 1e-6) break;
			controlOld = control;
		}

		clog<<"-------------\n";
	}

private:

	static RealType function1(SizeType i, const AuxiliaryType& a)
	{
		return trappedLogModif(a.logEta(1,i));
	}

	static RealType function2(SizeType i, const AuxiliaryType& a)
	{
		return a.logEta(0,i);
	}

	static RealType function11(SizeType i, const AuxiliaryType& a)
	{
		assert(a.n > 0);
		RealType tmp1 = trappedLogModif(a.logEta(a.n-1,i));
		if (a.n+1 < a.logEta.n_row())
			tmp1 += trappedLogModif(a.logEta(a.n+1,i));
		return tmp1;
	}

	static RealType function12(SizeType i, const AuxiliaryType& a)
	{
		assert(a.n > 0);
		return a.logEta(a.n,i);
	}

	RealType calcEta1()
	{
		AuxiliaryType a1(0,logEta_);
		Convolution conv1(1,function1,a1,mesh_);
		Convolution conv2(2,function2,a1,mesh_);

		RealType sum = 0;
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			RealType k = mesh_.x(i);
			logEta_(0,i) = minusTwoJOverT_/(1.0 + k*k) + conv1(i) - conv2(i);
			sum += logEta_(0,i);
		}

		return fabs(sum*mesh_.step());
	}

	RealType calcEtaN(SizeType n)
	{
		assert(n > 0);
		AuxiliaryType a11(n,logEta_);
		Convolution conv11(1,function11,a11,mesh_);
		Convolution conv12(2,function12,a11,mesh_);

		RealType sum = 0;
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			logEta_(n,i) = conv11(i) - conv12(i);
			sum += logEta_(n,i);
		}

		return fabs(sum*mesh_.step());
	}

	static RealType trappedLogModif(RealType x)
	{
		if (x>25) return x;
		return log(1+exp(x));
	}

	MeshType mesh_;
	MatrixRealType logEta_;
	RealType minusTwoJOverT_;
}; // class Heisenberg
} // namespace BetheAnsatz

#endif // BETHE_HEISENBERG_H

