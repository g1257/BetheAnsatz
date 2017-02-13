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
#ifndef BETHE_LOGETA_H
#define BETHE_LOGETA_H
#include "../../Engine/Mesh.h"
#include "Convolution.h"

namespace BetheAnsatz {

template<typename ParametersType_>
class LogEta {

public:

	typedef ParametersType_ ParametersType;
	typedef typename ParametersType::RealType RealType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;
	typedef Mesh<RealType> MeshType;
	typedef typename MeshType::VectorRealType VectorRealType;

	struct Auxiliary {
		Auxiliary(SizeType n_,
		          const MatrixRealType& logEta_,
		          const MeshType& mesh_)
		    : n(n_),logEta(logEta_),mesh(mesh_)
		{}

		SizeType n;
		const MatrixRealType& logEta;
		const MeshType& mesh;
	};

	typedef Auxiliary AuxiliaryType;
	typedef RealType (BracketableFunctionType)(SizeType, const AuxiliaryType&);
	typedef Convolution<BracketableFunctionType,AuxiliaryType,MeshType>
	ConvolutionType;

	LogEta(const ParametersType& params,
	           RealType temperature)
	    : mesh_(params.meshLambdaTotal,-params.infty,2.0*params.infty/params.meshLambdaTotal),
	      logEta_(params.nMax,mesh_.total()),
	      minusTwoJOverT_(-2.0*params.J/temperature)
	{
		//clog<<"#LogEta T="<<temperature<<" J="<<params.J<<"\n";
		RealType controlOld = 0;
		for (SizeType it = 0; it < params.iterations; ++it) {
			RealType control = calcEta1();
			for (SizeType n = 1; n < params.nMax; ++n) {
				control += calcEtaN(n);
			}

			control = maxLogEta(params.nMax - 1);
			//clog<<it<<" "<<control<<"\n";
			if (fabs(1.0 - controlOld/control) < params.errorRelative) break;
			controlOld = control;
		}

		//clog<<"RemainderCutoff="<<remainderCutoff()<<"\n";
		//clog<<"#LogEta-------------\n";
	}

	const MeshType& mesh() const { return mesh_; }

	const RealType& operator()(SizeType n, SizeType i) const
	{
		return logEta_(n,i);
	}

	template<typename T>
	friend std::ostream& operator<<(std::ostream&, const LogEta<T>&);

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

	static RealType trappedLogModif(RealType x)
	{
		if (x>25) return x;
		return log(1+exp(x));
	}

	RealType calcEta1()
	{
		AuxiliaryType a1(0,logEta_,mesh_);
		ConvolutionType conv1(1,function1,a1);
		ConvolutionType conv2(2,function2,a1);

		RealType sum = 0;
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			RealType k = mesh_.x(i);
			logEta_(0,i) = minusTwoJOverT_/(1.0 + k*k) + conv1(i) - conv2(i);
			if (fabs(logEta_(0,i)) > sum)
				sum = fabs(logEta_(0,i));
		}

		return sum;
	}

	RealType calcEtaN(SizeType n)
	{
		assert(n > 0);
		AuxiliaryType a11(n,logEta_,mesh_);
		ConvolutionType conv11(1,function11,a11);
		ConvolutionType conv12(2,function12,a11);

		RealType sum = 0;
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			logEta_(n,i) = conv11(i) - conv12(i);
			if (fabs(logEta_(n,i)) > sum)
				sum = fabs(logEta_(n,i));
		}

		return sum;
	}

	RealType remainderCutoff() const
	{
		assert(logEta_.n_row() > 0);
		SizeType last = logEta_.n_row() - 1;
		RealType sum = 0;
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			sum += fabs(logEta_(last,i));
		}

		return sum;
	}

	RealType maxLogEta(SizeType n) const
	{
		RealType max = 0.0;
		for (SizeType i = 0; i < mesh_.total(); ++i)
			if (fabs(logEta_(n,i)) > max) max = logEta_(n,i);

		return max;
	}

	MeshType mesh_;
	MatrixRealType logEta_;
	RealType minusTwoJOverT_;
}; // class LogEta

template<typename T>
std::ostream& operator<<(std::ostream& os, const LogEta<T>& logEta)
{
	os<<"LogEta\n";
	os<<logEta.mesh();
	os<<logEta.logEta_;
	os<<"minusTwoJOverT="<<logEta.minusTwoJOverT_<<"\n";
	return os;
}

} // namespace BetheAnsatz

#endif // BETHE_LOGETA_H

