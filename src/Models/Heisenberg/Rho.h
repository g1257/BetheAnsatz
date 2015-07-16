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

#ifndef BETHE_RHO_H
#define BETHE_RHO_H

namespace BetheAnsatz {

template<typename LogEtaType>
class Rho {

	typedef typename LogEtaType::RealType RealType;
	typedef typename LogEtaType::MatrixRealType MatrixRealType;
	typedef typename LogEtaType::MeshType MeshType;

	struct Auxiliary {
		Auxiliary(SizeType n_,
		          const LogEtaType& logEta_,
		          const MatrixRealType& rho_)
		    : n(n_),logEta(logEta_),rho(rho_)
		{}

		SizeType n;
		const LogEtaType& logEta;
		const MatrixRealType& rho;
	};

	typedef Auxiliary AuxiliaryType;
	typedef RealType (BracketableFunctionType)(SizeType, const AuxiliaryType&);
	typedef Convolution<BracketableFunctionType,AuxiliaryType,MeshType>
	ConvolutionType;

public:

	Rho(const LogEtaType& logEta) : logEta_(logEta),mesh_(logEta_.mesh())
	{}

private:

	static RealType function1(SizeType i, const AuxiliaryType& a)
	{
		return 0.0;
	}

	static RealType function2(SizeType i, const AuxiliaryType& a)
	{
		return a.logEta(0,i);
	}

	RealType calRho1()
	{
		AuxiliaryType a1(0,logEta_,rho_);
		ConvolutionType conv1(1,function1,a1,mesh_);
		ConvolutionType conv2(2,function2,a1,mesh_);

		RealType sum = 0;
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			RealType k = mesh_.x(i);
			RealType n1k = exp(logEta_(i));
			RealType tmp = M_PI*(k*k + 1)*(n1k + 1);
			rho_(0,i) = 1.0/tmp + conv1(i) - conv2(i);
			sum += log(rho_(0,i));
		}

		return fabs(sum*mesh_.step());
	}

	const LogEtaType& logEta_;
	const MeshType& mesh_;
	MatrixRealType rho_;
}; // class Rho

} // namespace BetheAnsatz

#endif // BETHE_RHO_H

