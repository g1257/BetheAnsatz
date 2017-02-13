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

	typedef typename LogEtaType::ParametersType ParametersType;
	typedef typename LogEtaType::RealType RealType;
	typedef typename LogEtaType::MatrixRealType MatrixRealType;
	typedef typename LogEtaType::MeshType MeshType;

	struct Auxiliary {
		Auxiliary(SizeType n_,
		          const LogEtaType& logEta_,
		          const MatrixRealType& rho_,
		          const MeshType& mesh_)
		    : n(n_),lastRho(0.0),logEta(logEta_),rho(rho_),mesh(mesh_)
		{}

		SizeType n;
		RealType lastRho;
		const LogEtaType& logEta;
		const MatrixRealType& rho;
		const MeshType& mesh;
	};

	typedef Auxiliary AuxiliaryType;
	typedef RealType (BracketableFunctionType)(SizeType, const AuxiliaryType&);
	typedef Convolution<BracketableFunctionType,AuxiliaryType,MeshType>
	ConvolutionType;

public:

	Rho(const ParametersType& params,
	    RealType temperature,
	    const LogEtaType& logEta)
	    : logEta_(logEta),mesh_(logEta_.mesh()),rho_(params.nMax,mesh_.total())
	{
		//clog<<"#Rho T="<<temperature<<" J="<<params.J<<"\n";
		RealType controlOld = 0;
		for (SizeType it = 0; it < params.iterations; ++it) {
			RealType control = calcRho1();
			for (SizeType n = 1; n < params.nMax; ++n) {
				control += calcRhoN(n);
			}

			control = (params.nMax - 1)*integralRho(params.nMax - 1);
			//clog<<it<<" "<<control<<"\n";
			if (fabs(1.0 - controlOld/control) < params.errorRelative) break;
			controlOld = control;
		}

		//clog<<"#Rho-------------\n";
	}

	const RealType& operator()(SizeType n, SizeType i) const
	{
		return rho_(n,i);
	}

	static RealType integralSz(SizeType n,
	                           const MatrixRealType& rho,
	                           const MeshType& mesh)
	{
		RealType sum = 0;
		for (SizeType i = 0; i < mesh.total(); ++i)
			sum += rho(n,i);

		return sum*mesh.step();
	}

	const MatrixRealType& matrix() const { return rho_; }

	template<typename T>
	friend	std::ostream& operator<<(std::ostream&, const Rho<T>&);

private:

	static RealType function1(SizeType i, const AuxiliaryType& a)
	{
		return a.rho(1,i)*exp(a.logEta(1,i));
	}

	static RealType function2(SizeType i, const AuxiliaryType& a)
	{
		return a.rho(0,i)*(1.0 + exp(a.logEta(0,i)));
	}

	static RealType function11(SizeType i, const AuxiliaryType& a)
	{
		SizeType n = a.n;
		assert(n > 0);
		RealType tmp1 = a.rho(n-1,i)*exp(a.logEta(n-1,i));
		RealType tmp2 = 0;
		if (n + 1 < a.rho.n_row())
			tmp2 = a.rho(n+1,i)*exp(a.logEta(n+1,i));
		else
			tmp2 = a.lastRho;
		return tmp1 + tmp2;
	}

	static RealType function12(SizeType i, const AuxiliaryType& a)
	{
		SizeType n = a.n;
		assert(n > 0);
		return a.rho(n,i)*(exp(a.logEta(n,i))+1.0);
	}

	RealType findLastRho() const
	{
		RealType sz = 0.5;
		SizeType nMax = rho_.n_row();
		for (SizeType n = 0; n < nMax; ++n) {
			sz -= (n+1)*integralSz(n,rho_,mesh_);
		}

		return sz/(nMax*rho_.n_col());
	}

	RealType calcRho1()
	{
		AuxiliaryType a1(0,logEta_,rho_,mesh_);
		ConvolutionType conv1(1,function1,a1);
		ConvolutionType conv2(2,function2,a1);

		RealType sum = 0;
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			RealType k = mesh_.x(i);
			RealType n1k = exp(logEta_(0,i));
			RealType factor = 1.0/(n1k + 1);
			RealType tmp = M_PI*(k*k + 1);
			rho_(0,i) = factor*(1.0/tmp + conv1(i) - conv2(i));
			sum += rho_(0,i);
		}

		return fabs(sum*mesh_.step());
	}

	RealType calcRhoN(SizeType n)
	{
		assert(n > 0);
		AuxiliaryType a11(n,logEta_,rho_,mesh_);
		a11.lastRho = 0.0; //findLastRho();
		ConvolutionType conv11(1,function11,a11);
		ConvolutionType conv12(2,function12,a11);

		RealType sum = 0;
		for (SizeType i = 0; i < mesh_.total(); ++i) {
			RealType nnk = exp(logEta_(n,i));
			RealType factor = 1.0/(nnk + 1);
			rho_(n,i) = factor*(conv11(i) - conv12(i));
			sum += rho_(0,i);
		}

		return fabs(sum*mesh_.step());
	}

	RealType integralRho(SizeType n) const
	{
		RealType sum = 0.0;
		for (SizeType i = 0; i < mesh_.total(); ++i)
			sum += rho_(n,i);

		return sum;
	}

	const LogEtaType& logEta_;
	const MeshType& mesh_;
	MatrixRealType rho_;
}; // class Rho

template<typename T>
std::ostream& operator<<(std::ostream& os, const Rho<T>& rho)
{
	os<<"LogRho\n";
	os<<rho.rho_;
	return os;
}

} // namespace BetheAnsatz

#endif // BETHE_RHO_H

