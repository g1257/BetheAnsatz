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

#ifndef GRANDPOTENTIAL_H
#define GRANDPOTENTIAL_H
#include "Grounded.h"
#include "Matrix.h"

namespace BetheAnsatz {

template<typename ParametersType>
class GrandPotential {

	typedef typename ParametersType::RealType RealType;

public:

	typedef Grounded<ParametersType> GroundedType;
	typedef typename GroundedType::MeshType MeshType;
	typedef typename MeshType::VectorRealType VectorRealType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;

	GrandPotential(const ParametersType& params,
	               const GroundedType& grounded,
	               RealType mu,
	               RealType T,
	               std::ostream& clog)
	    : grounded_(grounded),
	      mu_(mu),
	      T_(T),
	      result_(0.0),
	      Ep_(params.nMax,grounded_.lambdaIndex().total()),
	      Em_(params.nMax,Ep_.n_col()),
	      ep_(params.nMax,Ep_.n_col()),
	      em_(params.nMax,Ep_.n_col()),
	      kappa_(grounded_.kIndex().total(),0.0)
	{
		RealType xplus = (2*grounded.U()-mu_)/T_;
		RealType constant = T_*log(cosh(xplus));

		MatrixRealType bmatrix(2,params.nMax);
		for (SizeType i = 0; i < params.nMax; ++i) {
			RealType tmp = sinh(xplus)/sinh((i+2)*xplus);
			bmatrix(0,i) = tmp*tmp;
			tmp = 1.0/(i+2);
			bmatrix(1,i) = tmp*tmp;
		}

		RealType constant2 = -T_*log(bmatrix(1,1));

		clog<<"T="<<T_<<" mu="<<mu_<<"\n";
		for (SizeType it = 0; it < params.iterations; ++it) {
			iterate(it,params.nMax,constant,bmatrix);
			RealType result = updateResult(constant2);
			if (fabs(result-result_) < params.errorRelative && it > 0) break;
			result_ = result;
			clog<<it<<" "<<result_<<"\n";
		}

		clog<<"-------------\n";
	}

	RealType operator()() const
	{
		return result_;
	}

private:

	void iterate(SizeType,
	             SizeType nMax,
	             RealType constant,
	             const MatrixRealType& bmatrix)
	{
		updateKappa(constant);
		for (SizeType n = 0; n < nMax; ++n) {
			updateEpsilon(ep_,bmatrix,0,n);
			updateEpsilon(em_,bmatrix,1,n);
			updateE(Ep_,0,n);
			updateE(Em_,1,n);
		}

	}

	void updateEpsilon(MatrixRealType& epsilon,
	                   const MatrixRealType& bmatrix,
	                   SizeType plusOrMinus,
	                   SizeType n)
	{
		RealType b = bmatrix(plusOrMinus,n);
		RealType oneOverT = 1.0/T_;
		SizeType meshTotal = grounded_.lambdaIndex().total();
		for (SizeType j = 0; j < meshTotal; ++j) {
			RealType en = (plusOrMinus == 0) ? Ep_(n,j) : Em_(n,j);
			epsilon(n,j) = T_*log(b+(1-b)*exp(oneOverT*en));
		}
	}

	void updateE(MatrixRealType& E,
	             SizeType plusOrMinus,
	             SizeType n)
	{
		if (n == 0) return updateEfirst(E,plusOrMinus);

		SizeType meshTotal = grounded_.lambdaIndex().total();
		for (SizeType i = 0; i < meshTotal; ++i) {
			RealType lambda = grounded_.lambdaIndex().x(i);
			RealType sum = 0.0;
			for (SizeType j = 0; j < meshTotal; ++j) {
				RealType lambdaPrime = grounded_.lambdaIndex().x(j);
				RealType tmp1 = (plusOrMinus == 0) ? ep_(n-1,j) : em_(n-1,j);
				RealType tmp2 = 0;
				if (n + 1 < ep_.n_row())
					tmp2 = (plusOrMinus == 0) ? ep_(n+1,j) : em_(n+1,j);
				else
					tmp2 = 0.0;

				sum += s(lambda-lambdaPrime)*(tmp1 + tmp2);
			}

			E(n,i) = sum * grounded_.lambdaIndex().step();
		}
	}

	void updateEfirst(MatrixRealType& E,
	                  SizeType plusOrMinus)
	{
		SizeType meshTotal = grounded_.lambdaIndex().total();
		SizeType meshKtotal = grounded_.kIndex().total();
		RealType pm = (plusOrMinus == 0) ? 1.0 : -1.0;

		for (SizeType i = 0; i < meshTotal; ++i) {
			RealType lambda = grounded_.lambdaIndex().x(i);
			RealType sum = 0.0;
			for (SizeType j = 0; j < meshTotal; ++j) {
				RealType lambdaPrime = grounded_.lambdaIndex().x(j);
				RealType tmp = (plusOrMinus == 0) ? ep_(1,j) : em_(1,j);

				sum += s(lambda-lambdaPrime)*tmp;
			}

			sum *= grounded_.lambdaIndex().step();

			RealType sum2 = 0.0;
			for (SizeType j = 0; j < meshKtotal; ++j) {
				RealType k = grounded_.kIndex().x(j);
				sum2 += cos(k)*s(lambda-sin(k))*G(pm*kappa_[j]);
			}

			sum2 *= grounded_.kIndex().step();

			E(0,i) = sum - sum2;
		}
	}

	void updateKappa(RealType constant)
	{
		SizeType meshTotal = grounded_.kIndex().total();
		for (SizeType i = 0; i < meshTotal; ++i) {
			kappa_[i] = grounded_.kappa0(i) + constant + equation2point3(i);
		}
	}

	RealType equation2point3(SizeType ind) const
	{
		SizeType meshTotal = grounded_.lambdaIndex().total();
		RealType sum = 0.0;
		RealType k = grounded_.kIndex().x(ind);
		RealType sink = sin(k);
		for (SizeType j = 0; j < meshTotal; ++j) {
			RealType lambda = grounded_.lambdaIndex().x(j);
			sum += s(lambda-sink)*(ep_(0,j)-em_(0,j));
		}

		return sum*grounded_.lambdaIndex().step();
	}

	RealType updateResult(RealType constant) const
	{
		RealType sum = 0.0;
		SizeType meshTotal = grounded_.kIndex().total();
		for (SizeType i = 0; i < meshTotal; ++i) {
			sum += grounded_.rho0(i)*G(kappa_[i]);
		}

		sum *= grounded_.kIndex().step();

		RealType sum2 = 0.0;
		SizeType meshLambdaTotal = grounded_.lambdaIndex().total();
		for (SizeType j = 0; j < meshLambdaTotal; ++j) {
			RealType tmp = em_(0,j) + constant;
			sum2 += grounded_.sigma0(j)*tmp;
		}

		sum2 *= grounded_.lambdaIndex().step();

		return grounded_.e0() - (mu_ + sum + sum2);
	}

	RealType G(RealType y) const
	{
		return T_*log(1.0 + exp(y/T_));
	}

	RealType s(RealType lambda) const
	{
		RealType U = grounded_.U();
		RealType factor1 = 0.25/U;
		RealType factor2 = acos(-1)*0.5/U;
		return factor1/cosh(factor2*lambda);
	}

	const GroundedType& grounded_;
	RealType mu_;
	RealType T_;
	RealType result_;
	MatrixRealType Ep_, Em_, ep_, em_;
	VectorRealType kappa_;
}; // class GrandPotential

} // namespace BetheAnsatz
#endif // GRANDPOTENTIAL_H

