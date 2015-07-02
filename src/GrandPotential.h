#ifndef GRANDPOTENTIAL_H
#define GRANDPOTENTIAL_H
#include "Grounded.h"
#include "Matrix.h"

namespace BetheAnsatz {

template<typename RealType>
class GrandPotential {

public:

	typedef Grounded<RealType> GroundedType;
	typedef typename GroundedType::MeshType MeshType;
	typedef typename MeshType::VectorRealType VectorRealType;
	typedef PsimagLite::Matrix<RealType> MatrixRealType;

	GrandPotential(const GroundedType& grounded,
	               RealType mu,
	               RealType T,
	               RealType infty,
	               SizeType meshLambdaTotal,
	               SizeType nMax)
	    : grounded_(grounded),
	      mu_(mu),
	      T_(T),
	      result_(0.0),
	      lambdaMesh_(-infty,meshLambdaTotal,2.0*infty/meshLambdaTotal),
	      sigma0_(lambdaMesh_.total(),0.0),
	      Ep_(nMax,meshLambdaTotal),
	      Em_(nMax,meshLambdaTotal),
	      ep_(nMax,meshLambdaTotal),
	      em_(nMax,meshLambdaTotal),
	      kappa_(grounded_.kIndex().total(),0.0)
	{
		RealType xplus = (2*grounded.U()-mu)/T;
		RealType constant = T_*log(cosh(xplus));
		initSigmaZero();
		SizeType iterations = 1;

		MatrixRealType bmatrix(2,nMax);
		for (SizeType i = 0; i < nMax; ++i) {
			RealType tmp = sinh(xplus)/sinh((i+2)*xplus);
			bmatrix(0,i) = tmp*tmp;
			tmp = 1.0/(i+2);
			bmatrix(1,i) = tmp*tmp;
		}

		for (SizeType it = 0; it < iterations; ++it) {
			iterate(it,nMax,constant,bmatrix);
		}

		RealType constant2 = -T_*log(bmatrix(1,1));
		result_ = updateResult(constant2);
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
		}

	}

	void updateEpsilon(MatrixRealType& epsilon,
	                   const MatrixRealType& bmatrix,
	                   SizeType plusOrMinus,
	                   SizeType n)
	{
		RealType b = bmatrix(plusOrMinus,n);
		RealType oneOverT = 1.0/T_;
		SizeType meshTotal = lambdaMesh_.total();
		for (SizeType j = 0; j < meshTotal; ++j) {
			RealType en = (plusOrMinus == 0) ? Ep_(n,j) : Em_(n,j);
			epsilon(n,j) = T_*log(b+(1-b)*exp(oneOverT*en));
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
		SizeType meshTotal = lambdaMesh_.total();
		RealType sum = 0.0;
		RealType k = grounded_.kIndex().x(ind);
		RealType sink = sin(k);
		for (SizeType j = 0; j < meshTotal; ++j) {
			RealType lambda = lambdaMesh_.x(j);
			sum += s(lambda-sink)*(ep_(0,j)-em_(0,j));
		}

		return sum*lambdaMesh_.step();
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
		SizeType meshLambdaTotal = lambdaMesh_.total();
		for (SizeType j = 0; j < meshLambdaTotal; ++j) {
			RealType tmp = em_(0,j) + constant;
			sum2 += sigma0_[j]*tmp;
		}

		sum2 *= lambdaMesh_.step();
		// FIXME: ADD e0 here
		return -(mu_ + sum + sum2);
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

	void initSigmaZero()
	{
		SizeType meshLambdaTotal = lambdaMesh_.total();
		for (SizeType j = 0; j < meshLambdaTotal; ++j) {
			sigma0_[j] = grounded_.sigma0(lambdaMesh_.x(j));
		}
	}

	const GroundedType& grounded_;
	RealType mu_;
    RealType T_;
	RealType result_;
	MeshType lambdaMesh_;
	VectorRealType sigma0_;
	MatrixRealType Ep_, Em_, ep_, em_;
	VectorRealType kappa_;
}; // class GrandPotential

} // namespace BetheAnsatz
#endif // GRANDPOTENTIAL_H

