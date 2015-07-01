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
	      Ep_(nMax,meshLambdaTotal),
	      Em_(nMax,meshLambdaTotal),
	      ep_(nMax,meshLambdaTotal),
	      em_(nMax,meshLambdaTotal),
	      kappa_(grounded_.kIndex().total(),0.0)
	{
		RealType xplus = (2*grounded.U()-mu)/T;
		RealType constant = T_*log(cosh(xplus));
		updateKappa(constant);
		result_ = updateResult();
	}

	RealType operator()() const
	{
		return result_;
	}

private:

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
			sum += s(lambda-sink)*(ep_(1,j)-em_(1,j));
		}

		return sum*lambdaMesh_.step();
	}

	RealType updateResult() const
	{
		RealType sum = 0.0;
		SizeType meshTotal = grounded_.kIndex().total();
		for (SizeType i = 0; i < meshTotal; ++i) {
			sum += grounded_.rho0(i)*G(kappa_[i]);
		}
		// FIXME: ADD MORE
		return sum*grounded_.kIndex().step();
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
	MeshType lambdaMesh_;
	MatrixRealType Ep_, Em_, ep_, em_;
	VectorRealType kappa_;
}; // class GrandPotential

} // namespace BetheAnsatz
#endif // GRANDPOTENTIAL_H

