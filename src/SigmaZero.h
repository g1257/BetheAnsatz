#ifndef SIGMAZERO_H
#define SIGMAZERO_H

namespace BetheAnsatz {

template<typename RealType>
class SigmaZero {

	class SigmaZeroInternal {

	public:

		SigmaZeroInternal(RealType U, RealType delta)
		    : delta_(delta), factor1_(0.25/U), factor2_(0.5*acos(-1)/U)
		{}

		RealType operator()(RealType k) const
		{
			return s(delta_-sin(k));
		}

	private:

		RealType s(RealType delta) const
		{
			return factor1_/cosh(factor2_*delta);
		}

		RealType delta_;
		RealType factor1_;
		RealType factor2_;

	}; // class SigmaZeroInternal

public:

	SigmaZero(RealType U) : U_(U)
	{}

	RealType operator()(RealType delta) const
	{
		SigmaZeroInternal sigmaZeroInternal(U_,delta);
		return 0.0;
	}

private:

	RealType U_;
}; // class SigmaZero

} // namespace BetheAnsatz
#endif // SIGMAZERO_H

