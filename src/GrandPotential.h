#ifndef GRANDPOTENTIAL_H
#define GRANDPOTENTIAL_H

namespace BetheAnsatz {

template<typename RealType>
class GrandPotential {

public:

	GrandPotential(RealType, RealType)
	{}

	RealType at(RealType)
	{
		return 0.0;
	}
}; // class GrandPotential

} // namespace BetheAnsatz
#endif // GRANDPOTENTIAL_H

