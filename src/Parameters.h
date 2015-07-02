#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "Vector.h"

namespace BetheAnsatz {

template<typename RealType_, typename InputType>
struct Parameters {

	typedef RealType_ RealType;

	Parameters(InputType& io)
	    : meshKtotal(1000),
	      meshLambdaTotal(10000),
	      nMax(50),
	      tt(0.0),
	      iterations(1),
	      infty(1e6),
	      tb(0.0),
	      ts(0.0),
	      U(0.0)
	{
		io.readline(nMax,"Nmax=");
		io.readline(tt,"TemperatureTotal=");
		io.readline(iterations,"Iterations=");
		io.readline(tb,"TemperatureBegin=");
		io.readline(ts,"TemperatureStep=");
		io.readline(U,"U=");

		try {
			io.readline(meshKtotal,"MeshKTotal=");
		} catch (std::exception&) {}

		try {
			io.readline(meshLambdaTotal,"MeshLambdaTotal=");
		} catch (std::exception&) {}

		try {
			io.readline(infty,"Infty=");
		} catch (std::exception&) {}
	}

	SizeType meshKtotal;
	SizeType meshLambdaTotal;
	SizeType nMax;
	SizeType tt;
	SizeType iterations;
	RealType infty;
	RealType tb;
	RealType ts;
	RealType U;
}; // struct Parameters

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, const Parameters<T1,T2>& params)
{
	os<<"MeshKTotal="<<params.meshKtotal<<"\n";
	os<<"MeshLambdaTotal="<<params.meshLambdaTotal<<"\n";
	os<<"Nmax="<<params.nMax<<"\n";
	os<<"TemperatureTotal="<<params.tt<<"\n";
	os<<"Iterations="<<params.iterations<<"\n";
	os<<"Infty="<<params.infty<<"\n";
	os<<"TemperatureBegin="<<params.tb<<"\n";
	os<<"TemperatureStep="<<params.ts<<"\n";
	os<<"U="<<params.U<<"\n";

	return os;
}
} // namespace BetheAnsatz
#endif // PARAMETERS_H

