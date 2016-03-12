#ifndef FOURSPINONHEISENBERG_H
#define FOURSPINONHEISENBERG_H
#include "TwoSpinonHeisenberg.h"

namespace BetheAnsatz {

// J.-S. Caux, R. Hagemans, J. Stat. Mech. 0612:P12013 (2006)
template<typename RealType_>
class FourSpinonHeisenberg {

	typedef typename PsimagLite::Vector<RealType_>::Type VectorRealType;
	typedef TwoSpinonHeisenberg<RealType_> TwoSpinonHeisenbergType;
	typedef Mesh<RealType_> MeshType;

	class InnerIntegrand {

	public:

		enum {F1, F2};

		struct Params {
			Params(RealType_ k_,
			       RealType_ w_,
			       RealType_ K_,
			       const MeshType& kmesh,
			       const MeshType& emesh)
			    : k(k_),w(w_),K(K_),tsh(kmesh,emesh)
			{}

			RealType_ k;
			RealType_ w;
			RealType_ K;
			TwoSpinonHeisenbergType tsh;
		};

		typedef RealType_ RealType;

		InnerIntegrand(RealType k,
		               RealType w,
		               RealType K,
		               const MeshType& kmesh,
		               const MeshType& emesh)
		    : params_(k,w,K,kmesh,emesh)
		{}

		Params& params() { return params_; }

		static RealType function(RealType W, void* ptr)
		{
			Params* paramsPtr = static_cast<Params*>(ptr);
			RealType k = paramsPtr->k;
			RealType w = paramsPtr->w;
			RealType K = paramsPtr->K;

			RealType tmp = w2u(K);
			RealType den = (tmp*tmp-W*W);
			RealType tmp2 = w2u(k-K);
			RealType tmp3 = w - W;
			den *= (tmp2*tmp2-tmp3*tmp3);
			return jFunction(k,w,K,W,paramsPtr->tsh)/sqrt(den);
		}

	private:

		static RealType w2u(RealType k)
		{
			return M_PI*sin(k*0.5);
		}

		static RealType jFunction(RealType k,
		                          RealType w,
		                          RealType K,
		                          RealType W,
		                          const TwoSpinonHeisenbergType& tsh)
		{
			VectorRealType p(4,0.0);
			// get p1 and p2
			getPlow(p,K,W);
			// get p3 and p4
			getPhigh(p,k-K,w-W);
			// get rho
			VectorRealType rho(4,0.0);
			RealType factor = 0.5/M_PI;
			for (SizeType i = 0; i < p.size(); ++i) {
				RealType tmp = cot(p[i]);
				rho[i] = tmp + sqrt(tmp*tmp + 1.0);
				rho[i] *= factor;
			}

			return jOfRho(rho, tsh);
		}

		static RealType jOfRho(const VectorRealType& rho, const TwoSpinonHeisenbergType& tsh)
		{
			RealType prod = 1.0;
			for (SizeType i = 0; i < rho.size(); ++i) {
				for (SizeType j = i + 1; j < rho.size(); ++j) {
					RealType rhoij = rho[i] - rho[j];
					prod *= tsh.eMinusIfunction(rhoij);
				}
			}

			RealType sum2 = 0.0;
			for (SizeType i = 0; i < rho.size(); ++i)
				sum2 += computeGl2(i,rho);

			return prod*sum2;
		}

		static void getPlow(VectorRealType& p, RealType K, RealType W)
		{
			RealType x = -0.5*K;
			RealType y = acos(W/w2u(K));
			p[0] = x + y;
			p[1] = x - y;
		}

		static void getPhigh(VectorRealType& p, RealType kdiff, RealType wdiff)
		{
			bool sector1 = (kdiff < 0);
			RealType x = -0.5*kdiff;
			if (sector1) {
				x -= M_PI;
				kdiff *= (-1.0);
			}

			RealType y = acos(wdiff/w2u(kdiff));
			p[2] = x + y;
			p[3] = x - y;
		}

		static RealType cot(RealType x)
		{
			return 1.0/tan(x);
		}

		Params params_;
	}; // class InnerIntegrand

	class Eq35Integrand {

	public:

		enum {F1, F2};

		struct Params {
			Params(RealType_ k_,
			       RealType_ w_,
			       const MeshType& kmesh_,
			       const MeshType& emesh_)
			    : k(k_),w(w_),kmesh(kmesh_),emesh(emesh_)
			{}

			RealType_ k;
			RealType_ w;
			const MeshType& kmesh;
			const MeshType& emesh;
		};

		typedef RealType_ RealType;

		Eq35Integrand(RealType k,
		              RealType w,
		              const MeshType& kmesh,
		              const MeshType& emesh)
		    : params_(k,w,kmesh,emesh)
		{}

		Params& params() { return params_; }

		static RealType function(RealType K, void* ptr)
		{
			Params* paramsPtr = static_cast<Params*>(ptr);
			RealType k = paramsPtr->k;
			RealType w = paramsPtr->w;

			VectorRealType pts(2,0.0);
			pts[0] = omegaL(k,w,K);
			pts[1] = omegaU(k,w,K);
			InnerIntegrand innerIntegrand(k,w,K,paramsPtr->kmesh,paramsPtr->emesh);
			PsimagLite::Integrator<InnerIntegrand> integrator(innerIntegrand);
			return integrator(pts);
		}

	private:

		static RealType omegaL(RealType k, RealType w, RealType K)
		{
			RealType x1 = M_PI*0.5*fabs(sin(K));
			RealType x2 = w - M_PI*fabs(sin(0.5*(k-K)));
			return std::max(x1,x2);
		}

		static RealType omegaU(RealType k, RealType w, RealType K)
		{
			RealType x1 = M_PI*sin(K*0.5);
			RealType x2 = w - M_PI*0.5*fabs(sin(k-K));
			return std::min(x1,x2);
		}

		Params params_;
	}; // class Eq35Integrand

public:

	typedef RealType_ RealType;

	FourSpinonHeisenberg(const BetheAnsatz::Mesh<RealType>& kmesh,
	                     const BetheAnsatz::Mesh<RealType>& emesh)
	    : kmesh_(kmesh),
	      emesh_(emesh),
	      m_(kmesh.total(),emesh.total())
	{
		// fixme: C4
		RealType c4 = 1.0;
		for (SizeType i = 0; i < kmesh.total(); ++i) {
			RealType k = kmesh.x(i);
			for (SizeType j = 0; j < emesh.total(); ++j) {
				RealType w = emesh.x(j);
				// Eq. (1.7) of Ref [1]
				m_(i,j) = c4*mainFunction(k,w);
			}
		}
	}

private:

	RealType mainFunction(RealType k,RealType w)
	{
		Eq35Integrand eq5Integrand(k,w,kmesh_,emesh_);
		VectorRealType pts(2,0.0);
		pts[1] = 2.0*M_PI; // FIXME: limits of integration on many sectors
		PsimagLite::Integrator<Eq35Integrand> integrator(eq5Integrand);
		return integrator(pts);
	}

	const BetheAnsatz::Mesh<RealType>& kmesh_;
	const BetheAnsatz::Mesh<RealType>& emesh_;
	PsimagLite::Matrix<RealType> m_;
}; // class FourSpinonHeisenberg
} // namespace BetheAnsatz

#endif // FOURSPINONHEISENBERG_H
