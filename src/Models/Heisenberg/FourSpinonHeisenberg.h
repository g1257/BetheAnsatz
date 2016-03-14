#ifndef FOURSPINONHEISENBERG_H
#define FOURSPINONHEISENBERG_H
#include "TwoSpinonHeisenberg.h"
#include "SpecialFunctions.h"

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
					RealType rhoij = fabs(rho[i] - rho[j]);
					if (rhoij == 0) return 0.0;
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

		static RealType computeGl2(SizeType l, const VectorRealType& rho)
		{
			std::complex<RealType> sum1 = 0.0;
			for (SizeType j = 0; j < 4; ++j)
				sum1 += cosh(2.0*M_PI*rho[j])*computeGl2Aux(l,rho,j);

			std::complex<RealType> sum = std::conj(sum1)*sum1;
			return (l&1) ? -std::real(sum) : std::real(sum);
		}

		static std::complex<RealType> computeGl2Aux(SizeType l,
		                                            const VectorRealType& rho,
		                                            SizeType j)
		{
			std::complex<RealType> sum1 = 0.0;
			std::complex<RealType> summand = 0.0;
			std::complex<RealType> prev = 0.0;
			SizeType start = (j <= l) ? 0 : 1;
			SizeType mMax = 100;
			for (SizeType m = start; m < mMax; ++m) {
				prev = summand;
				summand = prodOne(l,rho,j,m)*prodTwo(rho,j,m);
				sum1 += summand;
			}

			if (std::norm(prev-summand)>1e-6)
				throw PsimagLite::RuntimeError("computeGl2Aux failed\n");

			return sum1;
		}

		static std::complex<RealType> prodOne(SizeType l,
		                                      const VectorRealType& rho,
		                                      SizeType j,
		                                      SizeType m)
		{
			std::complex<RealType> prod1 = 1.0;
			for (SizeType i = 0; i < 4; ++i) {
				if (i == l && i == j) continue;
				std::complex<RealType> termNum(m,rho[j] - rho[i]);
				if (l > i) termNum -= 0.5;
				RealType termDen = sinh(M_PI*(rho[j] - rho[i]));
				std::complex<RealType> term = (i != l) ?  termNum : 1.0;
				if (i != j) term /= termDen;
				prod1 *= term;
			}

			return prod1;
		}

		static std::complex<RealType> prodTwo(const VectorRealType& rho,
		                                      SizeType j,
		                                      SizeType m)
		{
			std::complex<RealType> prod1 = 1.0;
			for (SizeType i = 0; i < 4; ++i) {
				RealType rhoji = rho[j] - rho[i];
				std::complex<RealType> g1 = PsimagLite::LnGammaFunction(
				            std::complex<RealType>(m-0.5,rhoji));

				RealType r1 = std::real(g1);
				RealType a1 = std::imag(g1);
				std::complex<RealType> c1(r1*cos(a1),r1*sin(a1));

				std::complex<RealType> g2 = PsimagLite::LnGammaFunction(
				            std::complex<RealType>(m+1.0,rhoji));

				RealType r2 = std::real(g2);
				RealType a2 = std::imag(g2);
				std::complex<RealType> c2(r2*cos(a2),r2*sin(a2));
				prod1 *= exp(c1/c2);
			}

			return prod1;
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

			if (inValidRegion(k,w,K)) return 0.0;

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

		static bool inValidRegion(RealType k, RealType w, RealType K)
		{
			RealType piSinkOver2 = M_PI*sin(0.5*k);
			RealType piCoskOver2 = M_PI*cos(0.5*k);
			RealType piCoskOver4 = M_PI*cos(k*0.25);
			RealType piSinkOver4 = M_PI*sin(k*0.25);

			RealType x = 0.5*k + M_PI;
			RealType y = 2.0*acos(0.5*w/piCoskOver4);
			RealType k1am = x - y;
			RealType k1ap = x + y;

			x = 0.5*k;
			y = 2.0*acos(0.5*w/piSinkOver4);
			RealType k1bm = x - y;
			RealType k1bp = x + y;

			if (w > piSinkOver2 && w <= 2.0*piSinkOver4) {
				bool b1 = (K >= k1am && K <= k1ap);
				bool b2 = (K >= k1bm && K <= k1bp);
				if (!b1 && !b2) return true;
			}

			RealType piOver2Sink = M_PI*0.5*sin(k);

			x = 0.5*(k+M_PI);
			y = acos(w/piCoskOver2);
			RealType k2cm = x - y;
			RealType k2cp = x + y;

			if (w >= piOver2Sink && w <= piSinkOver2) {
				bool b1 = (K >= k2cm && K <= k2cp);
				bool b2 = (K >= k2cm + M_PI && K <= k2cp+M_PI);
				if (b1 || b2) return true;
			}

			x  = 0.5*k;
			y = acos(w/piSinkOver2);
			RealType k2dm = x - y;
			RealType k2dp = x + y;
			if (w >= piOver2Sink && w <= piCoskOver2) {
				bool b1 = (K >= k2dm && K <= k2dp);
				bool b2 = (K >= k2dm + M_PI && K <= k2dp+M_PI);
				if (b1 || b2) return true;
			}

			return false;
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
				if (k <= 0 || k >= M_PI) continue;
				if (w <= w4l(k) || w >= w4u(k)) continue;
				// Eq. (1.7) of Ref [1]
				m_(i,j) = c4*mainFunction(k,w);
			}
		}
	}

	void toGnuplot(std::ostream& os, RealType cutoff)
	{
		for (SizeType j = 0; j < emesh_.total(); ++j) {
			RealType w = emesh_.x(j);
			for (SizeType i = 0; i < kmesh_.total(); ++i) {
				RealType k = kmesh_.x(i);
				RealType value = m_(i,j);
				if (cutoff > 0)
					value = (value > cutoff) ? cutoff : value;
				os<<k<<" "<<w<<" "<<value<<"\n";
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

	RealType w4l(RealType k) const
	{
		return 0.5*M_PI*fabs(sin(k));
	}


	RealType w4u(RealType k) const
	{
		return M_PI*sqrt(2.0*(1.0+fabs(cos(0.5*k))));
	}

	const BetheAnsatz::Mesh<RealType>& kmesh_;
	const BetheAnsatz::Mesh<RealType>& emesh_;
	PsimagLite::Matrix<RealType> m_;
}; // class FourSpinonHeisenberg
} // namespace BetheAnsatz

#endif // FOURSPINONHEISENBERG_H
