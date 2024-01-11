#include "../../include/tinynurbs/tinynurbs.h"

#include <algorithm>
#include <vector>

namespace zlz{
	template <typename T>
	int findSpanLinear(unsigned int degree, const std::vector<T>& knots, T u) {
        auto span = degree + 1;  // Knot span index starts from zero
        auto num_ctrlpts= knots.size() - degree - 1;
        while (span < num_ctrlpts and knots[span] <= u){
            span ++;
        }
        return span - 1;
	}

	template <typename T>
	std::vector<T> bsplineBasis(unsigned int deg, int span, const std::vector<T>& knots, T u) {
		std::vector<T> N(deg+1,1),left(deg+1,0), right(deg+1,0);
        for(auto j=1;j<deg+1;j++){
            left[j] = u - knots[span+1-j];
            right[j] = knots[span+j] - u;
            T saved = 0;
            for(auto r=0;r<j;r++){
                T temp = N[r] / (right[r+1] + left[j-r]);
                N[r] = saved + right[r+1] * temp;
                saved = left[j-r] * temp;
            }
            N[j] = saved;
        }

		return N;
	}

	template <typename T>
	std::vector<std::vector<T>> bsplineDerBasis(unsigned int deg, int span, const std::vector<T>& knots, T u,
		int num_ders) {
        auto left= std::vector<T>(deg+1,0),right= std::vector<T>(deg+1,0);
        auto ndu= std::vector<std::vector<T>>(deg+1, std::vector<T>(deg+1,0));

        for(auto j=1;j<deg+1;j++){
            left[j] = u - knots[span+1-j];
            right[j] = knots[span+j] - u;
            T saved = 0;
            for(auto r=0;r<j;r++){
                ndu[j][r] = right[r+1] + left[j-r];
                T temp = ndu[r][j-1] / ndu[j][r];
                ndu[r][j] = saved + right[r+1] * temp;
                saved = left[j-r] * temp;
            }
            ndu[j][j] = saved;
        }

        auto ders= std::vector<std::vector<T>>(num_ders+1, std::vector<T>(deg+1,0));
        for(auto j=0;j<deg+1;j++){
            ders[0][j] = ndu[j][deg];
        }

        auto a= std::vector<std::vector<T>>(2, std::vector<T>(deg+1,0));
        for(auto r=0;r<deg+1;r++){
            auto s1=0, s2=1;
            a[0][0] = 1.0;
            for(auto k=1;k<num_ders+1;k++){
                T d = 0.0;
                auto rk = r - k;
                auto pk = deg - k;
                if(r >= k){
                    a[s2][0] = a[s1][0] / ndu[pk+1][rk];
                    d = a[s2][0] * ndu[rk][pk];
                }
                auto j1 = 0, j2 = 0;
                if(rk >= -1){
                    j1 = 1;
                }else{
                    j1 = -rk;
                }
                if(r-1 <= pk){
                    j2 = k - 1;
                }else{
                    j2 = deg - r;
                }
                for(auto j=j1;j<j2+1;j++){
                    a[s2][j] = (a[s1][j] - a[s1][j-1]) / ndu[pk+1][rk+j];
                    d += a[s2][j] * ndu[rk+j][pk];
                }
                if(r <= pk){
                    a[s2][k] = -a[s1][k-1] / ndu[pk+1][r];
                    d += a[s2][k] * ndu[r][pk];
                }
                ders[k][r] = d;
                auto temp = s1;
                s1 = s2;
                s2 = temp;
            }
        }

        double r=deg;
        for(auto k=1;k<num_ders+1;k++){
            for(auto j=0;j<deg+1;j++){
                ders[k][j] *= r;
            }
            r *= static_cast<T>(deg - k);
        }
        return ders;
	}

    template <int dim, typename T>
    glm::vec<dim, T> curvePoint(unsigned int degree, const std::vector<T>& knots,
        const std::vector<glm::vec<dim, T>>& control_points, T u) {
        // Initialize result to 0s
        glm::vec<dim, T> point(T(0));

        // Find span and corresponding non-zero basis functions
        int span = zlz::findSpanLinear(degree, knots, u);
        std::vector<T> N = zlz::bsplineBasis(degree, span, knots, u);

        // Compute point
        for(unsigned int j = 0; j <= degree; j++) {
            point += static_cast<T>(N[j]) * control_points[span - degree + j];
        }
        return point;
    }

    template <int dim, typename T>
    std::vector<glm::vec<dim, T>> curveDerivatives(unsigned int degree, const std::vector<T>& knots,
        const std::vector<glm::vec<dim, T>>& control_points,
        int num_ders, T u) {

        auto du = std::min(degree, static_cast<unsigned>(num_ders));
        auto CK = std::vector<glm::vec<dim, T>>(du + 1, glm::vec<dim, T>(0));
        auto span= zlz::findSpanLinear(degree, knots, u);
        auto ders= zlz::bsplineDerBasis(degree, span, knots, u, du);

        for(auto k=0;k<du+1;k++){
            for(auto j=0;j<degree+1;j++){
                CK[k] += static_cast<T>(ders[k][j]) * control_points[span - degree + j];
            }
        }
    }
}