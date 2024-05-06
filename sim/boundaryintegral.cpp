// adapted from Fang Da's 'Surface-only Liquids' source code BoundaryIntegral.cpp

#include "boundaryintegral.h"
#include <map>

const std::vector<Eigen::Vector2d> BoundaryIntegral::gaussian_quadrature(int N){

    static std::map<int, std::vector<Eigen::Vector2d> > s_quadratures;

    std::map<int, std::vector<Eigen::Vector2d> >::iterator i = s_quadratures.find(N);
    
    if (i == s_quadratures.end())
    {
        std::vector<Eigen::Vector2d> q;
        switch (N)
        {
            // N = 1~4: Gaussian quadrature rules
            case 1:
                q.push_back(Eigen::Vector2d(0, 2));
                break;
            case 2:
                q.push_back(Eigen::Vector2d(-sqrt(1.0 / 3.0), 1));
                q.push_back(Eigen::Vector2d( sqrt(1.0 / 3.0), 1));
                break;
            case 3:
                q.push_back(Eigen::Vector2d(-sqrt(3.0 / 5.0), 5.0/9.0));
                q.push_back(Eigen::Vector2d(               0, 8.0/9.0));
                q.push_back(Eigen::Vector2d( sqrt(3.0 / 5.0), 5.0/9.0));
                break;
            case 4:
                q.push_back(Eigen::Vector2d(-sqrt(3.0/7.0+2.0/7.0*sqrt(6.0/5.0)), (18.0-sqrt(30))/36.0));
                q.push_back(Eigen::Vector2d(-sqrt(3.0/7.0-2.0/7.0*sqrt(6.0/5.0)), (18.0+sqrt(30))/36.0));
                q.push_back(Eigen::Vector2d( sqrt(3.0/7.0-2.0/7.0*sqrt(6.0/5.0)), (18.0+sqrt(30))/36.0));
                q.push_back(Eigen::Vector2d( sqrt(3.0/7.0+2.0/7.0*sqrt(6.0/5.0)), (18.0-sqrt(30))/36.0));
                break;
            default:
                // N > 4: uniform grid of N points (midpoint rule)
                for (int j = 0; j < N; j++)
                    q.push_back(Eigen::Vector2d((j + 0.5) / N * 2 - 1, 2.0 / N));
                break;
        }
        
        std::pair<std::map<int, std::vector<Eigen::Vector2d> >::iterator, bool> res = s_quadratures.insert(std::pair<int, std::vector<Eigen::Vector2d> >(N, q));
        assert(res.second);
        i = res.first;
    }

    return i->second;
}
const std::vector<Eigen::Vector2d> BoundaryIntegral::tanh_sinh_quadrature(int N){
    static std::map<int, std::vector<Eigen::Vector2d> > s_quadratures;

    assert(N%2==1);

    std::map<int, std::vector<Eigen::Vector2d> >::iterator i = s_quadratures.find(N);
    
    if (i == s_quadratures.end())
    {
        int n = (N-1)/2;
        double h = 0.3;

        std::vector<Eigen::Vector2d> q;

        for (int j=-n; j<=n; j++){
            //q.push_back(Eigen::Vector2d((1.0/2.0)*(1+tanh((M_PI/2.0)*sinh(j*h))), h*(M_PI/4.0)*cosh(j*h)/(pow(cosh((M_PI/2.0)*sinh(j*h)),2)))); //[0,1]
            q.push_back(Eigen::Vector2d(tanh((M_PI/2.0)*sinh(j*h)), (2.0)*h*(M_PI/4.0)*cosh(j*h)/(pow(cosh((M_PI/2.0)*sinh(j*h)),2))));
        }
        
        std::pair<std::map<int, std::vector<Eigen::Vector2d> >::iterator, bool> res = s_quadratures.insert(std::pair<int, std::vector<Eigen::Vector2d> >(N, q));
        assert(res.second);
        i = res.first;
    }

    return i->second;
}