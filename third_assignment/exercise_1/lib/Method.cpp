#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include <fstream>
#include <iomanip>

namespace Method
{
    
    static double fix_point_worker(double M, double e, double E_i, int i)
    {
        
        double E_ipo = M + e * sin(E_i);
        if(std::abs(E_ipo - E_i) <= 1e-9)
        {
            std::cout<<i << " iterations" << std::endl;
            return E_ipo;
        }
        return fix_point_worker(M, e, E_ipo, i + 1);
    }

    static double fix_point(double M, double e)
    {
        if (e < 0.8)
        {
            return fix_point_worker(M, e, M, 0);
        }
        else
        {
            return fix_point_worker(M, e, M_PI, 0);
        }
    }

    

     static double newton_raphson_worker(double M, double e, double E_i, int i)
    {
        
        double g = E_i - e * sin(E_i) - M;
        double dg = 1 - e * cos(E_i);
        double E_ipo = E_i - g / dg;
        if(std::abs(E_ipo - E_i) <= 1e-9)
        {
            return E_ipo;
        }
        return newton_raphson_worker(M, e, E_ipo, i + 1);
    }

    static double newton_raphson(double M, double e)
    {
        if (e < 0.8)
        {
            return newton_raphson_worker(M, e, M, 0);
        }
        else
        {
            return newton_raphson_worker(M, e, M_PI, 0);
        }
    }

   

}