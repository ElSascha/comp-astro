#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include <fstream>

namespace Method
{
    static double fix_point_worker(int max_ite, double M, double e, double E_i, int i)
    {
        if (i == max_ite)
        {
            return E_i;
        }
        double E_ipo = M + e * sin(E_i);
        return fix_point_worker(max_ite, M, e, E_ipo, i + 1);
    }

    static double fix_point(int max_ite, double M, double e)
    {
        if (e < 0.8)
        {
            return fix_point_worker(max_ite, M, e, M, 0);
        }
        else
        {
            return fix_point_worker(max_ite, M, e, M_PI, 0);
        }
    }

    

     static double newton_raphson_worker(int max_ite, double M, double e, double E_i, int i)
    {
        if (i == max_ite)
        {
            return E_i;
        }
        double g = E_i - e * sin(E_i) - M;
        double dg = 1 - e * cos(E_i);
        double E_ipo = E_i - g / dg;
        return newton_raphson_worker(max_ite, M, e, E_ipo, i + 1);
    }

    static double newton_raphson(int max_ite, double M, double e)
    {
        if (e < 0.8)
        {
            return newton_raphson_worker(max_ite, M, e, M, 0);
        }
        else
        {
            return newton_raphson_worker(max_ite, M, e, M_PI, 0);
        }
    }

   

}