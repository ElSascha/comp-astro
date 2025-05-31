#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include <fstream>
#include <iomanip>

namespace Distance{
    // Data of earth and mars on the 01.01.2000
    const double earth_a = 1.0; // semi-major axis in AU
    const double earth_e = 0.0167; // eccentricity 
    const double earth_phi_0 =  102.95 * M_PI/180.0 ; // mean anomaly in radians
    const double earth_lambda = 100.46 * M_PI/180.0; // mean longitude of the orbit in radians
    const double mars_a = 1.524; // semi-major axis in AU
    const double mars_e = 0.0934; // eccentricity
    const double mars_phi_0 =  336.04 * M_PI/180.0; // mean anomaly in radians
    const double mars_lambda = 355.45 * M_PI/180.0; // mean longitude of the orbit in radians
    const double earth_P = 365; // orbital period of earth in days
    const double mars_P = 687; // orbital period of mars in days

    // Calculate M_0 on the 01.01.1985
    const double earth_M_0 = earth_lambda -earth_phi_0 - ((15.0*365.0) * (2.0*M_PI)/(earth_P));
    const double mars_M_0 = mars_lambda - mars_phi_0 - ((15.0*365.0) * (2.0*M_PI)/(mars_P));

}