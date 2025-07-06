#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include <fstream>
using namespace std;

namespace Methods{
    vector<double> centered_differencing(const std::vector<double>&, double);
    vector<double> upwind(const std::vector<double>&, double);
    vector<double> lax_wendroff(const std::vector<double>&, double);
}