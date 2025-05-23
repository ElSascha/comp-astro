#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include <fstream>
#include "lib/Body.hpp"
void print_data_in_file(Body& body, const int max_ite , const std::string& filename){
    std::ofstream outfile("data/"+filename);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    outfile << "time (days); x ; y ; r ; f ; E ; M ; P" << std::endl;
    for (int i = 0; i < max_ite; ++i) {
                outfile << body.get_t() << " ; "
                << body.get_x() << " ; "
                << body.get_y() << " ; "
                << body.get_r() << " ; "
                << body.get_f() << " ; "
                << body.get_E() << " ; "
                << body.get_M() << " ; "
                << body.get_P() << std::endl;
        body.update(max_ite, true);
    }
    outfile.close();
}



int main(){
    Body mercury = Body(0.205, 0.39, 0);
    Body halleys_comet = Body(0.967, 17.8, 0);
    print_data_in_file(mercury,1000, "mercury_1000_ite.txt");

    return 0;
}