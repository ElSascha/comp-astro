#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include <fstream>
#include "lib/Body.hpp"
#include "lib/Distance.cpp"
#include <chrono>
#include <format> // For std::format (C++20)
#include <iomanip>


void print_data_in_file(Body& body,const int max_time , const std::string& filename, bool newton){
    std::ofstream outfile("data/"+filename);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    outfile << "time (days); x ; y ; r ; f ; E ; M ; P" << std::endl;
    for (int i = 0; i < max_time; ++i) {
        body.update(newton);
                outfile << body.get_t() << " ; "
                << body.get_x() << " ; "
                << body.get_y() << " ; "
                << body.get_r() << " ; "
                << body.get_f() << " ; "
                << body.get_E() << " ; "
                << body.get_M() << " ; "
                << body.get_P() << std::endl;
        
    }
    outfile.close();
}

void print_distance_in_file(Body& earth, Body& mars, const std::string& filename){
    using namespace std::chrono;
    std::ofstream outfile("data/"+filename);
    if (!outfile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    outfile << "Date (days); distance (AU)" << std::endl;
    std::tm tm_start = {};
    tm_start.tm_mday = 1;
    tm_start.tm_mon = 0; // January (0-based)
    tm_start.tm_year = 1985 - 1900;
    std::tm tm_end = {};
    tm_end.tm_mday = 29;
    tm_end.tm_mon = 4; // May (0-based)
    tm_end.tm_year = 2025 - 1900;
    std::time_t time_start = std::mktime(&tm_start);
    std::time_t time_end = std::mktime(&tm_end);
    for(std::time_t t = time_start; t <= time_end; t += 24*60*60) {
        std::tm* tm_ptr = std::localtime(&t);
        char buf[11];
        std::strftime(buf, sizeof(buf), "%d.%m.%Y", tm_ptr);
        earth.update(true, Distance::earth_M_0);
        mars.update(true, Distance::mars_M_0);
        double distance = std::sqrt(std::pow(earth.get_x() - mars.get_x(), 2) + std::pow(earth.get_y() - mars.get_y(), 2));
        outfile << buf << " ; " << distance << std::endl;
        
    }
}


int main(){
    Body mercury = Body(0.205, 0.39, 0, 256);
    Body halleys_comet = Body(0.967, 17.8, 0, 256);
    Body earth = Body(Distance::earth_e, Distance::earth_a);
    Body mars = Body(Distance::mars_e, Distance::mars_a);
    // std::cout << "Mercury fixed point:" << std::endl;
    // print_data_in_file(mercury,mercury.get_data_points(), "mercury_fixed_point.txt", false);
    // std::cout << "Halley's Comet fixed point :" << std::endl;
    // print_data_in_file(halleys_comet,halleys_comet.get_data_points(), "halleys_comet_fixed_point.txt",false);
    // std::cout << "Mercury Newton-Raphson:" << std::endl;
    // print_data_in_file(mercury,mercury.get_data_points(), "mercury_newton_raphson.txt", true);
    // std::cout << "Halley's Comet Newton-Raphson:" << std::endl;
    // print_data_in_file(halleys_comet,halleys_comet.get_data_points(), "halleys_comet_newton_raphson.txt", true);
    print_distance_in_file(earth, mars, "distance_earth_mars.txt");
    return 0;
}