#ifndef BODY_H
#define BODY_H

#include <string>
#include <array>
#include <iostream>
#include <cmath>
#include <vector>

class Body {
public:
    Body(const double& e, const double&a, const double& t_0, const int& data_points);
    Body(const double& e, const double&a);
    void update(bool newton);
    void update(bool newton, const double& M_0);
    double get_x();
    double get_y();
    double get_r();
    double get_f();
    double get_E();
    double get_M();
    double get_P();
    double get_t();
    int get_data_points();



private:
    double x,y,e,a,M,P,t,t_0,f,r,E,dt;
    int data_points;
    double calc_P();
    void update_M();
    void update_M(const double& M_0);
    void update_x();
    void update_y();
    void update_f();
    void update_r();
    void update_E(bool newton);

    
};

#endif // BODY_H