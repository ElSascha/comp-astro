#ifndef BODY_H
#define BODY_H

#include <string>
#include <array>
#include <iostream>
#include <cmath>
#include <vector>

class Body {
public:
    Body(const double& e, const double&a, const double& t_0);
    void update(int max_ite, bool newton);
    double get_x();
    double get_y();
    double get_r();
    double get_f();
    double get_E();
    double get_M();
    double get_P();
    double get_t();



private:
    double x,y,e,a,M,P,t,t_0,f,r,E;
    double calc_P();
    void update_M();
    void update_x();
    void update_y();
    void update_f();
    void update_r();
    void update_E(int max_ite, bool newton);
    
};

#endif // BODY_H