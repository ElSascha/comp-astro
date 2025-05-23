#include "Body.hpp"
#include "Method.cpp" // Include the Method class for calculations

double Body::calc_P()
{
    return sqrt(pow(this->a, 3)) * 365; // in days
}

void Body::update_M()
{
    this->M = 2 * M_PI * P * (this->t - this->t_0);
}

void Body::update_x()
{
    this->x = this->r * cos(this->f);
}
void Body::update_y()
{
    this->y = this->r * sin(this->f);
}
void Body::update_f()
{
    this->f = 2 * atan(sqrt((1 + this->e) / (1 - this->e)) * tan(this->E / 2));
}
void Body::update_r()
{
    this->r = this->a * (1 - pow(this->e, 2)) / (1 + this->e * cos(this->f));
}
void Body::update_E(int max_ite, bool newton)
{
    if (newton)
    {
        this->E = Method::newton_raphson(max_ite, this->M, this->e);
    }
    else
    {
        this->E = Method::fix_point(max_ite, this->M, this->e);
    }
}

Body::Body(const double &e, const double &a, const double &t_0)
{
    this->e = e;
    this->a = a;
    this->t_0 = t_0;
    this->x = 0;
    this->y = 0;
    this->P = calc_P();
    this->M = 0;
    this->f = 0;
    this->r = 0;
    if (e < 0.8)
    {
        this->E = 0;
    }
    else
    {
        this->E = M_PI;
    }
}

void Body::update(int max_ite, bool newton)
{
    this->t = this->t + 1;
    update_M();
    update_E(max_ite, newton);
    update_f();
    update_r();
    update_x();
    update_y();
}
double Body::get_x()
{
    return this->x;
}
double Body::get_y()
{
    return this->y;
}
double Body::get_r()
{
    return this->r;
}
double Body::get_f()
{
    return this->f;
}
double Body::get_E()
{
    return this->E;
}
double Body::get_M()
{
    return this->M;
}
double Body::get_P()
{
    return this->P;
}
double Body::get_t()
{
    return this->t;
}

// The Body class represents a celestial body in an elliptical orbit.
// It contains methods to calculate its position and update its state over time.
// The class uses the Kepler's equation to calculate the position of the body in its orbit.
// The class also uses the Method class to perform the calculations for the position of the body.
// The class contains private member variables to store the state of the body, such as its eccentricity, semi-major axis, mean anomaly, and position.
// The class also contains public methods to update the state of the body and get its position.
// The class is designed to be used in a simulation of celestial bodies in an elliptical orbit.
