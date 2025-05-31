#include "Body.hpp"
#include "Method.cpp" // Include the Method class for calculations

double Body::calc_P()
{
    return sqrt(pow(this->a, 3)) * 365; // in days
}

void Body::update_M()
{
    this->M = 2 * M_PI * (this->t - this->t_0) / P;
}
void Body::update_M(const double &M_0)
{
    this->M = (2 * M_PI * (this->t) / P) + M_0;
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
void Body::update_E(bool newton)
{
    if (newton)
    {
        this->E = Method::newton_raphson(this->M, this->e);
    }
    else
    {
        this->E = Method::fix_point(this->M, this->e);
    }
}

Body::Body(const double &e, const double &a, const double &t_0, const int &data_points)
{
    this->e = e;
    this->a = a;
    this->t_0 = t_0;
    this-> data_points = data_points;
    this->P = calc_P();
    this->dt = this->P / this->data_points;
    this->t = this->t_0;
    this->M = 0;
    this->f = 0;
    this->r = 0;
    this->x = 0;
    this->y = 0;
    // Initialize E based on eccentricity
    // If eccentricity is less than 0.8, we can start with M, otherwise we start with M_PI
    if (e < 0.8)
    {
        this->E = 0;
    }
    else
    {
        this->E = M_PI;
    }
}

Body::Body(const double &e, const double &a)
{
    this->e = e;
    this->a = a;
    this->t_0 = 0;
    this-> data_points = 0;
    this->P = calc_P();
    this->dt = 0.0;
    this->t = 0.0;
    this->M = 0.0;
    this->f = 0.0;
    this->r = 0.0;
    this->x = 0.0;
    this->y = 0.0;
    // Initialize E based on eccentricity
    // If eccentricity is less than 0.8, we can start with M, otherwise we start with M_PI
    if (e < 0.8)
    {
        this->E = 0;
    }
    else
    {
        this->E = M_PI;
    }
}



void Body::update(bool newton)
{
    this->t = this->t + this->dt;
    update_M();
    update_E(newton);
    update_f();
    update_r();
    update_x();
    update_y();
}
void Body::update(bool newton,const double &M_0)
{
    this->t = this->t + 1;
    update_M(M_0);
    update_E(newton);
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
    std::cout << "M: " << this->M << std::endl;
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
int Body::get_data_points()
{
    return this->data_points;
}

// The Body class represents a celestial body in an elliptical orbit.
// It contains methods to calculate its position and update its state over time.
// The class uses the Kepler's equation to calculate the position of the body in its orbit.
// The class also uses the Method class to perform the calculations for the position of the body.
// The class contains private member variables to store the state of the body, such as its eccentricity, semi-major axis, mean anomaly, and position.
// The class also contains public methods to update the state of the body and get its position.
// The class is designed to be used in a simulation of celestial bodies in an elliptical orbit.
