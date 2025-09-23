#include <iostream>
#include<vector>
#include<cmath>
#include <numeric>
#include <fstream>   //for saving to computer
#include <algorithm>
#include <random>
#include <any>

const float G =1;    // gravitational constant, mass units, and distance units are all normalised
const float pi = 3.14159265358979323846;   //pi

/**
 * @brief  Function for finding the cartesian acceleration of a star in a galaxy, when feeling a force 
 * from the galaxy centre of mass and a point perturber.
 * @param  x:  x position of the star
 * @param  y:  y position of the star
 * @param x_c: x position of the galaxy centre of mass
 * @param y_c: y position of the galaxy centre of mass
 * @param x_p: x position of the Moon
 * @param y_p: y position of the Moon
 * @param M: mass of the galaxy
 * @param m_p: mass of the point perturber
 * @param e: softening constant
 * @param x_dir: if True, calculates the acceleration in the x direction, otherwise calculates the acceleration in the y-direction
 * 
 * @return the acceleration of the star in either the x or y direction
 */
float accel(float x, float y,float x_c,float y_c, float x_p, float y_p, float M, float m_p, float e, bool x_dir){
    float D_C = sqrt(pow(x-x_c,2) + pow(y-y_c,2) + pow(e,2));
    float D_P = sqrt(pow(x-x_p,2) + pow(y-y_p,2) + pow(e,2));

    float acceleration;

    if (x_dir == true){
        acceleration = -G*M*(x-x_c)/pow(D_C,3) - G*m_p*(x-x_p)/pow(D_P,3);
    }
    else{
        acceleration = -G*M*(y-y_c)/pow(D_C,3) - G*m_p*(y-y_p)/pow(D_P,3);
    }
    return acceleration;
}

/**
 * 
 * @brief Uses the taylor expansion method to find the new position and speed of a star at each time-step
 * 
 * @param x: x coordinate of the star
 * @param y: y coordinate of the star
 * @param x_c: x coordinate of the galaxy centre of mass
 * @param y_c: y coordinate of the galaxy centre of mass
 * @param x_p: x coordinate of the point perturber
 * @param y_p: y coordinate of the point perturber
 * @param v_x: x component of the star velocity
 * @param v_y: y component of the star velocity
 * @param M: mass of the galaxy
 * @param m_p: mass of the point perturber
 * @param a: time-step
 * @param e: the softening constant
 * 
 * @return x_new, the updated x coordinate of the star. y_new (float), the updated y coordinate of the 
 * star. v_x_new, the updated x component of the star velocity and v_y_new (float), the updated y 
 * component of the star velocity
 */
std::vector<float> taylor_one_step(float x,float y, float x_c, float y_c, float x_p, float y_p, float v_x, float v_y, float M, float m_p, float a, float e){

    float x_accel = accel(x,y,x_c, y_c, x_p, y_p,M, m_p, e, true);
    float y_accel = accel(x,y,x_c, y_c, x_p, y_p,M, m_p, e, false);

    float x_new = x+a*v_x + (pow(a,2)/2)*x_accel;
    float y_new = y+a*v_y + (pow(a,2)/2)*y_accel;

    float v_x_new = v_x + a*x_accel; 
    float v_y_new = v_y + a*y_accel;

    std::vector<float> results = {x_new, y_new, v_x_new, v_y_new};
    return results;
}

/**
 * 
 * @brief Uses the leapfrog method of integration to find the new position and speed of a star at each time-step
 * 
 * @param x: x coordinate of the star
 * @param y: y coordinate of the star
 * @param x_c: x coordinate of the galaxy centre of mass
 * @param y_c: y coordinate of the galaxy centre of mass
 * @param x_p: x coordinate of the point perturber
 * @param y_p: y coordinate of the point perturber
 * @param v_x: x component of the star velocity
 * @param v_y: y component of the star velocity
 * @param M: mass of the galaxy
 * @param m_p: mass of the point perturber
 * @param a: time-step
 * @param e: the softening constant
 * 
 * @return x_new, the updated x coordinate of the star. y_new (float), the updated y coordinate of the 
 * star. v_x_new, the updated x component of the star velocity and v_y_new (float), the updated y 
 * component of the star velocity
 */
std::vector<float> leapfrog_one_step(float x,float y, float x_c, float y_c, float x_p, float y_p, float v_x, float v_y, float M, float m_p, float a, float e){

    float x_accel = accel(x,y,x_c, y_c, x_p, y_p,M, m_p, e, true);
    float y_accel = accel(x,y,x_c, y_c, x_p, y_p,M, m_p, e, false);

    float x_new = x+a*v_x + (pow(a,2)/2)*x_accel;
    float y_new = y+a*v_y + (pow(a,2)/2)*y_accel;

    float x_accel_new = accel(x_new,y_new,x_c, y_c, x_p, y_p,M, m_p, e, true);
    float y_accel_new = accel(x_new,y_new,x_c, y_c, x_p, y_p,M, m_p, e, false);

    float v_x_new = v_x + 0.5*a*(x_accel+x_accel_new); 
    float v_y_new = v_y + 0.5*a*(y_accel+y_accel_new);

    std::vector<float> results = {x_new, y_new, v_x_new, v_y_new};
    return results;
}

/**
 * 
 * @brief At each time-step calculates the new coordinates and velocities of each star using either the taylor expansion 
 * or leapfrog method, then recalculates the position of the barycentre of the galaxy. Assumes the perturber moves diagonally, 
 * feeling no force from the galaxy, and the galaxy is initially stationary.
 * 
 * @param x: starting x coordinates of the stars
 * @param y: starting y coordinates of the stars
 * @param v_x: starting x components of the stars's velocities
 * @param v_y: starting y components of the stars's velocities
 * @param v_p: velocity of the perturber in both the x and y direction
 * @param m: masses of the stars
 * @param m_p: the mass of the point perturber
 * @param a: time-step
 * @param e: the softening constant
 * @param n: the number of time-steps
 * @param initial: starting [x,y] position of the point perturber
 * @param leapfrog: if true uses the leapfrog method, else uses the Taylor expansion method
 * 
 * @return an tuple containing the positions of each star at each time-step, the position of the perturber
 * at each time-step, and the position of the barycentre of the galaxy at each time-step.
 */
std::tuple< 
std::vector<std::vector<std::vector<float>>>, std::vector<std::vector<float>>, std::vector<std::vector<float>>
> 
all_stars(std::vector<float> x,std::vector<float> y, std::vector<float> v_x, std::vector<float> v_y, float v_p, std::vector<float> m, float m_p, float a, float e, int n, std::vector<float> initial, bool leapfrog){
   
    int num_stars = x.size();
    //to store the data we create an array, with a row for each time step so [time step][stars][x and y].
    std::vector<std::vector<std::vector<float>>> positions_array(n, std::vector<std::vector<float>>(num_stars, std::vector<float>(2))
    );
    std::vector<std::vector<float>> perturber_positions;    //create vector to store perturber positions
    std::vector<std::vector<float>> barycentre;      //create vector to store galaxy centre of mass positions

    float M = std::accumulate(m.begin(), m.end(), 0.0f);   //find the mass of the galaxy

    float x_c = 0.0f;   //the galaxy is initially centered at 0.
    float y_c = 0.0f;
    
    for (int i = 0; i<n; i++){
        barycentre.push_back({x_c, y_c});

        //calculate perturber position when moving in a linear trajectory
        float x_p = initial[0] + v_p*i*a;
        float y_p = initial[1] + v_p*i*a;

        perturber_positions.push_back({x_p, y_p});

        for(int j = 0; j<num_stars; j++){
            
            std::vector<float> results;
            if (leapfrog == true){
                results = leapfrog_one_step(x[j],y[j],x_c, y_c, x_p, y_p, v_x[j], v_y[j], M, m_p, a, e);
            } else{
                results = taylor_one_step(x[j],y[j],x_c, y_c, x_p, y_p, v_x[j], v_y[j], M, m_p, a, e);
            }
        
            //update the list of star positions and velocities with the new values
            x[j] = results[0];
            y[j] = results[1];
            v_x[j] = results[2];
            v_y[j] = results[3];   
            
            positions_array[i][j][0] = results[0];
            positions_array[i][j][1] = results[1];
        }

    //update the galaxy centre of mass
    x_c = std::inner_product(x.begin(), x.end(), m.begin(), 0.0f)/M;  
    y_c = std::inner_product(y.begin(), y.end(), m.begin(), 0.0f)/M;
    }

    return std::make_tuple(positions_array, perturber_positions, barycentre);
}
            
/**
 * 
 * @brief creates an analogue of the numpy linspace function
 * 
 * @param start: starting value of the vector
 * @param end: end value of the vector
 * @param q: number of values in the vector
 * 
 * @return the vector result
 * 
 */
std::vector<float> linspace(float start, float end, int q) {
    std::vector<float> result;
    if (q == 0) return result;
    if (q == 1) {
        result.push_back(start);
        return result;
    }
    float step = (end - start) / (q - 1);
    for (int i = 0; i < q; i++) {
        result.push_back(start + step * i);
    }
    return result;
}

/**
 * 
 * @brief creates a vector of random numbers of length q, picked a uniformly distributed range of numbers
 * 
 * @param start: lower bound of range
 * @param end: upper bound of range
 * @param q: number of values in the vector
 * 
 * @return vector of random numbers
 */
std::vector<float> random_number(float start, float end, int q){
    std::vector<float> number_vector;
    std::random_device rd;
    std::uniform_real_distribution<float> dist(start, end);

    for (int i = 0; i < q; ++i) {
        number_vector.push_back(dist(rd));
    }
    return number_vector;
}

/**
 * @brief Given certain parameters of a galaxy, returns the initial conditions of the stars.
 * 
 * @param p: number of stars
 * @param r_max: radius of the galaxy
 * 
 * @return a tuple containing the initial x and y coordinates of the stars, the x and y components of the
 * initial velocities of the stars, the masses of the stars, the total mass of the galaxy, the softening 
 * constant, the shortest and longest orbital periods of stars in the galaxy, and the speed of the slowest
 * orbiting star in the galaxy
 */
std::tuple<
std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, std::vector<float>, float, float, float, float, float
>
initialise_values(int p, float r_max){

    // distribute stars according to the exponential distribution
    float R_d = r_max/3;    //set scale length of the galaxy to be a third of the radius
    std::vector<float> u = random_number(0,1,p);

    std::vector<float> r;  //find circular velocities of the stars
    for(int i = 0; i<p; i++){
        r.push_back(-R_d*std::log(1-u[i]*(1-std::exp(-r_max/R_d))));
    }

    std::vector<float> theta = random_number(0, 2*pi, p);
    std::vector<float> m(p,1);

    float r_max_val = *std::max_element(r.begin(), r.end());
    float r_min_val = *std::min_element(r.begin(), r.end());

    float e = r_max/sqrt(p);    //define the softening constant
    std::cout<< "Softening constant is " << e <<std::endl;

    float M = std::accumulate(m.begin(), m.end(), 0.0f);   //total mass of the galaxy

    std::vector<float> v;  //find circular velocities of the stars
    for(int i = 0; i<r.size(); i++){
        v.push_back(sqrt(abs(G*M/r[i])));
    }

    float v_max_val = *std::max_element(v.begin(), v.end());
    float v_min_val = *std::min_element(v.begin(), v.end());

    float T_max = 2*pi*r_max_val/v_min_val;
    float T_min = 2*pi*r_min_val/v_max_val;

    std::cout << "Shortest and longest orbital periods are " <<  T_min << " " << T_max << std::endl;

    //find the x and y components of the initial positions and velocities
    std::vector<float> x;
    std::vector<float> y;
    std::vector<float> v_x;
    std::vector<float> v_y;

    for(int i = 0; i<r.size(); i++){
        x.push_back(r[i]*cos(theta[i]));
        y.push_back(r[i]*sin(theta[i]));
        v_x.push_back(v[i]*sin(theta[i]));
        v_y.push_back(-v[i]*cos(theta[i]));
    }

    return std::make_tuple(x,y,v_x, v_y, m, M, e, T_min, T_max, v_min_val);
}

int main(){

    /* Here the ratios of the perturber mass and speed to the galaxy can be altered to see the effect on
    the galaxy. The initial position of the perturber can also be altered the change the distance of
    closest approach */

    int p = 1000;    //number of stars
    float mu = 5;    //mass ratio of perturber to galaxy
    float V = 1.2;   //ratio of perturber velocity to slowest star's circular velocity
    float r_max = 5;    //radius of galaxy
    int n = 50000;    //number of time steps
    std::vector<float> initial = {-20, 10};

    auto values = initialise_values(p, r_max);

    float v_p = V*std::get<9>(values);
    std::cout << "Perturber speed is " << v_p << std::endl;

    // set the time period of the simulation to be triple the time taken for the perturber to reach the galaxy 
    float T = 3*sqrt(pow(initial[0],2)+pow(initial[1],2))/v_p;
    std::cout << "Time scale of simulation is " << T << std::endl;

    float m_p = mu*std::get<5>(values);
    float a = T/n;
    std::cout << "Time-step is " << a <<std::endl;

    auto results = all_stars(std::get<0>(values),std::get<1>(values), std::get<2>(values), std::get<3>(values), v_p, std::get<4>(values), m_p, a, std::get<6>(values), n, initial, true);
    
    std::cout<< "Simulation run, saving data now" << std::endl;

    std::vector<std::vector<std::vector<float>>> stars = std::get<0>(results);
    std::vector<std::vector<float>> perturber = std::get<1>(results);   
    std::vector<std::vector<float>> barycentre = std::get<2>(results); 

    //save the stars, perturber, and barycentre data
    std::ofstream out("simp_ex2_stars.txt");
    for (int i = 0; i < stars.size(); i++){
        for (int j =0; j < stars[0].size(); j++){
            out << stars[i][j][0] << "\t" << stars[i][j][1] << "\n";
        }
    }

    std::ofstream out1("simp_ex2_perturber.txt");
    for (int i = 0; i < perturber.size(); ++i) {
        out1 << perturber[i][0] << "\t" << perturber[i][1] << "\n";
    }

    std::ofstream out2("simp_ex2_barycentre.txt");
    for (int i = 0; i < barycentre.size(); ++i) {
        out2 << barycentre[i][0] << "\t" << barycentre[i][1] << "\n";
    }

    //save the data containing the parameters
    std::ofstream out3("simp_ex2_parameters.txt");
        out3 << p << "\t" << mu << "\t" << V << "\t" << n ;

    out.close();
    out1.close();
    out2.close();
    out3.close();
}