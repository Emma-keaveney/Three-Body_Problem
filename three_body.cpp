#include <iostream>
#include<vector>
#include<cmath>
#include <fstream>   //for saving to computer

const double G =6.6726*pow(10,-11);    // gravitational constant
const double m_E = 5.9742*pow(10,24);     // mass of the earth
const double m_M = 7.35*pow(10,22);       // mass of the moon
const double d = 3.844*pow(10,8);         // distance from the moon to the earth
const double pi = 3.14159265358979323846;   //pi

const double r_E = d*m_M/(m_M+m_E);    //distance of the Earth from the Earth-Moon barycentre
const double r_M = d*m_E/(m_M+m_E);    // distance of the Moon from the Earth-Moon barycentre

const double T = sqrt(4*pow(pi,2)*pow(d,3)/(G*(m_M+m_E)));    //time period of Earth-Moon orbit, found using Kepler's law
const double a = 2; // set the time step

/**
 * @brief  Gives the acceleration in the x or y direction for a satellite orbiting the Earth-Moon system
 * 
 * @param  x:  x position of the satellite
 * @param  y:  y position of the satellite
 * @param x_E: x position of the Earth
 * @param y_E: y position of the Earth
 * @param x_M: x position of the Moon
 * @param y_M: y position of the Moon
 * @param x_dir: if True, calculates the acceleration in the x direction, otherwise calculates the acceleration in the y-direction
 * 
 * @return the acceleration of the satellite in either the x or y direction
 */
double accel(double x, double y,double x_E,double y_E,double x_M,double y_M,bool x_dir){
    double D_E = sqrt(pow(x-x_E,2) + pow(y-y_E,2));
    double D_M = sqrt(pow(x-x_M,2) + pow(y-y_M,2));

    double acceleration;

    if (x_dir == true){
        acceleration = -G*m_E*(x-x_E)/pow(D_E,3) - G*m_M*(x-x_M)/pow(D_M,3);
    }
    else{
        acceleration = -G*m_E*(y-y_E)/pow(D_E,3) - G*m_M*(y-y_M)/pow(D_M,3);
    }
    return acceleration;
}

/**
 * @brief  an analogue to the numpy linspace function
 * 
 * Creates a vector of doubles of length n, with defined starting and ending points
 * 
 * @param  start:  starting value of the vector
 * @param  end:  end value of the vector
 * @param n: the length of the vector
 * 
 * @return the vector defined
 */
std::vector<double> linspace(double start, double end, int n) {
    std::vector<double> result;
    if (n == 0) return result;
    if (n == 1) {
        result.push_back(start);
        return result;
    }
    double step = (end - start) / (n - 1);
    for (int i = 0; i < n; i++) {
        result.push_back(start + step * i);
    }
    return result;
}

/**
 * @brief Uses the taylor approximation to find positions of a satellite orbiting the Earth-Moon barycentre
 * 
 * Calculates the cartesian positions of the Earth and Moon over a given time period. Then uses the taylor
 * approximation to find the position of a satellite orbiting the system.
 * 
 * @param  r:  the initial distance of the satellite to the Earth-Moon barycentre
 * @param  v:  the initial speed of the satellite
 * @param a: the size of the time-step
 * @param r_E: the distance from the Earth to the Earth-Moon barycenre
 * @param r_M: the distance from the Moon to the Earth-Moon barycenre
 * 
 * @return positions, a vector of doubles containing the cartesian coordinates of the satellite, Earth and Moon at each time-step
 */
auto taylor(double r,double v, double a, double T, double r_E, double r_M){
    
    // initialise the starting position and speed of the satellite
    double x = r;
    double y = 0;
    double v_x = 0;
    double v_y = v;

    // calculate the number of time-steps required, and time at each step
    int n = int(T/a)+1;
    std::vector<double> times = linspace(0,T,n);

    //create arrays to store the positions of the satellite, Earth, and Moon
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> x_E;
    std::vector<double> y_E;
    std::vector<double> x_M;
    std::vector<double> y_M;

    for (int i = 0; i<n; i++){

        xs.push_back(x);
        ys.push_back(y);

        // calculate cartesian coordinates of the Earth and Moon using circular motion equations
        x_E.push_back(r_E*cos(2*pi/T*times[i]));
        y_E.push_back(r_E*sin(2*pi/T*times[i]));
        x_M.push_back(r_M*cos(2*pi/T*times[i]));
        y_M.push_back(r_M*sin(2*pi/T*times[i]));

        double x_new = x+a*v_x + (pow(a,2)/2)*accel(x,y,x_E[i],y_E[i],x_M[i],y_M[i], true);
        v_x = v_x + a*accel(x,y,x_E[i],y_E[i],x_M[i],y_M[i], true);
        
        double y_new = y+a*v_y + (pow(a,2)/2)*accel(x,y,x_E[i],y_E[i],x_M[i],y_M[i], false); 
        v_y = v_y + a*accel(x,y,x_E[i],y_E[i],x_M[i],y_M[i], false);

        x = x_new;
        y = y_new;
    }

    std::vector<std::vector<double>> positions = {xs,ys, x_E, y_E, x_M, y_M};
    return positions;
}

/**
 * @brief Uses the Runge Kutta approximation to find positions of a satellite orbiting the Earth-Moon barycentre
 * 
 * Calculates the cartesian positions of the Earth and Moon over a given time period. Then uses the fourth order Runge-Kutta
 * approximation to find the position of a satellite orbiting the system
 * 
 * @param  r:  the initial distance of the satellite to the Earth-Moon barycentre
 * @param  v:  the initial speed of the satellite
 * @param a: the size of the time-step
 * @param r_E: the distance from the Earth to the Earth-Moon barycenre
 * @param r_M: the distance from the Moon to the Earth-Moon barycenre
 * 
 * @return positions, a vector of doubles containing the cartesian coordinates of the satellite, Earth and Moon at each time-step
 */
auto RK(double r,double v, double a, double T, double r_E, double r_M){

    // initialise the starting position and speed of the satellite
    double x = r;
    double y = 0;
    double v_x = 0;
    double v_y = v;

    // calculate the number of time-steps required, and time at each step
    int n = int(T/a)+1;
    std::vector<double> times = linspace(0,T,n);

    //create arrays to store the positions of the satellite, Earth, and Moon
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> x_E;
    std::vector<double> y_E;
    std::vector<double> x_M;
    std::vector<double> y_M;

    for (int i = 0; i<n; i++){
        xs.push_back(x);
        ys.push_back(y); 

        // calculate cartesian coordinates of the Earth and Moon using circular motion equations
        x_E.push_back(r_E*cos(2*pi/T*times[i]));
        y_E.push_back(r_E*sin(2*pi/T*times[i]));
        x_M.push_back(r_M*cos(2*pi/T*times[i]));
        y_M.push_back(r_M*sin(2*pi/T*times[i]));


        double accel_x = accel(x,y,x_E[i],y_E[i],x_M[i],y_M[i], true);
        double accel_y = accel(x,y,x_E[i],y_E[i],x_M[i],y_M[i], false);
        
        double z_x_1 = x + (a/2)*v_x;
        double z_vx_1 = v_x + (a/2)*accel_x;
        double z_y_1 = y + (a/2)*v_y;
        double z_vy_1 = v_y + (a/2)*accel_y;
        
        double z_accx_1 = accel(z_x_1,z_y_1,x_E[i],y_E[i],x_M[i],y_M[i], true);
        double z_accy_1 = accel(z_x_1,z_y_1,x_E[i],y_E[i],x_M[i],y_M[i], false);
        
        double z_x_2 = x + (a/2)*z_vx_1;
        double z_vx_2 = v_x + (a/2)*z_accx_1;
        double z_y_2 = y + (a/2)*z_vy_1;
        double z_vy_2 = v_y + (a/2)*z_accy_1;
        
        double z_accx_2 = accel(z_x_2,z_y_2,x_E[i],y_E[i],x_M[i],y_M[i], true);
        double z_accy_2 = accel(z_x_2,z_y_2,x_E[i],y_E[i],x_M[i],y_M[i], false);
        
        double z_x_3 = x + a*z_vx_2;
        double z_vx_3 = v_x + a*z_accx_2;
        double z_y_3 = y + a*z_vy_2;
        double z_vy_3 = v_y + a*z_accy_2;
        
        double z_accx_3 = accel(z_x_3,z_y_3,x_E[i],y_E[i],x_M[i],y_M[i], true);
        double z_accy_3 = accel(z_x_3,z_y_3,x_E[i],y_E[i],x_M[i],y_M[i], false);
        
        x = x + (a/6)*(v_x+2*z_vx_1+2*z_vx_2+z_vx_3);
        v_x = v_x + (a/6)*(accel_x + 2*z_accx_1 + 2*z_accx_2 + z_accx_3);
        y = y + (a/6)*(v_y+2*z_vy_1+2*z_vy_2+z_vy_3);
        v_y = v_y + (a/6)*(accel_y + 2*z_accy_1 + 2*z_accy_2 + z_accy_3);
    }
    std::vector<std::vector<double>> positions = {xs,ys, x_E, y_E, x_M, y_M};
    return positions;
}


double L_2 = d * pow(m_M / (3 * m_E), double(1.0/3.0)); 
double r = L_2*1.075377 + r_M;  //starting position of the satellite for Taylor known to result in stable orbit
double v = 2*pi*r/T;   //starting speed of the satellite found using circular motion

double r_RK = L_2*1.075383 + r_M;    //starting position of the satellite for RK4 known to result in stable orbit
double v_RK = 2*pi*r_RK/T;     // starting speed of the satellite found using circular motion

int main(){

    std::cout << "Period of orbit is " << T/3600/24 <<" days"<<std::endl;
    std::cout << "Step size is " << a << " s"<< std::endl;

    std::vector<std::vector<double>> satellite_taylor = taylor(r, v, a, T, r_E, r_M);  //run Taylor
    std::cout << "Taylor run" << std::endl;

    std::vector<std::vector<double>> satellite_RK = RK(r_RK, v_RK,a, T, r_E, r_M);   //run RK
    std::cout << "Runge Kutta run" << std::endl;

    //save the data
    std::ofstream out("satellite_taylor.txt");
    for (int i = 0; i < satellite_taylor[0].size(); ++i) {
        out << satellite_taylor[0][i] << "\t" << satellite_taylor[1][i] << "\n";
    }

    std::ofstream out1("earth.txt");
    for (int i = 0; i < satellite_taylor[0].size(); ++i) {
        out1 << satellite_taylor[2][i] << "\t" << satellite_taylor[3][i] << "\n";
    }

    std::ofstream out2("moon.txt");
    for (int i = 0; i < satellite_taylor[0].size(); ++i) {
        out2 << satellite_taylor[4][i] << "\t" << satellite_taylor[5][i] << "\n";
    }

    std::ofstream out3("satellite_RK.txt");
    for (int i = 0; i < satellite_RK[0].size(); ++i) {
       out3 << satellite_RK[0][i] << "\t" << satellite_RK[1][i] << "\n";
    }

out.close();
out1.close();
out2.close();
out3.close();

}
