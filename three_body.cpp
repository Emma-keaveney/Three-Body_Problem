#include <iostream>
#include<vector>
#include<cmath>
#include <fstream>   //for saving to computer

const double G =6.6726*pow(10,-11);    // gravitational constant
const double m_E = 5.9742*pow(10,24);     // mass of the earth
const double m_M = 7.35*pow(10,22);       // mass of the moon
const double d = 3.844*pow(10,8);         // distance from the moon to the earth
const double pi = 3.14159265358979323846;   //pi

const double r_E = d*m_M/(m_M+m_E);
const double r_M = d*m_E/(m_M+m_E);

const double T = sqrt(4*pow(pi,2)*pow(d,3)/(G*(m_M+m_E)));
const double a = 2; // set the time step
int n = int(T/a);  // number of time steps

// Gives the acceleration in the x and y directions for a particle given positions of moon and earth
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

// Uses the taylor approximation to find all the positions across n time steps of size a
auto taylor(double r,double v, int n, double a, std::vector<double> x_E, std::vector<double> y_E, std::vector<double> x_M, std::vector<double> y_M){
    double x = r;
    double y = 0;
    double v_x = 0;
    double v_y = v;

    std::vector<double> xs;
    std::vector<double> ys;

    for (int i = 0; i<n; i++){
        xs.push_back(x);
        ys.push_back(y);

        double x_new = x+a*v_x + (pow(a,2)/2)*accel(x,y,x_E[i],y_E[i],x_M[i],y_M[i], true);
        v_x = v_x + a*accel(x,y,x_E[i],y_E[i],x_M[i],y_M[i], true);
        
        double y_new = y+a*v_y + (pow(a,2)/2)*accel(x,y,x_E[i],y_E[i],x_M[i],y_M[i], false); 
        v_y = v_y + a*accel(x,y,x_E[i],y_E[i],x_M[i],y_M[i], false);

        x = x_new;
        y = y_new;
    }

    std::vector<std::vector<double>> x_and_y = {xs,ys};
    return x_and_y;
}

auto RK(double r,double v, int n, double a, std::vector<double> x_E, std::vector<double> y_E, std::vector<double> x_M, std::vector<double> y_M){
    double x = r;
    double y = 0;
    double v_x = 0;
    double v_y = v;

    std::vector<double> xs;
    std::vector<double> ys;

    for (int i = 0; i<n; i++){
        xs.push_back(x);
        ys.push_back(y); 

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
    std::vector<std::vector<double>> x_and_y = {xs,ys};
    return x_and_y;
}

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

auto earth_moon_positions(double r_E, double r_M, double T, int n){
    std::vector<double> x_E;
    std::vector<double> y_E;
    std::vector<double> x_M;
    std::vector<double> y_M;

    std::vector<double> times = linspace(0,T,n);

    for(int i = 0; i<n; i++){
        x_E.push_back(r_E*cos(2*pi/T*times[i]));
        y_E.push_back(r_E*sin(2*pi/T*times[i]));
        x_M.push_back(r_M*cos(2*pi/T*times[i]));
        y_M.push_back(r_M*sin(2*pi/T*times[i]));
    }

    std::vector<std::vector<double>> moon_earth = {x_M,y_M, x_E, y_E};
    return moon_earth;
}  

double L_2 = d * pow(m_M / (3 * m_E), double(1.0/3.0));      //starting distance
double r = L_2*1.075377 + r_M;
double v = 2*pi*r/T;                      //starting speed

double r_RK = L_2*1.075383 + r_M;
double v_RK = 2*pi*r_RK/T;

int main(){

    std::cout << "Period of orbit is " << T/3600/24 <<" days"<<std::endl;
    std::cout << "Step size is " << a << " s"<< std::endl;

    std::vector<std::vector<double>> moon_earth = earth_moon_positions(r_E, r_M, T, n);
    std::vector<std::vector<double>> satellite_taylor = taylor(r, v, n+1, a, moon_earth[2], moon_earth[3], moon_earth[0], moon_earth[1]);
    std::cout << "Taylor run" << std::endl;

    std::vector<std::vector<double>> satellite_RK = RK(r_RK, v_RK, n+1, a, moon_earth[2], moon_earth[3], moon_earth[0], moon_earth[1]);
    std::cout << "Runge Kutta run" << std::endl;

    std::ofstream out("satellite_taylor.txt");
    for (int i = 0; i < satellite_taylor[0].size(); ++i) {
        out << satellite_taylor[0][i] << "\t" << satellite_taylor[1][i] << "\n";
    }
    std::ofstream out1("moon.txt");
    for (int i = 0; i < moon_earth[0].size(); ++i) {
       out1 << moon_earth[0][i] << "\t" << moon_earth[1][i] << "\n";
    }
    std::ofstream out2("earth.txt");
    for (int i = 0; i < moon_earth[0].size(); ++i) {
       out2 << moon_earth[2][i] << "\t" << moon_earth[3][i] << "\n";
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