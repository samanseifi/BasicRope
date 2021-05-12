/* Simulating a hanging rope:
 *      L:  initial length of rope
 *      N:  number of mass points
 */

#include <iostream>
#include <vector>
#include "Vec.h"

using namespace std;

const double PI = 3.1415926535;     // pi number
const double G = -9.8;              // earth's gravity $kg/(m sec^2)$
const int N = 10;                   // Number of mass points
const double L = 10;                // Length of rope

/* Calculates the internal forces due to contraction and expansion of the springs */
vector<Vec> internal_forces(vector<Vec> pos, double init_length, double k);

/* Calculates the viscous forces to induce damping */
vector<Vec> viscous_forces(vector<Vec> vel, double c);

/* Calculates the gravity forces on each mass */
vector<Vec> gravity_forces(double mass,  const Vec& g);

/* Print out the vector values */
void print_vector(const vector<Vec>& v) {
    for (auto & vec_i : v) {
        cout << vec_i.x << " " << vec_i.y << endl;
    }
}

int main() {

    const double k = 1.0;             // Spring stiffness [N/m]
    const double m = 0.01;            // Mass [kg]
    const double c = 0.001;           // Viscous damping coeff. [N sec/m]
    const Vec g(0.0, G);              // Gravity acceleration vector $\mathbf{g} = G\,\hat{\mathbf{j}}$
    const double dt = 0.001;          // time step [sec]

    const double t_final = 300;       // final time of simulation [sec]

    // Check the numerical instability for Euler-forward explicit scheme
    try {
        if (dt < PI * sqrt(m / k)) {
            cout << "PASS: Stable numerical scheme." << endl;
        } else throw (505);
    }
    catch (int err) {
        cout << "ERROR: dt should be less than " << PI * sqrt(m / k) << endl;
    }

    /* Initialize the position vector with N mass points equally distanced with init_ID parameter: */
    //  init_DEG = 0             :   Slope of the initial configuration in DEG
    float init_DEG = PI/4.0;

    // Create the vector of positions
    vector<Vec> pos(N);
    double init_length = L / (N - 1);
    for (int i = 0; i < N; i++) {
        pos[i] = Vec(cos(init_DEG) * init_length * i, -sin(init_DEG) * init_length * i);
    }

    // Declaration and initialization the velocity and acceleration vectors (STL) of vectors (2D-Cartesians)
    vector<Vec> disp(N, Vec(0,0));
    vector<Vec> vel(N, Vec(0,0));
    vector<Vec> acc(N, Vec(0,0));

    // Vector of forces
    vector<Vec> forces(N, Vec(0,0));

    // Vector of applied velocity
    vector<Vec> applied_vel(N, Vec(0.01, 0.01));

    // time zero!
    double t = 0.0;

    // Time stepping!
    while (t < t_final) {

        /* loop through all of the mass points */
        for (int i = 1; i < N; i++) {
            /* updating the forces */
            forces[i] = internal_forces(pos, init_length, k)[i] + viscous_forces(vel, c)[i] + gravity_forces(m, g)[i];

            /* calculate accelerations at time t+dt: a(t + dt) = 1/m * F(t) */
            acc[i] = forces[i].setToScale(1.0/m);

            /* calculate velocities at t+dt: v(t + dt) = v(t) + dt * a(t + dt) */
            if (t > 200) {
                vel[i] = vel[i] + applied_vel[i] + acc[i].setToScale(dt);
            } else {
                vel[i] = vel[i] + acc[i].setToScale(dt);
            }

            /* update the positions: p(t + dt) = p(t) + disp(t + dt) with disp(t) = dt * v(t + dt) */
            disp[i] = vel[i].setToScale(dt);
            pos[i] = pos[i] + disp[i];

        }

        /* Keep the fixed boundary conditions */
        acc[0] = acc[N-1] = Vec(0.0, 0.0);
        vel[0] = vel[N-1] = Vec(0.0, 0.0);
        pos[0] = Vec(0.0, 0.0);
        pos[N-1] = Vec(cos(init_DEG) * init_length * (N-1), -sin(init_DEG) * init_length * (N-1));
        
        /* Update time */
        t += dt;
    }

    /* print out final positions */
    print_vector(pos);

    return 0;
}

/*  This function calculate all the internal forces due to spring contraction and expansion  */
vector<Vec> internal_forces(vector<Vec> pos, const double init_length, const double k) {

    vector<Vec> F_int(N, Vec(0.0, 0.0));

    // Structural Springs
    Vec f0 = (pos[1] - pos[0]).UnitVec().setToScale(k*((pos[1] - pos[0]).Magnitude() - init_length));
    F_int[0] = (f0);

    for (int i = 1; i < N - 1; i++) {
        Vec f_i1 = (pos[i-1] - pos[i]).UnitVec().setToScale(k*((pos[i-1] - pos[i]).Magnitude() - init_length));
        Vec f_i2 = (pos[i+1] - pos[i]).UnitVec().setToScale(k*((pos[i+1] - pos[i]).Magnitude() - init_length));
        F_int[i] = (f_i1 + f_i2);
    }

    Vec f_N = (pos[N-2] - pos[N-1]).UnitVec().setToScale(k*((pos[N-2] - pos[N-1]).Magnitude() - init_length));
    F_int[N-1] = (f_N);

    // Bending Springs (the initial length is twice the distance between two nodes at rest)
    if (N > 2) {
        Vec f_b0 = (pos[2] - pos[0]).UnitVec().setToScale(k*((pos[2] - pos[0]).Magnitude() - 2*init_length));
        F_int[0] += f_b0;
        Vec f_b1 = (pos[3] - pos[1]).UnitVec().setToScale(k*((pos[3] - pos[1]).Magnitude() - 2*init_length));
        F_int[1] += f_b1;

        for (int i = 2; i < N - 2; i++) {
            Vec f_bi1 = (pos[i-2] - pos[i]).UnitVec().setToScale(k*((pos[i-2] - pos[i]).Magnitude() - 2*init_length));
            Vec f_bi2 = (pos[i+2] - pos[i]).UnitVec().setToScale(k*((pos[i+2] - pos[i]).Magnitude() - 2*init_length));
            F_int[i] += f_bi1;
            F_int[i] += f_bi2;
        }

        Vec f_bN_1 = (pos[N-4] - pos[N-2]).UnitVec().setToScale(k*((pos[N-4] - pos[N-2]).Magnitude() - 2*init_length));
        F_int[N-2] += f_bN_1;
        Vec f_bN = (pos[N-2] - pos[N-1]).UnitVec().setToScale(k*((pos[N-2] - pos[N-1]).Magnitude() - 2*init_length));
        F_int[N-1] += f_bN;
    }

    return F_int;
}

/* calculate the viscous forces */
vector<Vec> viscous_forces(vector<Vec> vel, const double c) {

    vector<Vec> F_v(N, Vec(0.0, 0.0));

    for (int i = 0; i < N; i++) {
        F_v[i] = vel[i].setToScale(-c);
    }

    return F_v;
}

/* calculate the gravity forces */
vector<Vec> gravity_forces(double m, const Vec& g) {
    vector<Vec> F_g(N, Vec(0.0, 0.0));

    for (int i = 0; i < N; i++) {
        F_g[i] = g.setToScale(m);
    }

    return F_g;
}