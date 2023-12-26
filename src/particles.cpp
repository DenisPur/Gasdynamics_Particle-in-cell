#pragma once
#include "array.cpp"

class Particles {
public:
    Particles() {
        _n = 0;
        _data = Array();
    }

    Particles(int n) {
        _n = n;
        _data = Array(n, 3);
    }

    int len() {
        return _n;
    }

    double& operator()(int i, int p) {
        return _data(i, p);
    }

    void set_mass_for_each(double m) {
        for (int i = 0; i < _n; i++) {
            _data(i, 2) = m;
        }
    }

    void evenly_distribute(double x_max, double y_max) {
        for (int i = 0; i < _n; i++) {
            _data(i, 0) = x_max * rand() / RAND_MAX;
            _data(i, 1) = y_max * rand() / RAND_MAX;
        }
    }

    void move_particles(Array vx, Array vy) {
        //
    }

    Array& get_data() {
        return _data; 
    }

private:
    int _n;
    Array _data; // [ [x0, y0, m0], [x1, y1, m1], ... ]
};