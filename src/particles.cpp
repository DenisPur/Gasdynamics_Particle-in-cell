#pragma once
#include <cmath>
#include "array.cpp"

class Particles {
public:
    Particles() {
        _n = 0;
        _data = Array();
    }

    Particles(int n) {
        _n = n;
        _data = Array(n, 6);
    }

    double& operator()(int number, int property) {
        return _data(number, property);
    }

    void set_mass_for_each(double m) {
        for (int i = 0; i < _n; i++) {
            _data(i, 2) = m;
        }
    }

    void set_borders(double x_max, double y_max, double size) {
        _x_max = x_max;
        _y_max = y_max;
        _size = size;
    }

    void evenly_distribute() {
        for (int i = 0; i < _n; i++) {
            _data(i, 0) = _x_max * rand() / RAND_MAX;
            _data(i, 1) = _y_max * rand() / RAND_MAX;
        }
    }

    void set_energies(Array inner_energy_array, Array vx_tilda_array, Array vy_tilda_array) {
        for (int i = 0; i < _n; i++) {
            int nx = std::floor(_data(i, 0) / _size);
            int ny = std::floor(_data(i, 1) / _size);
            double vx = vx_tilda_array(ny, nx);
            double vy = vy_tilda_array(ny, nx);
            _data(i, 5) = _data(i, 2) * (inner_energy_array(ny, nx) + (vx*vx + vy*vy) / 2);
        }
    }

    void move_particles(Array vx_array, Array vy_array, Array vx_tilda_array, Array vy_tilda_array, double tau) {
        for (int i = 0; i < _n; i++) {
            double x = _data(i, 0);
            double y = _data(i, 1);
            int nx = std::floor(x / _size);
            int ny = std::floor(y / _size);
            double vx = vx_array(ny, nx);
            double vy = vy_array(ny, nx);
            double vx_tilda = vx_tilda_array(ny, nx);
            double vy_tilda = vy_tilda_array(ny, nx);

            x += tau * vx;
            y += tau * vy;

            if (x < 0) {
                x = - x;
                vx = - vx;
            } else if (x > _x_max) {
                x = 2 * _x_max - x;
                vx = - vx;
            }
            
            if (y < 0) {
                y = - y;
                vy = - vy;
            } else if (y > _y_max) {
                y = 2 * _y_max - y;
                vy = - vy;
            }

            _data(i, 0) = x;
            _data(i, 1) = y;
            _data(i, 3) = vx_tilda;
            _data(i, 4) = vy_tilda;
        }
    }

    int len() {
        return _n;
    }

    Array& get_data() {
        return _data; 
    }

private:
    int _n;
    Array _data; // [ [x0, y0, m0, vx, vy, eps], [x1, y1, m1, vx, vy, eps], ... ]
    double _x_max;
    double _y_max;
    double _size;
};