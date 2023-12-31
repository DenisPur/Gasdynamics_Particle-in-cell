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

    void evenly_distribute(int i_from = -1, int i_to = -1) {
        if (i_from == -1) {
            for (int i = 0; i < _n; i++) {
                _data(i, 0) = _x_max * rand() / RAND_MAX;
                _data(i, 1) = _y_max * rand() / RAND_MAX;
            }
        } else {
            for (int i = i_from; i < i_to; i++) {
                _data(i, 0) = _x_max * rand() / RAND_MAX;
                _data(i, 1) = _y_max * rand() / RAND_MAX;
            }        
        }
    }

    void add_random_movements(double v_max, int i_from = -1, int i_to = -1) {
        if (i_from == -1) {
            for (int i = 0; i < _n; i++) {
                _data(i, 3) += v_max * 2 * (rand() / RAND_MAX - 0.5);
                _data(i, 4) += v_max * 2 * (rand() / RAND_MAX - 0.5);
            }
        } else {
            for (int i = i_from; i < i_to; i++) {
                _data(i, 3) += v_max * 2 * (rand() / RAND_MAX - 0.5);
                _data(i, 4) += v_max * 2 * (rand() / RAND_MAX - 0.5);
            }
        }
    }

    void set_coordinates(int i, double x, double y) {
        _data(i, 0) = x;
        _data(i, 1) = y;
    }

    void set_velocity(int i, double vx, double vy) {
        _data(i, 3) = vx;
        _data(i, 4) = vy;
    }

    //------------------------------------------------------------------------

    void set_energies_and_move(Array vx_array, Array vy_array, 
                               Array vx_tilda_array, Array vy_tilda_array,
                               Array inner_energy_array, double tau) {
        for (int i = 0; i < _n; i++) {
            double x = _data(i, 0);
            double y = _data(i, 1);
            int nx = std::floor(x / _size);
            int ny = std::floor(y / _size);
            double vx = vx_array(ny, nx);
            double vy = vy_array(ny, nx);
            double vx_tilda = vx_tilda_array(ny, nx);
            double vy_tilda = vy_tilda_array(ny, nx);

            _data(i, 5) = _data(i, 2) * (inner_energy_array(ny, nx) + (vx_tilda*vx_tilda + vy_tilda*vy_tilda) / 2);

            x += tau * vx;
            y += tau * vy;

            if (x < 0) {
                x = - x;
                vx = - vx;
            } else if (x >= _x_max) {
                x = 2 * _x_max - x;
                vx = - vx;
                if (x == _x_max) {
                    x = x - 0.001;
                }
            }

            if (y < 0) {
                y = - y;
                vy = - vy;
            } else if (y >= _y_max) {
                y = 2 * _y_max - y;
                vy = - vy;
                if (y == _y_max) {
                    y = y - 0.001;
                }
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