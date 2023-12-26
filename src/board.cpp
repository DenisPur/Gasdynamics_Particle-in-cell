#pragma once
#include <iostream>
#include <cmath>
#include "array.cpp"
#include "particles.cpp"

class Board_2IdealGases {
public:
    Board_2IdealGases(int n, double cell_size) : Board_2IdealGases(n, n, cell_size) {}

    Board_2IdealGases(int ny, int nx, double cell_size) {
        _ny = ny;
        _nx = nx;
        _size = cell_size;

        _massesA = Array(_ny, _nx);
        _massesB = Array(_ny, _nx);
        _energiesA = Array(_ny, _nx);
        _energiesB = Array(_ny, _nx);

        _pressures = Array(_ny, _nx);
        _presures_shifted_x = Array(_ny, _nx + 1);
        _presures_shifted_y = Array(_ny + 1, _nx);

        _vx = Array(_ny, _nx); 
        _vy = Array(_ny, _nx);

        _vx_tilted = Array(_ny, _nx); 
        _vy_tilted = Array(_ny, _nx);

        _w_energies = Array(_ny, _nx); 
        _w_energies_tilted = Array(_ny, _nx);

        _vx_shifted = Array(_ny, _nx + 1);
        _vy_shifted = Array(_ny + 1, _nx);
    }

    //------------------------------------------------------------------------

    void add_particlesA(int number, double mass) {
        _particlesA = Particles(number);
        _particlesA.set_mass_for_each(mass);
        _particlesA.evenly_distribute(_size * _nx, _size * _ny);
    }
    
    void add_particlesB(int number, double mass) {
        _particlesB = Particles(number);
        _particlesB.set_mass_for_each(mass);
        _particlesB.evenly_distribute(_size * _nx, _size * _ny);
    }

    //------------------------------------------------------------------------

    void initiate_energy_test_function() {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                _energiesA(y, x) = x * y + 1;
                _energiesB(y, x) = (_nx - x) * (_ny - y) + 1;
            }
        }
    }

    //------------------------------------------------------------------------

    void re_mass() {
        _massesA.zeros();
        for (int i = 0; i < _particlesA.len(); i++) {
            int nx = std::floor(_particlesA(i, 0) / _size);
            int ny = std::floor(_particlesA(i, 1) / _size);
            _massesA(ny, nx) += _particlesA(i, 2);
        }

        _massesB.zeros();
        for (int i = 0; i < _particlesB.len(); i++) {
            int nx = std::floor(_particlesB(i, 0) / _size);
            int ny = std::floor(_particlesB(i, 1) / _size);
            _massesB(ny, nx) += _particlesB(i, 2);
        }
    }

    void re_pressure() {
        double eps = 1e-9;

        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                double sA = _massesA(y, x) * _energiesA(y, x) + eps;
                double sB = _massesB(y, x) * _energiesB(y, x) + eps;
                
                double sigmaA = sA / (sA + sB);
                _pressures(y, x) = sA / (sigmaA * _size * _size);
            }
        }

        for (int y = 0; y < _ny; y++) {
            _presures_shifted_x(y, 0) = _pressures(y, 0);
            for (int x = 1; x < _nx; x++) {
                _presures_shifted_x(y, x) = (_pressures(y, x - 1) + _pressures(y, x)) / 2;
            }
            _presures_shifted_x(y, _nx) = _pressures(y, _nx - 1);
        }

        for (int x = 0; x < _nx; x++) {
            _presures_shifted_y(0, x) = _pressures(0, x);
        }
        for (int y = 1; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                _presures_shifted_y(y, x) = (_pressures(y - 1, x) + _pressures(y, x)) / 2;
            }
        }
        for (int x = 0; x < _nx; x++) {
            _presures_shifted_y(_ny, x) = _pressures(_ny - 1, x);
        }
    }

    void re_v_tilted(double tau) {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                double rho = (_massesA(y, x) + _massesB(y, x)) / (_size * _size);
                _vx_tilted(y, x) = _vx(y, x) + tau / rho / _size * (_presures_shifted_x(y, x) - _presures_shifted_x(y, x + 1));
                _vy_tilted(y, x) = _vy(y, x) + tau / rho / _size * (_presures_shifted_y(y, x) - _presures_shifted_y(y + 1, x));
            }
        }
    }

    void re_w_energies() {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                double potential = _massesA(y, x) * _energiesA(y, x) + _massesB(y, x) * _energiesB(y,x);
                double kinetic = (_massesA(y, x) + _massesB(y, x)) * (_vx(y, x)*_vx(y, x) + _vy(y,x)*_vy(y,x)) / 2;
                _w_energies(y, x) = (potential + kinetic) / (_size*_size);
            }
        }
    }

    void re_v_shifted() {
        _vx_shifted.zeros();
        for (int y = 0; y < _ny; y++) {
            // _vx_shifted(y, 0) = (_vx(y, 0) + _vx_tilted(y, 0)) / 2;
            _vx_shifted(y, 0) = 0;
            for (int x = 1; x < _nx; x++) {
                _vx_shifted(y, x) = (_vx(y, x - 1) + _vx(y, x) + _vx_tilted(y, x - 1) + _vx_tilted(y, x)) / 4;
            }
            // _vx_shifted(y, _nx) = (_vx(y, _nx - 1) + _vx_tilted(y, _nx - 1)) / 2;
            _vx_shifted(y, _nx) = 0;
        }

        _vy_shifted.zeros();
        for (int x = 0; x < _nx; x++) {
            // _vy_shifted(0, x) = (_vy(0, x) + _vy_tilted(0, x)) / 2;
            _vy_shifted(0, x) = 0;
        }
        for (int y = 1; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                _vy_shifted(y, x) = (_vy(y - 1, x) + _vy(y, x) + _vy_tilted(y - 1, x) + _vy_tilted(y, x)) / 4;
            }
        }
        for (int x = 0; x < _nx; x++) {
            // _vy_shifted(_ny, x) = (_vy(_ny - 1, x) + _vy_tilted(_ny - 1, x)) / 2;
            _vy_shifted(_ny, x) = 0;
        }
    }

    void re_w_energies_tilted(double tau) {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                _w_energies_tilted(y, x) = _w_energies(y, x) + tau / _size * ( \
                    _vx_shifted(y, x) * _presures_shifted_x(y, x) + \
                    _vy_shifted(y, x) * _presures_shifted_y(y, x) - \
                    _vx_shifted(y, x + 1) * _presures_shifted_x(y, x + 1) - \
                    _vy_shifted(y + 1, x) * _presures_shifted_y(y + 1, x));
            }
        }
    }

    //------------------------------------------------------------------------

    void print_masses() {
        std::cout << "# masses A:\n";
        _massesA.print();
        std::cout << "# masses B:\n";
        _massesB.print();
    }

    void print_pressures(bool shifted = false) {
        std::cout << "# pressure :\n";
        _pressures.print();
        if (shifted) {
            std::cout << "# shifted x :\n";
            _presures_shifted_x.print();
            std::cout << "# shifted y\n";
            _presures_shifted_y.print();
        }
    }

    void print_v_tilted() {
        std::cout << "vx tilted :\n";
        _vx_tilted.print();
        std::cout << "vy tilted :\n";
        _vy_tilted.print();
    }

    void print_v(bool shifted = false) {
        std::cout << "# vx :\n";
        _vx.print();
        std::cout << "# vy :\n";
        _vy.print();
        if (shifted) {
            std::cout << "# shifted x :\n";
            _vx_shifted.print();
            std::cout << "# shifted y :\n";
            _vy_shifted.print();
        }
    }

    void print_energies_tilted() {
        std::cout << "# energies :\n";
        _w_energies.print();
        std::cout << "# energies tilted :\n";
        _w_energies_tilted.print();
    }

private:
    Particles _particlesA;
    Array _massesA;
    Array _energiesA;

    Particles _particlesB;
    Array _massesB;
    Array _energiesB;

    Array _pressures;
    Array _vx;
    Array _vy;

    Array _presures_shifted_x;
    Array _presures_shifted_y;

    Array _vx_tilted;
    Array _vy_tilted;

    Array _w_energies;
    Array _w_energies_tilted;

    Array _vx_shifted;
    Array _vy_shifted;

    // Array _energiesA_tilted;
    // Array _energiesB_tilted;

    double _size;
    int _nx;
    int _ny;
};