#pragma once
#include <iostream>
#include <cmath>
#include <string>
#include "array.cpp"
#include "particles.cpp"

class Board_2IdealGases {
public:
    Board_2IdealGases(int n, double cell_size) : Board_2IdealGases(n, n, cell_size) {}

    Board_2IdealGases(int ny, int nx, double cell_size) {
        _ny = ny;
        _nx = nx;
        _size = cell_size;

        _mass_A = Array(_ny, _nx);
        _energy_A = Array(_ny, _nx);
        _energy_A_tilda = Array(_nx, _ny); 

        _mass_B = Array(_ny, _nx);
        _energy_B = Array(_ny, _nx);
        _energy_B_tilda = Array(_nx, _ny); 

        _pressure = Array(_ny, _nx);
        _pressure_shifted_x = Array(_ny, _nx + 1);
        _pressure_shifted_y = Array(_ny + 1, _nx);

        _w_energy = Array(_ny, _nx); 
        _w_energy_tilda = Array(_ny, _nx);

        _vx = Array(_ny, _nx); 
        _vy = Array(_ny, _nx);
        _vx_tilda = Array(_ny, _nx); 
        _vy_tilda = Array(_ny, _nx);
        _vx_shifted = Array(_ny, _nx + 1);
        _vy_shifted = Array(_ny + 1, _nx);
    }

    //------------------------------------------------------------------------

    Array* get_masses_A() {
        return &_mass_A;
    }

    Array* get_masses_B() {
        return &_mass_B;
    }

    //------------------------------------------------------------------------

    void add_particles_A_evenly_distribute(int particles_number, double particle_mass, double adiabatic_index) {
        _gamma_A = adiabatic_index;
        _particles_A = Particles(particles_number);
        _particles_A.set_mass_for_each(particle_mass);
        _particles_A.set_borders(_size * _nx, _size * _ny, _size);
        _particles_A.evenly_distribute();
        _particles_A.add_random_movements(0.05);
    }
    
    void add_particles_B_evenly_distribute(int particles_number, double particle_mass, double adiabatic_index) {
        _gamma_B = adiabatic_index;
        _particles_B = Particles(particles_number);
        _particles_B.set_mass_for_each(particle_mass);
        _particles_B.set_borders(_size * _nx, _size * _ny, _size);
        _particles_B.evenly_distribute();
        _particles_B.add_random_movements(0.05);
    }

    //------------------------------------------------------------------------

    void initiate_energy_test_function_01();

    void initiate_energy_test_function_02();

    //------------------------------------------------------------------------

    double get_tau_max() {
        double vx_max = _vx.absmax();
        double vy_max = _vy.absmax();
        double softening = 1.0;
        return _size / (vx_max + vy_max + softening);
    }

    void re_pressure() {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                double sA = (_gamma_A - 1) * _mass_A(y, x) * _energy_A(y, x);
                double sB = (_gamma_B - 1) * _mass_B(y, x) * _energy_B(y, x);
                
                _pressure(y, x) = (sA + sB) / (_size * _size);
            }
        }

        for (int y = 0; y < _ny; y++) {
            _pressure_shifted_x(y, 0) = _pressure(y, 0);
            for (int x = 1; x < _nx; x++) {
                _pressure_shifted_x(y, x) = (_pressure(y, x - 1) + _pressure(y, x)) / 2;
            }
            _pressure_shifted_x(y, _nx) = _pressure(y, _nx - 1);
        }

        for (int x = 0; x < _nx; x++) {
            _pressure_shifted_y(0, x) = _pressure(0, x);
        }
        for (int y = 1; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                _pressure_shifted_y(y, x) = (_pressure(y - 1, x) + _pressure(y, x)) / 2;
            }
        }
        for (int x = 0; x < _nx; x++) {
            _pressure_shifted_y(_ny, x) = _pressure(_ny - 1, x);
        }
    }

    void re_v_tilda(double tau) {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                double rho = (_mass_A(y, x) + _mass_B(y, x)) / (_size * _size) + _eps;
                _vx_tilda(y, x) = _vx(y, x) + \
                    tau / rho / _size * (_pressure_shifted_x(y, x) - _pressure_shifted_x(y, x + 1));
                _vy_tilda(y, x) = _vy(y, x) + \
                    tau / rho / _size * (_pressure_shifted_y(y, x) - _pressure_shifted_y(y + 1, x));
            }
        }
    }

    void re_w_energy() {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                double potential = _mass_A(y, x) * _energy_A(y, x) + _mass_B(y, x) * _energy_B(y, x);
                double kinetic = (_mass_A(y, x) + _mass_B(y, x)) * \
                    (_vx(y, x)*_vx(y, x) + _vy(y, x)*_vy(y, x)) / 2;
                _w_energy(y, x) = (potential + kinetic) / (_size*_size);
            }
        }
    }

    void re_v_shifted() {
        _vx_shifted.zeros();
        for (int y = 0; y < _ny; y++) {
            _vx_shifted(y, 0) = 0;
            for (int x = 1; x < _nx; x++) {
                _vx_shifted(y, x) = (_vx(y, x - 1) + _vx(y, x) + _vx_tilda(y, x - 1) + _vx_tilda(y, x)) / 4;
            }
            _vx_shifted(y, _nx) = 0;
        }

        _vy_shifted.zeros();
        for (int x = 0; x < _nx; x++) {
            _vy_shifted(0, x) = 0;
        }
        for (int y = 1; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                _vy_shifted(y, x) = (_vy(y - 1, x) + _vy(y, x) + _vy_tilda(y - 1, x) + _vy_tilda(y, x)) / 4;
            }
        }
        for (int x = 0; x < _nx; x++) {
            _vy_shifted(_ny, x) = 0;
        }
    }

    void re_w_energy_tilda(double tau) {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                _w_energy_tilda(y, x) = _w_energy(y, x) + tau / _size * ( \
                    _vx_shifted(y, x) * _pressure_shifted_x(y, x) + \
                    _vy_shifted(y, x) * _pressure_shifted_y(y, x) - \
                    _vx_shifted(y, x + 1) * _pressure_shifted_x(y, x + 1) - \
                    _vy_shifted(y + 1, x) * _pressure_shifted_y(y + 1, x));
            }
        }
    }

    void re_energy_tilda() {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                double rho = (_mass_A(y, x) + _mass_B(y, x)) / (_size * _size) + _eps;
                double energy_tilda = _w_energy_tilda(y, x) / rho - \
                    (_vx_tilda(y, x)*_vx_tilda(y, x) + _vy_tilda(y, x)*_vy_tilda(y, x)) / 2;
                double delta_energy = energy_tilda - \
                    (_energy_A(y, x) * _mass_A(y, x) + _energy_B(y, x) * _mass_B(y, x)) / (_mass_A(y, x) + _mass_B(y, x));
                _energy_A_tilda(y, x) = _energy_A(y, x) + delta_energy;                
                _energy_B_tilda(y, x) = _energy_B(y, x) + delta_energy;    
            }
        }
    }

    void move_particles(double tau) {
        _particles_A.set_energies(_energy_A_tilda, _vx_tilda, _vy_tilda);
        _particles_B.set_energies(_energy_B_tilda, _vx_tilda, _vy_tilda);

        _particles_A.move_particles(_vx, _vy, _vx_tilda, _vy_tilda, tau);
        _particles_B.move_particles(_vx, _vy, _vx_tilda, _vy_tilda, tau);
    }

    void re_v_and_mass() {
        Array impuls_x = Array(_ny, _nx);
        Array impuls_y = Array(_ny, _nx);

        _mass_A.zeros();
        for (int i = 0; i < _particles_A.len(); i++) {
            int x = std::floor(_particles_A(i, 0) / _size);
            int y = std::floor(_particles_A(i, 1) / _size);
            if (x == -2147483648) {
                std::cout << _particles_A(i, 0) << '\n';
            }

            _mass_A(y, x) += _particles_A(i, 2);
            impuls_x(y, x) += _particles_A(i, 3) * _particles_A(i, 2);
            impuls_y(y, x) += _particles_A(i, 4) * _particles_A(i, 2);
        }

        _mass_B.zeros();
        for (int i = 0; i < _particles_B.len(); i++) {
            int x = std::floor(_particles_B(i, 0) / _size);
            int y = std::floor(_particles_B(i, 1) / _size);

            _mass_B(y, x) += _particles_B(i, 2);
            impuls_x(y, x) += _particles_B(i, 3) * _particles_B(i, 2);
            impuls_y(y, x) += _particles_B(i, 4) * _particles_B(i, 2);
        }

        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                _vx(y, x) = impuls_x(y, x) / (_mass_A(y, x) + _mass_B(y, x) + _eps);
                _vy(y, x) = impuls_y(y, x) / (_mass_A(y, x) + _mass_B(y, x) + _eps);
            }
        }

        impuls_x.free_space();
        impuls_y.free_space();
    }

    void re_energy() {
        Array w_tmp(_ny, _nx);

        for (int i = 0; i < _particles_A.len(); i++) {
            int x = std::floor(_particles_A(i, 0) / _size);
            int y = std::floor(_particles_A(i, 1) / _size);

            w_tmp(y, x) += _particles_A(i, 5); // NOT DIVIDING BY H^2
        }

        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                _energy_A(y, x) = w_tmp(y, x) / (_mass_A(y, x) + _eps) - (_vx(y, x)*_vx(y, x) + _vy(y, x)*_vy(y, x)) / 2;
            }
        }

        w_tmp.zeros();

        for (int i = 0; i < _particles_B.len(); i++) {
            int x = std::floor(_particles_B(i, 0) / _size);
            int y = std::floor(_particles_B(i, 1) / _size);

            w_tmp(y, x) += _particles_B(i, 5); // NOT DIVIDING BY H^2
        }

        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                _energy_B(y, x) = w_tmp(y, x) / (_mass_B(y, x) + _eps) - (_vx(y, x)*_vx(y, x) + _vy(y, x)*_vy(y, x)) / 2;
            }
        }

        w_tmp.free_space();
    }

    //------------------------------------------------------------------------

    void print_masses() {
        std::cout << "# masses A:\n";
        _mass_A.print();
        std::cout << "# sum :" << _mass_A.sum() << "\n";
        std::cout << "# masses B:\n";
        _mass_B.print();
        std::cout << "# sum :" << _mass_B.sum() << "\n";
    }

    void print_pressures(bool shifted = false) {
        std::cout << "# pressure :\n";
        _pressure.print();
        if (shifted) {
            std::cout << "# shifted x :\n";
            _pressure_shifted_x.print();
            std::cout << "# shifted y\n";
            _pressure_shifted_y.print();
        }
    }

    void print_v_tilda() {
        std::cout << "vx tilda :\n";
        _vx_tilda.print();
        std::cout << "vy tilda :\n";
        _vy_tilda.print();
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

    void print_energies(bool full = false) {
        if (full) {
            std::cout << "# energies :\n";
            _w_energy.print();
        }
        std::cout << _w_energy.sum() << "\n";
    }

    //------------------------------------------------------------------------

    void write_energies(int shot) {
        std::string name;
        name = "./results/energy_a" + std::to_string(shot) + ".txt";
        _energy_A.write_in_file(name);
        name = "./results/energy_b" + std::to_string(shot) + ".txt";
        _energy_B.write_in_file(name);
    }

    void write_w_energies(int shot) {
        std::string name;
        name = "./results/energy" + std::to_string(shot) + ".txt";
        _w_energy.write_in_file(name);
        name = "./results/energy_tilda" + std::to_string(shot) + ".txt";
        _w_energy_tilda.write_in_file(name);
    }

    void write_masses(int shot) {
        std::string name;
        name = "./results/mass_a" + std::to_string(shot) + ".txt";
        _mass_A.write_in_file(name);
        name = "./results/mass_b" + std::to_string(shot) + ".txt";
        _mass_B.write_in_file(name);
    }

    void write_pressure(int shot) {
        std::string name;
        name = "./results/pressure" + std::to_string(shot) + ".txt";
        _pressure.write_in_file(name);
    }

    void write_v(int shot) {
        std::string name;
        name = "./results/vx" + std::to_string(shot) + ".txt";
        _vx.write_in_file(name);
        name = "./results/vy" + std::to_string(shot) + ".txt";
        _vy.write_in_file(name);
    }

    void write_v_tilda(int shot) {
        std::string name;
        name = "./results/vx_tilda" + std::to_string(shot) + ".txt";
        _vx_tilda.write_in_file(name);
        name = "./results/vy_tilda" + std::to_string(shot) + ".txt";
        _vy_tilda.write_in_file(name);
    }

private:
    double _size;
    double _eps = 1e-12;
    int _nx;
    int _ny;

    Particles _particles_A;
    Array _mass_A;
    Array _energy_A;
    Array _energy_A_tilda; 
    double _gamma_A;

    Particles _particles_B;
    Array _mass_B;
    Array _energy_B;
    Array _energy_B_tilda;
    double _gamma_B;

    Array _pressure;
    Array _pressure_shifted_x;
    Array _pressure_shifted_y;

    Array _w_energy;
    Array _w_energy_tilda;

    Array _vx;
    Array _vy;
    Array _vx_tilda;
    Array _vy_tilda;
    Array _vx_shifted;
    Array _vy_shifted;
};