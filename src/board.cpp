#pragma once
#include <iostream>
#include <cmath>
#include <algorithm>
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

    void experiment_01_initial_state() {
        _gamma_A = 1.67;
        _particles_A = Particles(1000000, _nx, _ny, _size);
        _particles_A.set_mass_for_each(1.0E-6);
        _particles_A.evenly_distribute();
        _particles_A.add_random_movements(0.05);

        _gamma_B = 1.4;
        _particles_B = Particles(1000000, _nx, _ny, _size);
        _particles_B.set_mass_for_each(2.0E-6);
        _particles_B.evenly_distribute();
        _particles_B.add_random_movements(0.05);

        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                _energy_A(y, x) = 0.05 * (x * y) + 2;
                _energy_B(y, x) = 0.05 * (_nx - x) * (_ny - y) + 2 ;
            }
        }

        set_minimal_mass(1.0E-6);
        prepare_for_run();
    }

    void experiment_02_initial_state() {
        // Diagonal shok wave with 25% of all particles in it

        _gamma_A = 1.67;
        _particles_A = Particles(1000000, _nx, _ny, _size);
        _particles_A.set_mass_for_each(1.0E-6);

        _gamma_B = 1.4;
        _particles_B = Particles(1000000, _nx, _ny, _size);
        _particles_B.set_mass_for_each(2.0E-6);

        // 750'000 particles distributed randomly
        _particles_A.evenly_distribute(0, 750000);
        _particles_A.add_random_movements(0.1, 0, 750000);
        // 250'000 pu on diagonal and given momentum
        for (int i = 750000; i < 1000000; i++) {
            double p;
            double x;
            double y;
            do {
                p = 1.0 * rand() / RAND_MAX;
                x = 5 * _size * 2 * (p - 0.5);
                y = _size * _ny * rand() / RAND_MAX;
                x = x + y;
            } while (x < 0 || x > _size * _nx);
            _particles_A.set_coordinates(i, x, y);
            _particles_A.set_velocity(i, p, -p);
        }

        // Same for "B" particles
        _particles_B.evenly_distribute(0, 750000);
        _particles_B.add_random_movements(0.1, 0, 750000);
        for (int i = 750000; i < 1000000; i++) {
            double p;   
            double x;
            double y;
            do {
                p = 1.0 * rand() / RAND_MAX;
                x = 5 * _size * 2 * (p - 0.5);
                y = _size * _ny * rand() / RAND_MAX;
                x = x + y;
            } while (x < 0 || x > _size * _nx);
            _particles_B.set_coordinates(i, x, y);
            _particles_B.set_velocity(i, p, -p);
        }

        // For all cells inner enrgy is roughly the same
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                _energy_A(y, x) = 1 + 0.01 * rand() / RAND_MAX;
                _energy_B(y, x) = 1 + 0.01 * rand() / RAND_MAX ;
            }
        }

        set_minimal_mass(1.0E-6);
        prepare_for_run();
    }

    //------------------------------------------------------------------------

    void set_minimal_mass(double minimal_mass) {
        _minim_mass = minimal_mass;
    }

    bool cell_is_empty(int y, int x) {
        return _mass_A(y, x) + _mass_B(y, x) < _minim_mass;
    }

    void prepare_for_run() {
        s08_v_and_mass();
    }

    double sqr_sum(double a, double b) {
        return a*a + b*b;
    }


    //------------------------------------------------------------------------

    double get_tau_max() {
        double vx_max = _vx.absmax();
        double vy_max = _vy.absmax();
        double softening = 1.0;
        return _size / (vx_max + vy_max + softening);
    }

    void s01_pressure() {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                double sA = (_gamma_A - 1) * _mass_A(y, x) * _energy_A(y, x);
                double sB = (_gamma_B - 1) * _mass_B(y, x) * _energy_B(y, x);
                
                _pressure(y, x) = std::max(0.0, (sA + sB) / (_size * _size));
            }
        }

        for (int y = 0; y < _ny; y++) {
            _pressure_shifted_x(y, 0) = _pressure(y, 0);
            for (int x = 1; x < _nx; x++) {
                if (_pressure(y, x - 1) * _pressure(y, x) <= 0) {
                    _pressure_shifted_x(y, x) = 0;
                } else {
                    _pressure_shifted_x(y, x) = (_pressure(y, x - 1) + _pressure(y, x)) / 2;
                }
            }
            _pressure_shifted_x(y, _nx) = _pressure(y, _nx - 1);
        }

        for (int x = 0; x < _nx; x++) {
            _pressure_shifted_y(0, x) = _pressure(0, x);
        }
        for (int y = 1; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                if (_pressure(y - 1, x) * _pressure(y, x) <= 0) {
                    _pressure_shifted_y(y, x) = 0;
                } else {
                    _pressure_shifted_y(y, x) = (_pressure(y - 1, x) + _pressure(y, x)) / 2;
                }
            }
        }
        for (int x = 0; x < _nx; x++) {
            _pressure_shifted_y(_ny, x) = _pressure(_ny - 1, x);
        }
    }

    void s02_v_tilda(double tau) {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                if (cell_is_empty(y, x)) {
                    _vx_tilda(y, x) = 0;
                    _vy_tilda(y, x) = 0;
                    continue;
                }
                double rho = (_mass_A(y, x) + _mass_B(y, x)) / (_size * _size);
                _vx_tilda(y, x) = _vx(y, x) + tau / rho / _size * (
                    _pressure_shifted_x(y, x) - _pressure_shifted_x(y, x + 1));
                _vy_tilda(y, x) = _vy(y, x) + tau / rho / _size * (
                    _pressure_shifted_y(y, x) - _pressure_shifted_y(y + 1, x));
            }
        }
    }

    void s03_w_energy() {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                double potential = _mass_A(y, x) * _energy_A(y, x) + _mass_B(y, x) * _energy_B(y, x);
                double kinetic = (_mass_A(y, x) + _mass_B(y, x)) * sqr_sum(_vx(y, x), _vy(y, x)) / 2;
                _w_energy(y, x) = (potential + kinetic) / (_size*_size);
            }
        }
    }

    void s04_v_shifted() {
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

    void s05_w_energy_tilda(double tau) {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                _w_energy_tilda(y, x) = _w_energy(y, x) + tau / _size * 
                    ( _vx_shifted(y, x) * _pressure_shifted_x(y, x)
                    + _vy_shifted(y, x) * _pressure_shifted_y(y, x)
                    - _vx_shifted(y, x + 1) * _pressure_shifted_x(y, x + 1)
                    - _vy_shifted(y + 1, x) * _pressure_shifted_y(y + 1, x) );
            }
        }
    }

    void s06_energy_tilda() {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                if (cell_is_empty(y, x)) {
                    _energy_A_tilda(y, x) = 0;
                    _energy_B_tilda(y, x) = 0;
                    _vx_tilda(y, x) = 0;
                    _vy_tilda(y, x) = 0;
                    continue;
                }

                double rho = (_mass_A(y, x) + _mass_B(y, x)) / (_size * _size);
                double energy_tilda = _w_energy_tilda(y, x) / rho - sqr_sum(_vx_tilda(y, x), _vy_tilda(y, x)) / 2;
                double energy = _w_energy(y, x) / rho - sqr_sum(_vx(y, x), _vy(y, x)) / 2;
                double delta_energy = energy_tilda - energy;

                _energy_A_tilda(y, x) = _energy_A(y, x) + delta_energy;                
                _energy_B_tilda(y, x) = _energy_B(y, x) + delta_energy;    
            }
        }
    }

    void s07_move_particles(double tau) {
        _particles_A.set_energies_and_move(_vx, _vy, _vx_tilda, _vy_tilda, _energy_A_tilda, tau);
        _particles_B.set_energies_and_move(_vx, _vy, _vx_tilda, _vy_tilda, _energy_B_tilda, tau);
    }

    void s08_v_and_mass() {
        Array impuls_x = Array(_ny, _nx);
        Array impuls_y = Array(_ny, _nx);

        _mass_A.zeros();
        for (int i = 0; i < _particles_A.len(); i++) {
            int x = std::floor(_particles_A(i, 0) / _size);
            int y = std::floor(_particles_A(i, 1) / _size);

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
                if (cell_is_empty(y, x)) {
                    _vx(y, x) = 0;
                    _vy(y, x) = 0;
                } else {
                    _vx(y, x) = impuls_x(y, x) / (_mass_A(y, x) + _mass_B(y, x));
                    _vy(y, x) = impuls_y(y, x) / (_mass_A(y, x) + _mass_B(y, x));
                }
            }
        }

        impuls_x.free_space();
        impuls_y.free_space();
    }

    void s09_energy() {
        Array w_tmp_A(_ny, _nx);
        Array w_tmp_B(_ny, _nx);

        for (int i = 0; i < _particles_A.len(); i++) {
            int x = std::floor(_particles_A(i, 0) / _size);
            int y = std::floor(_particles_A(i, 1) / _size);

            w_tmp_A(y, x) += _particles_A(i, 5); // NOT DIVIDING BY H^2
        }

        for (int i = 0; i < _particles_B.len(); i++) {
            int x = std::floor(_particles_B(i, 0) / _size);
            int y = std::floor(_particles_B(i, 1) / _size);

            w_tmp_B(y, x) += _particles_B(i, 5); // NOT DIVIDING BY H^2
        }

        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                if (_mass_A(y, x) < _minim_mass) {
                    _energy_A(y, x) = 0;
                } else {
                    _energy_A(y, x) = w_tmp_A(y, x) / _mass_A(y, x) - sqr_sum(_vx(y, x), _vy(y, x)) / 2;
                }

                if (_mass_B(y, x) < _minim_mass) {
                    _energy_B(y, x) = 0;
                } else {
                    _energy_B(y, x) = w_tmp_B(y, x) / _mass_B(y, x) - sqr_sum(_vx(y, x), _vy(y, x)) / 2;
                }
                
                if (_energy_A(y, x) < 0 || _energy_B(y, x) < 0) {
                    double energy = (w_tmp_A(y, x) + w_tmp_B(y, x)) / (_mass_A(y, x) + _mass_B(y, x)) 
                        - sqr_sum(_vx(y, x), _vy(y, x)) / 2;
                    _energy_A(y, x) = energy;
                    _energy_B(y, x) = energy;       
                }
            }
        }

        w_tmp_A.free_space();
        w_tmp_B.free_space();
    }

    void make_one_step(double tau) {
        s01_pressure();
        s02_v_tilda(tau);
        s03_w_energy();
        s04_v_shifted();
        s05_w_energy_tilda(tau);
        s06_energy_tilda();
        s07_move_particles(tau);
        s08_v_and_mass();
        s09_energy();
    }

    //------------------------------------------------------------------------

    void write_energies(std::string folder, int shot) {
        std::string name;
        name = "./" + folder + "/energy_a" + std::to_string(shot) + ".txt";
        _energy_A.write_in_file(name);
        name = "./" + folder + "/energy_b" + std::to_string(shot) + ".txt";
        _energy_B.write_in_file(name);
    }

    void write_w_energies(std::string folder, int shot) {
        std::string name;
        name = "./" + folder +"/energy" + std::to_string(shot) + ".txt";
        _w_energy.write_in_file(name);
        name = "./" + folder + "/energy_tilda" + std::to_string(shot) + ".txt";
        _w_energy_tilda.write_in_file(name);
    }

    void write_masses(std::string folder, int shot) {
        std::string name;
        name = "./" + folder + "/mass_a" + std::to_string(shot) + ".txt";
        _mass_A.write_in_file(name);
        name = "./" + folder + "/mass_b" + std::to_string(shot) + ".txt";
        _mass_B.write_in_file(name);
    }

    void write_pressure(std::string folder, int shot) {
        std::string name;
        name = "./" + folder + "/pressure" + std::to_string(shot) + ".txt";
        _pressure.write_in_file(name);
    }

    void write_v(std::string folder, int shot) {
        std::string name;
        name = "./" + folder + "/vx" + std::to_string(shot) + ".txt";
        _vx.write_in_file(name);
        name = "./" + folder + "/vy" + std::to_string(shot) + ".txt";
        _vy.write_in_file(name);
    }

    void write_everything(std::string results_folder, int shot) {
        write_masses(results_folder, shot);
        write_energies(results_folder, shot);
        write_w_energies(results_folder, shot);
        write_pressure(results_folder, shot);
        write_v(results_folder, shot);
    }


private:
    double _size;
    int _nx;
    int _ny;
    double _minim_mass = 0;

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