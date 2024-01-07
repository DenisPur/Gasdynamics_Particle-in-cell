#include <iostream>
#include <string>
#include <fstream>
#include "array.cpp"
#include "particles.cpp"
#include "board.cpp"

void Board_2IdealGases::initiate_energy_test_function_01() {
    add_particles_A_evenly_distribute(1000000, 0.000001, 1.67);
    add_particles_B_evenly_distribute(1000000, 0.000002, 1.4);

    for (int y = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++) {
            _energy_A(y, x) = 0.005 * (x * y) + 2;
            _energy_B(y, x) = 0.005 * (_nx - x) * (_ny - y) + 2 ;
        }
    }
}

void Board_2IdealGases::initiate_energy_test_function_02() {
    _gamma_A = 1.67;
    _particles_A = Particles(1000000);
    _particles_A.set_mass_for_each(0.000001);
    _particles_A.set_borders(_size * _nx, _size * _ny, _size);

    _gamma_B = 1.4;
    _particles_B = Particles(1000000);
    _particles_B.set_mass_for_each(0.000002);
    _particles_B.set_borders(_size * _nx, _size * _ny, _size);

    _particles_A.evenly_distribute(0, 800000);
    _particles_A.add_random_movements(0, 800000);

    _particles_B.evenly_distribute(0, 800000);
    _particles_B.add_random_movements(0, 800000);

    for (int i = 800000; i < 1000000; i++) {
        double x;
        double y;
        do {
            x = 0.3 * rand() / RAND_MAX - 0.15;
            y = 20.0 * rand() / RAND_MAX;
            x = x + y;
        } while (x < 0 || x > 20.0);
        _particles_A.set_coordinates(i, x, y);
        _particles_A.set_velocity(i, 1.0, -1.0);
    }

    for (int i = 800000; i < 1000000; i++) {
        double x;
        double y;
        do {
            x = 0.3 * rand() / RAND_MAX - 0.15;
            y = 20.0 * rand() / RAND_MAX;
            x = x + y;
        } while (x < 0 || x > 20.0);
        _particles_B.set_coordinates(i, x, y);
        _particles_B.set_velocity(i, 1.0, -1.0);
    }

    for (int y = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++) {
            _energy_A(y, x) = 0.1 + 0.03 * rand() / RAND_MAX;
            _energy_B(y, x) = 0.1 + 0.03 * rand() / RAND_MAX ;
        }
    }
}

void make_n_steps(Board_2IdealGases* b, int n) {
    double tau = 0.0001;
    double time = tau;
    int shot = 0;

    for (int i = 0; i < n; i++) {
        tau = b->get_tau_max();
        time += tau;        
        std::cout << "# tau :" << tau << "\n";

        b->re_pressure();
        b->re_v_tilda(tau);
        b->re_w_energy();
        b->re_v_shifted();
        b->re_w_energy_tilda(tau);
        b->re_energy_tilda();
        b->move_particles(tau);
        b->re_v_and_mass();
        b->re_energy();
        
        if (i % 1 == 0) {
            std::ofstream timestems("./results/timesteps.txt", std::ios::app);
            timestems << time << '\n';
            timestems.close();
            
            b->write_masses(shot);
            b->write_energies(shot);
            b->write_w_energies(shot);
            b->write_pressure(shot);
            b->write_v(shot);
            b->write_v_tilda(shot);
            shot++;
        }
    }
}

int main() {
    Board_2IdealGases b(200, 0.1);

    b.initiate_energy_test_function_02();
    b.re_v_and_mass();

    make_n_steps(&b, 10);

    return 0;
}
