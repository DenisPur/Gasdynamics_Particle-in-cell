#include <iostream>
#include <string>
#include <fstream>
#include "array.cpp"
#include "particles.cpp"
#include "board.cpp"

void Board_2IdealGases::initiate_energy_test_function_01() {
    for (int y = 0; y < _ny; y++) {
        for (int x = 0; x < _nx; x++) {
            _energy_A(y, x) = 0.005 * (x * y + 2);
            _energy_B(y, x) = 0.005 * ((_nx - x) * (_ny - y) + 2) ;
        }
    }
}

void make_n_steps(Board_2IdealGases* b, int n) {
    double tau = 0.0001;
    double time = tau;
    int shot = 0;

    std::ofstream timestems("./results/timesteps.txt");

    for (int i = 0; i < n; i++) {
        b->re_pressure();
        b->re_v_tilda(tau);

        tau = b->get_tau_max() / 2;
        time += tau;
        
        std::cout << "# tau :" << tau << "\n";
        
        b->re_w_energy();
        b->re_v_shifted();
        b->re_w_energy_tilda(tau);
        b->re_energy_tilda();
        b->move_particles(tau);
        b->re_v_and_mass();
        b->re_energy();
        
        if (i % 10 == 0) {
            timestems << time << '\n';
            b->write_masses(shot);
            b->write_energies(shot);
            b->write_pressure(shot);
            b->write_v(shot);
            shot++;
        }
    }
    timestems.close();
}

int main() {
    Board_2IdealGases b(200, 0.1);

    b.add_particles_A(1000000, 0.0000001, 1.67);
    b.add_particles_B(1000000, 0.0000002, 1.4);
    b.initiate_energy_test_function_01();
    b.re_mass();

    make_n_steps(&b, 101);

    return 0;
}
