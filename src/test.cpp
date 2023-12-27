#include <iostream>

#include "array.cpp"
#include "particles.cpp"
#include "board.cpp"

void make_n_steps(Board_2IdealGases* b, int n, double tau) {
    for (int i = 0; i < n; i++) {
        b->re_pressure();
        b->re_v_tilda(tau);
        b->re_w_energy();
        b->re_v_shifted();
        b->re_full_energy_tilda(tau);
        b->re_energy_tilda();
        b->move_partivles(tau);
        b->re_v_and_mass();
        b->re_energy();
    }
}

int main() {
    Board_2IdealGases b(5, 0.1);

    b.add_particles_A(100, 0.001, 1.67);
    b.add_particles_B(1000, 0.0001, 1.4);
    b.initiate_energy_test_function_01();
    b.re_mass();

    b.print_masses();

    make_n_steps(&b, 10, 0.01);

    b.print_masses();

    return 0;
}
