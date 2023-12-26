#include <iostream>

#include "array.cpp"
#include "particles.cpp"
#include "board.cpp"

int main() {
    Board_2IdealGases b(5, 0.1);

    b.add_particlesA(100, 1.0);
    b.add_particlesB(1000, 0.1);
    b.re_mass();

    b.initiate_energy_test_function();
    
    b.re_pressure();
    b.print_pressures();

    b.re_v_tilted(0.01);
    b.print_v_tilted();

    b.re_w_energies();

    b.re_v_shifted();

    b.re_w_energies_tilted(0.01);

    b.print_energies_tilted();

    return 0;
}
