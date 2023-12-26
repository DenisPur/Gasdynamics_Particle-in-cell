#include <iostream>

#include "array.cpp"
#include "particles.cpp"
#include "board.cpp"

int main() {
    Board b(10, 0.1);
    b.add_particles(1000, 1.0);
    // b.print_masses();

    return 0;
}
