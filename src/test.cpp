#include <iostream>

#include "board.cpp"

void print_tuple(std::tuple<int, int> tuple) {
    std::cout << " x=" << std::get<0>(tuple) << ",  y=" << std::get<1>(tuple) << "\n";
}

int main() {
    Board board(6, 6, 0.5);

    board.reserve_particles(10);
    board.add_default_particle(0, 0);
    board.add_default_particle(2.2, 0.2);
    board.add_default_particle(2.0, 0.9);
    board.add_default_particle(2.1, 0.8);
    board.add_default_particle(0.7, 2.1);

    for (int iy = 5; iy >= 0; iy --) {
        for (int ix = 0; ix < 6; ix ++) {
            std::cout << board.get_cell_by_number(ix, iy)->get_particles_number() << " ";
        }
        std::cout << "\n";
    }

    for (Particle part : *(board.get_particles())) {
        print_tuple(part.get_cell()->get_indexes());
    }

    Particle *particle_to_move = &(*board.get_particles()->begin());
    board.move_particle(particle_to_move, 3.0, 3.0);

    for (int iy = 5; iy >= 0; iy --) {
        for (int ix = 0; ix < 6; ix ++) {
            std::cout << board.get_cell_by_number(ix, iy)->get_particles_number() << " ";
        }
        std::cout << "\n";
    }

    for (Particle part : *(board.get_particles())) {
        print_tuple(part.get_cell()->get_indexes());
    }

    particle_to_move = &(*(board.get_particles()->begin() + 1));
    board.move_particle(particle_to_move, 2.5, 2.4);

    for (int iy = 5; iy >= 0; iy --) {
        for (int ix = 0; ix < 6; ix ++) {
            std::cout << board.get_cell_by_number(ix, iy)->get_particles_number() << " ";
        }
        std::cout << "\n";
    }

    for (Particle part : *(board.get_particles())) {
        print_tuple(part.get_cell()->get_indexes());
    }

    return 0;
}