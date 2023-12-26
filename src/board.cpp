#pragma once

#include <cmath>
#include "array.cpp"
#include "particles.cpp"

class Board {
public:
    Board(int n, double cell_size) {
        Board(n, n, cell_size);
    }

    Board(int nx, int ny, double cell_size) {
        _nx = nx;
        _ny = ny;
        _size = cell_size;

        _masses = Array(_nx, _ny);
        _energies = Array(_nx, _ny);
        _vx = Array(_nx, _ny); 
        _vy = Array(_nx, _ny);
    }

    void add_particles(int n, double m) {
        _particles = Particles(n);
        _particles.set_mass_for_each(m);
        _particles.evenly_distribute(_size * _nx, _size * _ny);
        re_mass();
    }

    void re_mass() {
        _masses.zeros();

        for (int i = 0; i < _particles.len(); i++) {
            int nx = std::floor(_particles(i, 0) / _size);
            int ny = std::floor(_particles(i, 1) / _size);
            _masses(nx, ny) += _particles(i, 2);
        }
    }

    void print_masses() {
        _masses.print();
    }

private:
    Particles _particles;

    Array _masses;
    Array _energies;
    Array _vx;
    Array _vy;

    double _size;
    int _nx;
    int _ny;
};