#pragma once
#include <vector>
#include <tuple>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include "pnc.cpp"

class Board {
public:
    Board(int x_size, int y_size, double cell_size) {
        _num_x = x_size;
        _num_y = y_size;
        _cell_size = cell_size;

        _cells.resize(_num_y);
        for (std::vector<Cell> &row : _cells) {
            row.resize(_num_x);
        }

        for (int ix = 0; ix < _num_x; ix++) {
            for (int iy = 0; iy < _num_y; iy++) {
                _cells[iy][ix].set_indexes(ix, iy);
            }
        }
    }

    Cell *get_cell_by_coordinates(double x, double y) {
        if ((x < 0) || (x > _num_x * _cell_size) || (y < 0) || (y > _num_y * _cell_size)) {
            std::cout << "x=" << x << " y=" << y << " not in ";
            std::cout << "[0, " << _num_x * _cell_size << "], [0, " << _num_y * _cell_size << "]" << "\n";
            throw std::invalid_argument( "get_cell() : (x, y) out of the bounds" );
        }
        if (x == _num_x * _cell_size) {
            x -= 1e-9;
        }
        if (y == _num_y * _cell_size) {
            y -= 1e-9;
        }

        int pos_x = std::floor(x / _cell_size);
        int pos_y = std::floor(y / _cell_size);
        return & _cells[pos_y][pos_x]; 
    }

    Cell *get_cell_by_number(int ix, int iy) {
        return & _cells[iy][ix];
    }

    void reserve_particles(int n) {
        _particles.reserve(n);
    }

    void add_default_particle(double x, double y) {
        Particle *particle = & _particles.emplace_back(Particle{});
        particle->set_coordinates(x, y);
        Cell *cell = get_cell_by_coordinates(x, y);
        particle->link_cell(cell);
        cell->add_particle(particle);
    }

    void move_particle(Particle* particle, double new_x, double new_y) {
        Cell* new_cell = get_cell_by_coordinates(new_x, new_y);
        particle->get_cell()->delete_particle(particle);
        new_cell->add_particle(particle);
        particle->link_cell(new_cell);
    }

    std::vector<Particle>* get_particles() {
        return & _particles;
    }

private:
    std::vector<std::vector<Cell>> _cells;
    std::vector<Particle> _particles;
    int _num_x;
    int _num_y;
    double _cell_size;
};
