#pragma once
#include <vector>
#include <tuple>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <algorithm>

class Cell;

class Particle {
public:
    void set_coordinates(double x, double y) {
        _x = x;
        _y = y;
    }

    void link_cell(Cell *cell) {
        _cell = cell;
    }

    Cell* get_cell() {
        return _cell;
    }

private:
    Cell *_cell;
    double _x;
    double _y;
    int _type;
    double _m;
};

class Cell {
public:
    void set_indexes(int x, int y) {
        _ix = x;
        _iy = y;
    }

    void set_speed(double vx, double vy) {
        _vx = vx;
        _vy = vy;
    }

    void add_particle(Particle *particle) {
        _particles.push_back(particle);
    }

    void delete_particle(Particle *particle) {
        std::vector<Particle*>::iterator it = std::find(_particles.begin(), _particles.end(), particle);
        if (it == _particles.end()) {
            throw std::invalid_argument( "delete_particle() : the Particle not found in the Cell" );
        }
        _particles.erase(it);
    }

    std::tuple<int, int> get_indexes() {
        return std::tuple<int, int> {_ix, _iy};
    }

    int get_particles_number() {
        return _particles.size();
    }

private:
    std::vector<Particle*> _particles;
    int _ix;
    int _iy;
    double _vx;
    double _vy;
};