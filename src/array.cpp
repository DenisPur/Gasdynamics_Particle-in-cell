#pragma once
#include <iostream>

class Array {
public:
    Array() {
        _nx = 0;
        _ny = 0;
        _ptr_zero = nullptr;
    }

    Array(int x, int y) {
        _nx = x;
        _ny = y;
        _ptr_zero = new double [x*y];
        zeros();
    }

    double& operator()(int x, int y) {
        return *(_ptr_zero + _ny * x  + y);
    }

    int len_x() {
        return _nx; 
    }

    int len_y() {
        return _ny;
    }

    void print() {
        for (int i = 0; i < _nx; i++) {
            for (int j = 0; j < _ny; j++) {
                std::cout << operator()(i, j) << ' ';
            }
            std::cout << '\n';
        }
    }

    void zeros() {
        for (int i = 0; i < _nx * _ny; i++) {
            *(_ptr_zero + i) = 0.0;
        }
    }

private:
    int _nx;
    int _ny;
    double* _ptr_zero;
};
