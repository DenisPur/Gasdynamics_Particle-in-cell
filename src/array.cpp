#pragma once
#include <iostream>

class Array {
public:
    Array() {
        _nx = 0;
        _ny = 0;
        _ptr_zero = nullptr;
    }

    Array(int y, int x) {
        _nx = x;
        _ny = y;
        _ptr_zero = new double [_nx * _ny];
        zeros();
    }

    double& operator()(int y, int x) {
        return *(_ptr_zero + _nx * y  + x);
    }

    int len_x() {
        return _nx; 
    }

    int len_y() {
        return _ny;
    }

    double sum() {
        double sum = 0;
        for (int i = 0; i < _nx * _ny; i++) {
            sum += *(_ptr_zero + i);
        }
        return sum;
    }

    void print() {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                std::cout << operator()(y, x) << ' ';
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
