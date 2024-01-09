#pragma once
#include <iostream>
#include <fstream>
#include <string>

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
        _ptr_zero = new double[_nx * _ny];
        zeros();
    }

    // ~Array() {
    //     delete [] _ptr_zero;
    // }

    void free_space() {
        delete [] _ptr_zero;
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

    double absmax() {
        double max = 0;
        for (int i = 0; i < _nx * _ny; i++) {
            if (std::abs(*(_ptr_zero + i)) > max) {
                max = std::abs(*(_ptr_zero + i));
            }
        }
        return max;
    }

    double min() {
        double min = 0;
        for (int i = 0; i < _nx * _ny; i++) {
            if (*(_ptr_zero + i) < min) {
                min = (*(_ptr_zero + i));
            }
        }
        return min;
    }

    //------------------------------------------------------------------------


    void zeros() {
        for (int i = 0; i < _nx * _ny; i++) {
            *(_ptr_zero + i) = 0.0;
        }
    }

    //------------------------------------------------------------------------

    void write_in_file(std::string name) {
        std::ofstream out(name);

        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx - 1; x++) {
                out << operator()(y, x) << ',';
            }
            out << operator()(y, _nx - 1);
            out << std::endl;
        }
        out.close();
    }

    void print() {
        for (int y = 0; y < _ny; y++) {
            for (int x = 0; x < _nx; x++) {
                std::cout << operator()(y, x) << ' ';
            }
            std::cout << '\n';
        }
    }

private:
    int _nx;
    int _ny;
    double* _ptr_zero;
};