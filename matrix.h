/*
    # Matrix class by Mikayel Kenjetsyan
    # Copyright(c) 2023 Mikayel Kenjetsyan. All rights reserved.
    # This code is provided for educational purposes and personal use only.
    # Redistribution and use in source and binary forms, with or without modification,
    # are not permitted without the express permission of the author.

            ███╗░░░███╗░█████╗░████████╗██████╗░██╗██╗░░██╗
            ████╗░████║██╔══██╗╚══██╔══╝██╔══██╗██║╚██╗██╔╝
            ██╔████╔██║███████║░░░██║░░░██████╔╝██║░╚███╔╝░
            ██║╚██╔╝██║██╔══██║░░░██║░░░██╔══██╗██║░██╔██╗░
            ██║░╚═╝░██║██║░░██║░░░██║░░░██║░░██║██║██╔╝╚██╗
            ╚═╝░░░░░╚═╝╚═╝░░╚═╝░░░╚═╝░░░╚═╝░░╚═╝╚═╝╚═╝░░╚═╝
    
*/
 
#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstddef>
#include <initializer_list>


// double f = 2.5;
//int a = 4;
//long u = a + f;

//stupencheti
//rankD
//obratni
//determinant 
//СЛУ
// A * X * C = B

template <typename T>
class Matrix
{
private:
    std::size_t _row;
    std::size_t _column;
    std::vector<std::vector<T>> _matrix;
public:

    /// Merged default constructor, Matrix(row, column) and Matrix(row, column, filler)
    /// (rows, columns, value) constructor, filling by value
    Matrix(std::size_t row = 0, std::size_t column = 0, std::size_t filler = 0) : 
        _row(row)
        , _column(column)
    {
        _matrix.resize(_row, std::vector<T>(_column, filler));
    }

    /// Copy constructor

    Matrix(Matrix<T>& other)
    {
        std::swap(_row, other._row);
        std::swap(_column, other._column);
        std::swap(_matrix, other._matrix);
    }

    /// Initializer list constructor
    Matrix(const std::initializer_list<std::initializer_list<T>>& list)
    {
        _row = list.size();
        if (_row > 0)
            _column = list.begin()->size();
        else
            _column = 0;
        _matrix.resize(_row, std::vector<T>(_column));
        auto it1 = list.begin();
        for (std::size_t i = 0; i < _row; ++i, ++it1)
        {
            auto it2 = it1->begin();
            for (std::size_t j = 0; j < _column; ++j, ++it2)
            {
                _matrix[i][j] = *it2;
            }
        }
    }

    ///Move constructor
    /*Matrix(Matrix&& other) : 
        _row(other._row)
        , _column(other._column)
    {
        for (std::size_t i = 0; i < other._row; ++i)
        {
            for(std::size_t u = )
        }
    }*/
    
    /// Destructor
    ~Matrix() = default;
    
    /// OPERATORS

    /// operator[]
    std::vector<T>& operator[](std::size_t index)
    {
        return _matrix[index];
    }
    const std::vector<T> operator[](std::size_t index) const
    {
        return _matrix[index];
    }

    /// operator+
    template <typename U>
    Matrix<double> operator+(const Matrix<U>& other) const
    {
        try {
            if (_row == other.get_row() && _column == other.get_column())
            {
                Matrix<double> new_matrix(_row, _column);
                for (std::size_t i = 0; i < _row; ++i)
                {
                    for (std::size_t j = 0; j < _column; ++j)
                    {
                        new_matrix[i][j] = _matrix[i][j] + other[i][j];
                    }
                }
                return new_matrix;
            }
            else
            {
                throw std::logic_error("Incompatible matrix dimensions to add\n");
            }
        }
        catch (const std::exception& e)
        {
            std::cout << e.what() << '\n';
        }
    }

    /// operator-
    template <typename U>
    Matrix<double> operator-(const Matrix<U>& other) const
    {
        try {
            if (_row == other.get_row() && _column == other.get_column())
            {
                Matrix<double> new_matrix(_row, _column);
                for (std::size_t i = 0; i < _row; ++i)
                {
                    for (std::size_t j = 0; j < _column; ++j)
                    {
                        new_matrix[i][j] = _matrix[i][j] - other[i][j];
                    }
                }
                return new_matrix;
            }
            else
            {
                throw std::logic_error("Incompatible matrix dimensions for sub");
            }
        }
        catch (const std::exception& e)
        {
            std::cout << e.what() << '\n';
        }
    }

    /// operator* for matrices
    template <typename U>
    Matrix<double> operator*(const Matrix<U>& other) const
    {
        try {
            if (_column != other._row)
            {
                throw std::logic_error("Incompatible matrix dimensions to multiplicate\n");
            }
            else {
                Matrix<double> new_matrix(_row, other._column);

                for (std::size_t i = 0; i < _row; ++i)
                {
                    for (std::size_t j = 0; j < other._column; ++j)
                    {
                        double sum = 0;
                        for (std::size_t u = 0; u < _column; ++u)
                        {
                            sum += _matrix[i][u] * other[u][j];
                        }
                        new_matrix[i][j] = sum;
                    }
                }
                return new_matrix;
            }
        }
        catch (const std::exception& e)
        {
            std::cout << e.what() << '\n';
        }
    }

    /// operator* to multiply by any number || THERE IS A PROBLEM
    /*template <typename U>
    Matrix<double> operator*(U value)
    {
        Matrix<T> new_matrix(_row, _column);
        for (std::size_t i = 0; i < _row; ++i)
        {
            for (std::size_t j = 0; j < _column; ++j)
            {
                new_matrix[i][j] = (_matrix[i][j]) * value;
            }
        }
        return new_matrix;
    }*/

    /// operator= to Matrix
    Matrix<T>& operator=(const Matrix<T>& other) {
        if (this != &other) {
            _row = other._row;
            _column = other._column;
            _matrix.resize(_row, std::vector<T>(_column));
            for (std::size_t i = 0; i < _row; ++i) {
                for (std::size_t j = 0; j < _column; ++j) {
                    _matrix[i][j] = other._matrix[i][j];
                }
            }
        }
        return *this;
    }

    //matr1 += matr2
    template <typename U>
    Matrix<double>& operator+=(const Matrix<U>& other)
    {
        (*this) = (*this) + other;
        return *this;
    }
    //matr1 -= matr2
    template <typename U>
    Matrix<double>& operator-=(const Matrix<U>& other)
    {
        (*this) = (*this) - other;
        return *this;
    }
    //matr1 *= matr2
    template <typename U>
    Matrix<double>& operator*=(const Matrix<U>& other)
    {
        (*this) = (*this) * other;
        return *this;
    }
    //matr1 *= num
    template <typename U>

    Matrix<double>& operator*=(U num)
    {
        (*this) = (*this) * num;
        return *this;
    }
    // matr(matr2);
    Matrix<int> operator() (const Matrix& mat)
    {
        (*this) = mat;
        return *this;
    }


    void swapRows(Matrix<T>& matrix, std::size_t row1, std::size_t row2) {
        std::swap(matrix[row1], matrix[row2]);
    }
    Matrix<T> transpose() const
    {
        Matrix<T> new_matrix(_column, _row);
        for (std::size_t i = 0; i < _row; ++i)
        {
            for (std::size_t j = 0; j < _column; ++j)
            {
                new_matrix[j][i] = _matrix[i][j];
            }
        }
        return new_matrix;
    }
    //stupencheti
    Matrix<T> upperTriangularForm() {
        for (std::size_t i = 0; i < _row; ++i)
        {



        }
    }
    T determinant()//tox poxeluc nshani poxum
    {
        try {
            if (_row == _column)
            {
                T determinant = 1;
                Matrix<T> new_matrix(_row, _column);
                new_matrix._matrix = _matrix;
                new_matrix = new_matrix.upperTriangularForm();
                for (std::size_t i = 0; i < _row; ++i)
                {
                    determinant *= new_matrix[i][i];
                }
                return determinant;
            }
            else
            {
                throw std::logic_error("Matrix is not square for calculating the determinant :(");
            }
        }
        catch (const std::exception& e)
        {
            std::cout << e.what() << '\n';
        }
    }
    T sled() {
        try {
            if (_column == _row)
            {
                T sled = 0;
                for (std::size_t i = 0; i < _row; ++i)
                {
                    sled += _matrix[i][i];
                }
                return sled;
            }
            else {
                throw std::logic_error("Matrix is not square :(");
            }
        }
        catch (const std::exception& e)
        {
            std::cout << e.what() << '\n';
        }
    }
    Matrix<T> inverse()  {

        try {
            if (_row != _column) {
                throw std::logic_error("Incompatible matrix dimensions for add");
            }
            else {
                // Create an augmented matrix [A|I]
                Matrix<T> augmented(_row, _column * 2);
                for (std::size_t i = 0; i < _row; ++i) {
                    for (std::size_t j = 0; j < _column; ++j) {
                        augmented[i][j] = _matrix[i][j];
                    }
                    augmented[i][_column + i] = 1;
                }

                // Apply Gauss-Jordan elimination
                for (std::size_t i = 0; i < _row; ++i) {
                    // Find pivot row
                    std::size_t pivotRow = i;
                    for (std::size_t j = i + 1; j < _row; ++j) {
                        if (std::abs(augmented[j][i]) > std::abs(augmented[pivotRow][i])) {
                            pivotRow = j;
                        }
                    }

                    // Swap rows
                    if (pivotRow != i) {
                        swapRows(augmented, i, pivotRow);
                    }

                    // Scale pivot row
                    T pivotElement = augmented[i][i];
                    for (std::size_t j = i; j < _column * 2; ++j) {
                        augmented[i][j] /= pivotElement;
                    }

                    // Elimination
                    for (std::size_t j = 0; j < _row; ++j) {
                        if (j != i) {
                            T factor = augmented[j][i];
                            for (std::size_t k = i; k < _column * 2; ++k) {
                                augmented[j][k] -= factor * augmented[i][k];
                            }
                        }
                    }
                }

                // Extract the inverted matrix from the augmented matrix
                Matrix<T> inverted(_row, _column);
                for (std::size_t i = 0; i < _row; ++i) {
                    for (std::size_t j = 0; j < _column; ++j) {
                        inverted[i][j] = augmented[i][_column + j];
                    }
                }

                return inverted;
            }

        }
        catch (const std::exception& e) {
            std::cout << e.what() << '\n';
        }
        


    }

    /// Getters
    std::size_t get_row() const
    {
        return _row;
    }
    std::size_t get_column() const
    {
        return _column;
    }
    T get_value(std::size_t i, std::size_t j) const
    {
        return _matrix[i][j];
    }
    std::vector<std::vector<T>> get_matrix() const
    {
        return _matrix;
    }

    /// Setter
    void value_setter(T num, std::size_t index1, std::size_t index2)
    {
        _matrix[index1][index2] = num;
    }

    std::pair<std::size_t, std::size_t> size()
    {
        return std::make_pair(_row, _column);
    }
};

































//------------------------
template <typename T>
bool symmetric(const Matrix<T>& matr)
{
    return matr == matr.transpose();
}
//------------------------
template <typename T>
bool operator==(const Matrix<T>& matr1, const Matrix<T>& matr2)
{

    if (matr1.get_column() != matr2.get_column() || matr1.get_row() != matr2.get_row())
    {
        return false;
    }
    for (std::size_t i = 0; i < matr1.get_row(); ++i)
    {
        for (std::size_t u = 0; u < matr1.get_column(); ++u)
        {
            if (matr1[i][u] != matr2[i][u])
                return false;
        }
    }
    return true;
}
template <typename T>
bool operator!=(const Matrix<T>& matr1, const Matrix<T>& matr2)
{
    return !(matr1 == matr2);
}
//------------------------
template <typename T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& matr)
{
    for (std::size_t i = 0; i < matr.get_row(); ++i)
    {
        for (std::size_t u = 0; u < matr.get_column(); ++u)
        {
            out << matr[i][u] << " ";
        }
        out << '\n';
    }
    return out;
}
template <typename T>
std::istream& operator>>(std::istream& in, Matrix<T>& matr)
{
    T num;
    for (std::size_t i = 0; i < matr.get_row(); ++i)
    {
        for (std::size_t u = 0; u < matr.get_column(); ++u)
        {
            in >> num;
            matr.value_setter(num, i, u);
        }
    }
    return in;
}
//------------------------
