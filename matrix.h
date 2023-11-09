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
#include <limits.h>
#include <limits>

///TODO
/*
    double f = 2.5;
    int a = 4;
    long u = a + f;

    In trace() it prints 0, when matrix is not square

    Iterators

    stupencheti
    rankD
    obratni
    determinant
    СЛУ
    A * X * C = B
*/

template <typename T>
class Matrix
{
private:
    std::size_t _row;
    std::size_t _column;
    std::vector<std::vector<T>> _matrix;
public:
    ///---------ITERATOR---------
    class Iterator
    {
    private:
        typename std::vector<std::vector<T>>::iterator it_row;
        typename std::vector<T>::iterator it_col;
    public:
        using value_type = T;
        using reference = value_type&;
        using pointer = value_type*;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::bidirectional_iterator_tag;

        Iterator& operator++()
        {
            ++it_col;
            if (it_col == it_row->end())
            {
                ++it_row;
                it_col = it_row->begin();
            }
            return *this;
        }

        Iterator operator++(int)
        {
            Iterator temp = *this;
            ++this;
            return temp;
        }
        Iterator& operator--()
        {
            ++it_col;
            if (it_col == it_row->end())
            {
                ++it_row;
                it_col = it_row->begin();
            }
            return *this;
        }

        bool operator!=(const Iterator &other) const
        {
            return it_row != other.row || it_col != other.col;
        }

        T &operator*()
        {
            return *it_col;
        }
    };

    ///---------ConstITERATOR---------
    class ConstIterator
    {
    private:
        typename std::vector<std::vector<T>>::iterator it_row;
        typename std::vector<T>::iterator it_col;
    public:
        using value_type = T;
        using reference = const value_type&;
        using pointer = const value_type*;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::bidirectional_iterator_tag;

        Iterator& operator++()
        {
            ++it_col;
            if (it_col == it_row->end())
            {
                ++it_row;
                it_col = it_row->begin();
            }
            return *this;
        }
        Iterator operator++(int)
        {
            Iterator temp = *this;
            ++this;
            return temp;
        }

        Iterator& operator--()
        {
            ++it_col;
            if (it_col == it_row->end())
            {
                ++it_row;
                it_col = it_row->begin();
            }
            return *this;
        }
        Iterator operator--(int)
        {
            Iterator temp = *this;
            --this;
            return temp;
        }

        bool operator!=(const Iterator &other) const
        {
            return it_row != other.row || it_col != other.col;
        }

        T &operator*() const
        {
            return *it_col;
        }
    };

    ///---------ForITERATORS----------
    Iterator begin() const
    {
        return {_matrix.begin(), _matrix.begin()->begin()};
    }
    Iterator end() const
    {
        auto last = _matrix.end();
        return {last, last->begin()};
    }

    ///---------CONSTRUCTORS---------

    /// Merged default constructor, Matrix(it_row, column) and Matrix(it_row, column, filler)
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
    Matrix(Matrix<T>&& other)
    {
        _row = std::move(other._row);
        _column = std::move(other._column);
        _matrix = std::move(other._matrix);
    }
    
    /// Destructor
    ~Matrix() = default;
    
    /// ---------OPERATORS---------

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
            std::cout << e.what();
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
                throw std::logic_error("Incompatible matrix dimensions to subtract\n");
            }
        }
        catch (const std::exception& e)
        {
            std::cout << e.what();
        }
    }

    /// operator* for matrix * matrix
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
            std::cout << e.what();
        }
    }

    /// operator* to multiply by any number
    template <typename U>
    friend Matrix<T> operator*(U value, const Matrix<T>& matrix) {
        Matrix<T> newMatrix(matrix.get_row(), matrix.get_column());
        for (std::size_t i = 0; i < matrix.get_row(); ++i) {
            for (std::size_t j = 0; j < matrix.get_column(); ++j) {
                newMatrix[i][j] = matrix.get_matrix()[i][j] * static_cast<T>(value);
            }
        }
        return newMatrix;
    }

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

    /// operator+=
    template <typename U>
    Matrix<T>& operator+=(const Matrix<U>& other)
    {
        for(std::size_t i = 0; i < _row; ++i)
        {
            for(std::size_t j = 0; j < _column; ++j)
            {
                _matrix[i][j] += other.get_matrix()[i][j];
            }
        }
        return *this;
    }

    /// operator-=
    template <typename U>
    Matrix<double>& operator-=(const Matrix<U>& other)
    {
        for(std::size_t i = 0; i < _row; ++i)
        {
            for(std::size_t j = 0; j < _column; ++j)
            {
                _matrix[i][j] -= other.get_matrix()[i][j];
            }
        }
        return *this;
    }

    /// operator*= for Matrix *= matrix;
    template <typename U>
    Matrix<double>& operator*=(const Matrix<U>& other)
    {
        _matrix = _matrix * other;
        return *this;
    }

    /// operator*= for Matrix *= any type of number
    template <typename U>
    Matrix<double>& operator*=(U num)
    {
        (*this) = (*this) * num;
        return *this;
    }

    /// operator==
    bool operator==(const Matrix<T>& rm) const
    {
        return _matrix == rm._matrix;
    }

    /// operator!=
    bool operator!=(const Matrix<T>& rm) const
    {
        return !(*this == rm);
    }

    /// ---------METHODS---------

    /// Transpose
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

    /// Upper Triangular Form
    Matrix<T> toUpper() const
    {
        for (std::size_t i = 0; i < _row; ++i)
        {



        }
    }

    /// Inverse of matrix
    Matrix<T> inverse() const
    {
        try {
            if (_row != _column) {
                throw std::logic_error("Matrix is not square to calculate the inverse\n");
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
                    // Find pivot it_row
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

                    // Scale pivot it_row
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
            std::cout << e.what();
        }
    }

    /// Determinant
    T determinant()
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
                throw std::logic_error("Matrix is not square for calculating the determinant\n");
            }
        }
        catch (const std::exception& e)
        {
            std::cout << e.what();
        }
    }

    /// Trace
    T trace() const
    {
        T trace = 0;
        try {
            if (_column == _row)
            {
                for (std::size_t i = 0; i < _row; ++i)
                {
                    trace += _matrix[i][i];
                }
                return trace;
            }
            else {
                throw std::logic_error("Matrix is not square to calculate the trace\n");
            }
        }
        catch (const std::exception& e)
        {
            std::cout << e.what();
        }
    }

    /// Size
    std::pair<std::size_t, std::size_t> size()
    {
        return std::make_pair(_row, _column);
    }

    /// IsSymmetric
    bool isSymmetric()
    {
        return (*this) == (*this).transpose();
    }

    /// ?
    void swapRows(Matrix<T>& matrix, std::size_t row1, std::size_t row2) {
        std::swap(matrix[row1], matrix[row2]);
    }

    /// ---------GETTERS----------

    /// Row
    std::size_t get_row() const
    {
        return _row;
    }

    /// Column
    std::size_t get_column() const
    {
        return _column;
    }

    /// Returns Matrix[i][j]
    T get_value(std::size_t i, std::size_t j) const
    {
        return _matrix[i][j];
    }

    /// Matrix
    std::vector<std::vector<T>> get_matrix() const
    {
        return _matrix;
    }

    /// ---------SETTERS---------

    /// Value setter
    void value_setter(T num, std::size_t index1, std::size_t index2)
    {
        _matrix[index1][index2] = num;
    }

};

/// operator* to multiply by any number
template <typename T, typename U>
Matrix<T> operator*(const Matrix<T>& matrix, U value) {
    Matrix<T> newMatrix(matrix.get_row(), matrix.get_column());
    for (std::size_t i = 0; i < matrix.get_row(); ++i) {
        for (std::size_t j = 0; j < matrix.get_column(); ++j) {
            newMatrix[i][j] = matrix.get_matrix()[i][j] * static_cast<T>(value);
        }
    }
    return newMatrix;
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& matr)
{

    for (std::size_t i = 0; i < matr.get_row(); ++i)
    {
        (i == 0) ? std::cout << "⎛" : ((i == matr.get_row() - 1) ? std::cout<< "⎝" : std::cout << "│");
        for (std::size_t u = 0; u < matr.get_column(); ++u)
        {
            (u == matr.get_column() - 1) ? (out << matr[i][u]) : (out << matr[i][u] << " ");
        }
        (i == 0) ? std::cout << "⎞" : ((i == matr.get_row() - 1) ? std::cout<< "⎠" : std::cout << "│");
        std::cout << std::endl;
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
