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
#include <limits>


///TODO
/*
    In trace() it prints 0, when matrix is not square

    Iterator works, but im not sure, that it correct
    ConstIterator
    ReverseIterator
    ConstReverseIterator

    stupencheti
    rankD
    СЛУ
    A * X * C = B

    write comments about working of functions 
    to improve them in future
*/

class Matrix
{
private:
    std::size_t _row;
    std::size_t _column;
    std::vector<std::vector<double>> _matrix;

    Matrix createMinorMatrix(size_t rowToRemove, size_t colToRemove) const {
        size_t n = _row;
        size_t m = _column;

        Matrix minor(n - 1, m - 1);

        for (size_t i = 0, new_i = 0; i < n; ++i) {
            if (i == rowToRemove) {
                continue;  // Skip the row to be removed
            }

            for (size_t j = 0, new_j = 0; j < m; ++j) {
                if (j == colToRemove) {
                    continue;  // Skip the column to be removed
                }

                minor[new_i][new_j] = _matrix[i][j];
                ++new_j;
            }

            ++new_i;
        }
        return minor;
    }
public:
    ///---------ITERATOR---------
    class Iterator
    {
    private:
        std::size_t index;
        std::vector<std::vector<double>>& data_;
    public:
        friend Matrix;
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type = double;
        using difference_type = std::ptrdiff_t;
        using pointer = value_type*;
        using reference = value_type&;

        /// Constructors
        Iterator(std::size_t ro, std::vector<std::vector<double>>& data)
                :index(ro)
                ,data_(data)
        {}
        Iterator(const Iterator& other)
                :index(other.index)
                ,data_(other.data_)
        {}

        /// Destructor
        ~Iterator() = default;

        /// Operators
        Iterator& operator=(const Iterator& other)
        {
            index = other.index;
            return *this;
        }

        Iterator& operator++()
        {
            ++index;
            return *this;
        }
        Iterator operator++(int)
        {
            Iterator tmp = *this;
            ++index;
            return tmp;
        }
        Iterator& operator--()
        {
            --index;
            return *this;
        }
        Iterator operator--(int)
        {
            Iterator tmp = *this;
            --index;
            return tmp;
        }
        reference operator*() const
        {
            std::size_t row = index / data_[0].size();
            std::size_t col = index % data_[0].size();
            return data_[row][col];
        }
        friend bool operator==(const Iterator& lhs, const Iterator& rhs)
        {
            return lhs.index == rhs.index && lhs.data_ == rhs.data_;
        }
        friend bool operator!=(const Iterator& lhs, const Iterator& rhs)
        {
            return !(lhs == rhs);
        }

    };
    ///---------ConstITERATOR---------
    class ConstIterator
    {
    private:
        std::size_t index;
        std::vector<std::vector<double>>& data_;
    public:
        friend Matrix;
        using iterator_category = std::bidirectional_iterator_tag;
        using value_type = const double;
        using difference_type = std::ptrdiff_t;
        using pointer = value_type*;
        using reference = const value_type&;

        /// Constructors
        ConstIterator(std::size_t ro, std::vector<std::vector<double>>& data)
                :index(ro)
                ,data_(data)
        {}
        ConstIterator(const ConstIterator& other)
                :index(other.index)
                ,data_(other.data_)
        {}

        /// Destructor
        ~ConstIterator() = default;

        /// Operators
        ConstIterator& operator=(const ConstIterator& other)
        {
            index = other.index;
            return *this;
        }

        ConstIterator& operator++()
        {
            ++index;
            return *this;
        }
        ConstIterator operator++(int)
        {
            ConstIterator tmp = *this;
            ++index;
            return tmp;
        }
        ConstIterator& operator--()
        {
            --index;
            return *this;
        }
        ConstIterator operator--(int)
        {
            ConstIterator tmp = *this;
            --index;
            return tmp;
        }
        reference operator*() const
        {
            std::size_t row = index / data_[0].size();
            std::size_t col = index % data_[0].size();
            return data_[row][col];
        }
        friend bool operator==(const ConstIterator& lhs, const ConstIterator& rhs)
        {
            return lhs.index == rhs.index && lhs.data_ == rhs.data_;
        }
        friend bool operator!=(const ConstIterator& lhs, const ConstIterator& rhs)
        {
            return !(lhs == rhs);
        }

    };
    ///---------ForITERATORS----------
    Iterator begin()
    {
        return Iterator{ 0, _matrix };
    }
    Iterator end()
    {
        return Iterator{ _matrix.size() * _matrix[0].size(), _matrix };
    }
    ConstIterator cbegin()
    {
        return ConstIterator{ 0, _matrix };
    }
    ConstIterator cend()
    {
        return ConstIterator{ _matrix.size() * _matrix[0].size(), _matrix };
    }
    Iterator rbegin()
    {
        // Calculate the index for the last element
        std::size_t lastIndex = _matrix.size() * _matrix[0].size() - 1;

        // Return an iterator pointing to the last element
        return Iterator(lastIndex, _matrix);
    }
    Iterator rend()
    {
        // Return an iterator pointing to one position before the first element
        return Iterator(-1, _matrix);
    }
    ///---------CONSTRUCTORS---------

    /// Merged default constructor, Matrix(it_row, column) and Matrix(it_row, column, filler)
    /// (rows, columns, value) constructor, filling by value
    Matrix(std::size_t row = 0, std::size_t column = 0, std::size_t filler = 0) : 
        _row(row)
        , _column(column)
    {
        _matrix.resize(_row, std::vector<double>(_column, filler));
    }

    /// Copy constructor
    Matrix(const Matrix& other) : 
        _row(other._row)
        ,_column(other._column)
        ,_matrix(other._matrix)
    {}

    /// Initializer list constructor
    Matrix(const std::initializer_list<std::initializer_list<double>>& list)
    {
        _row = list.size();
        if (_row > 0)
            _column = list.begin()->size();
        else
            _column = 0;
        _matrix.resize(_row, std::vector<double>(_column));
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
    Matrix( Matrix&& other) noexcept
    {
        _row = std::move(other._row);
        _column = std::move(other._column);
        _matrix = std::move(other._matrix);
    }
    
    /// Destructor
    ~Matrix() = default;
    
    /// ---------OPERATORS---------

    /// operator[]
    std::vector<double>& operator[](std::size_t index)
    {
        return _matrix[index];
    }
    const std::vector<double> operator[](std::size_t index) const
    {
        return _matrix[index];
    }

    /// operator+
    Matrix operator+(const Matrix& other) const
    {
        try {
            if (_row == other.row() && _column == other.column())
            {
                Matrix new_matrix(_row, _column);
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
    Matrix operator-(const Matrix& other) const
    {
        try {
            if (_row == other.row() && _column == other.column())
            {
                Matrix new_matrix(_row, _column);
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
    Matrix operator*(const Matrix& other) const
    {
        try {
            if (_column != other._row)
            {
                throw std::logic_error("Incompatible matrix dimensions to multiplicate\n");
            }
            else {
                Matrix new_matrix(_row, other._column);

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
    friend Matrix operator*(double value, const Matrix& matrix) {
        Matrix newMatrix(matrix.row(), matrix.column());
        for (std::size_t i = 0; i < matrix.row(); ++i) {
            for (std::size_t j = 0; j < matrix.column(); ++j) {
                newMatrix[i][j] = matrix[i][j] * value;
            }
        }
        return newMatrix;
    }

    /// operator= to Matrix
    Matrix& operator=(const Matrix& other) {
        _row = other.row();
        _column = other.column();
        _matrix.resize(_row, std::vector<double>(_column));
        for (std::size_t i = 0; i < _row; ++i) {
            for (std::size_t j = 0; j < _column; ++j) {
                _matrix[i][j] = other[i][j];
            }
        }
        return *this;
    }

    /// operator+=
    Matrix& operator+=(const Matrix& other)
    {
        for(std::size_t i = 0; i < _row; ++i)
        {
            for(std::size_t j = 0; j < _column; ++j)
            {
                _matrix[i][j] += other[i][j];
            }
        }
        return *this;
    }

    /// operator-=
    Matrix& operator-=(const Matrix& other)
    {
        for(std::size_t i = 0; i < _row; ++i)
        {
            for(std::size_t j = 0; j < _column; ++j)
            {
                _matrix[i][j] -= other[i][j];
            }
        }
        return *this;
    }

    /// operator*= for Matrix *= matrix;
    Matrix& operator*=(const Matrix& other)
    {
        (*this) = (*this) * other;
        return *this;
    }

    /// operator*= for Matrix *= any type of number
    Matrix& operator*=(double num)
    {
        (*this) = (*this) * num;
        return *this;
    }

    /// operator==
    bool operator==(const Matrix& rm) const {
        if (_row != rm._row || _column != rm._column) {
            return false;
        }

        for (std::size_t i = 0; i < _row; ++i) {
            for (std::size_t j = 0; j < _column; ++j) {
                if (std::abs(_matrix[i][j] - rm._matrix[i][j]) > std::numeric_limits<double>::epsilon() * 100) {
                    return false;
                }
            }
        }

        return true;
    }

    /// operator!=
    bool operator!=(const Matrix& rm) const
    {
        return !(*(this) == rm);
    }

    /// ---------METHODS---------

    /// Transpose
    const Matrix transpose() const
    {
        Matrix new_matrix(_column, _row);
        for (std::size_t i = 0; i < _row; ++i)
        {
            for (std::size_t j = 0; j < _column; ++j)
            {
                new_matrix[j][i] = _matrix[i][j];
            }
        }
        return new_matrix;
    }
    
    /// Inverse of matrix
    Matrix inverse() const {
        try {
            if (_row == _column && determinant() != 0)
            {
                Matrix augmentedMatrix = *this;

                std::size_t n = _row;
                std::size_t m = _column;

                // Augment the matrix with the identity matrix
                for (std::size_t i = 0; i < n; ++i) {
                    augmentedMatrix[i].resize(2 * n, 0.0);
                    augmentedMatrix[i][n + i] = 1.0;
                }

                // Perform Gaussian elimination
                for (std::size_t i = 0; i < n; ++i) {
                    // Find the pivot row
                    std::size_t pivotRow = i;
                    for (std::size_t k = i + 1; k < n; ++k) {
                        if (std::fabs(augmentedMatrix[k][i]) > std::fabs(augmentedMatrix[pivotRow][i])) {
                            pivotRow = k;
                        }
                    }

                    // Swap rows if necessary
                    if (pivotRow != i) {
                        std::swap(augmentedMatrix[i], augmentedMatrix[pivotRow]);
                    }

                    // Make the diagonal element 1
                    double divisor = augmentedMatrix[i][i];
                    for (std::size_t j = i; j < 2 * n; ++j) {
                        augmentedMatrix[i][j] /= divisor;
                    }

                    // Make the other elements in the column 0
                    for (std::size_t k = 0; k < n; ++k) {
                        if (k != i) {
                            double factor = augmentedMatrix[k][i];
                            for (std::size_t j = i; j < 2 * n; ++j) {
                                augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                            }
                        }
                    }
                }

                // Extract the inverse matrix from the augmented matrix
                Matrix inverseMatrix(n, n);
                for (std::size_t i = 0; i < n; ++i) {
                    for (std::size_t j = 0; j < n; ++j) {
                        inverseMatrix[i][j] = augmentedMatrix[i][n + j];
                    }
                }

                return inverseMatrix;
            }
            else
            {
                throw std::logic_error("Matrix is not square to calculate the inverse\n");
            }
        }
        catch (const std::exception& e)
        {
            std::cout << e.what();
        }
    }

    /// Determinant - WORKS PROPERLY
    double determinant() const {
        if (_row != _column) {
            throw std::invalid_argument("Matrix must be square to calculate the determinant.");
        }

        std::size_t n = _row;

        // Base case: a 2x2 matrix
        if (n == 2) {
            return _matrix[0][0] * _matrix[1][1] - _matrix[0][1] * _matrix[1][0];
        }

        double det = 0.0;

        for (std::size_t j = 0; j < n; ++j) {
            // Calculate the minor matrix without the current row and column
            Matrix minor = createMinorMatrix(0, j);

            // Calculate the determinant using Laplace expansion
            det += (j % 2 == 0 ? 1 : -1) * _matrix[0][j] * minor.determinant();
        }

        return det;
    }

    /// Trace
    double trace() const
    {
        
        try {
            if (_column == _row)
            {
                double trace_ = 0;
                for (std::size_t i = 0; i < _row; ++i)
                {
                    trace_ += _matrix[i][i];
                }
                return trace_;
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
    
    /// Row - WORKS PROPERLY
    std::size_t row() const
    {
        return _row;
    }

    /// Column - WORKS PROPERLY
    std::size_t column() const
    {
        return _column;
    }


    /// IsSymmetric - WORKS PROPERLY
    bool isSymmetric() const 
    {
        return *this == transpose();
    }
};

/// operator* to multiply by any number
Matrix operator*(const Matrix& matrix, double value) {
    Matrix newMatrix(matrix.row(), matrix.column());
    for (std::size_t i = 0; i < matrix.row(); ++i) {
        for (std::size_t j = 0; j < matrix.column(); ++j) {
            newMatrix[i][j] = matrix[i][j] * value;
        }
    }
    return newMatrix;
}

std::ostream& operator<<(std::ostream& out, const Matrix& matr)
{

    for (std::size_t i = 0; i < matr.row(); ++i)
    {
        (i == 0) ? std::cout << "|" : ((i == matr.row() - 1) ? std::cout<< "|" : std::cout << "|");
        for (std::size_t u = 0; u < matr.column(); ++u)
        {
            (u == matr.column() - 1) ? (out << matr[i][u]) : (out << matr[i][u] << " ");
        }
        (i == 0) ? std::cout << "|" : ((i == matr.row() - 1) ? std::cout<< "|" : std::cout << "|");
        std::cout << std::endl;
    }
    return out;
}
std::istream& operator>>(std::istream& in, Matrix& matr)
{
    double num;
    for (std::size_t i = 0; i < matr.row(); ++i)
    {
        for (std::size_t j = 0; j < matr.column(); ++j)
        {
            in >> num;
            matr[i][j] = num;
        }
    }
    return in;
}
