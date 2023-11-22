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
    double f = 2.5;
    int a = 4;
    long u = a + f;

    In trace() it prints 0, when matrix is not square

    Iterator works, but im not sure, that it correct
    ConstIterator
    ReverseIterator
    ConstReverseIterator

    Is it necessary , that return value of operators - + * is Matirx<double>?

    stupencheti
    rankD
    obratni
    determinant
    СЛУ
    A * X * C = B
*/


class Matrix
{
private:
    std::size_t _row;
    std::size_t _column;
    std::vector<std::vector<double>> _matrix;
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
    Matrix(Matrix& other) : 
        _row(other._row)
        ,_column(other._column)
        ,_matrix(other._matrix)
    {
        std::cout << "Copy called\n";
    }

    /// Initializer list constructor
    Matrix(const std::initializer_list<std::initializer_list<double>>& list)
    {
        std::cout << "Init called\n";
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
        std::cout << "Move called\n";
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
            if (_row == other.get_row_size() && _column == other.get_column_size())
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
    template <typename U>
    Matrix operator-(const Matrix& other) const
    {
        try {
            if (_row == other.get_row_size() && _column == other.get_column_size())
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
        Matrix newMatrix(matrix.get_row_size(), matrix.get_column_size());
        for (std::size_t i = 0; i < matrix.get_row_size(); ++i) {
            for (std::size_t j = 0; j < matrix.get_column_size(); ++j) {
                newMatrix[i][j] = matrix.get_value(i,j) * value;
            }
        }
        return newMatrix;
    }

    /// operator= to Matrix
    Matrix& operator=(const Matrix& other) {
        _row = other.get_row_size();
        _column = other.get_column_size();
        _matrix.resize(_row, std::vector<double>(_column));
        for (std::size_t i = 0; i < _row; ++i) {
            for (std::size_t j = 0; j < _column; ++j) {
                _matrix[i][j] = other.get_value(i,j);
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
                _matrix[i][j] += other.get_value(i,j);
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
                _matrix[i][j] -= other.get_value(i,j);
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
    Matrix transpose() const
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

    /// Gauss form
    /*
    std::vector<double> rowTransform(std::vector<double>& row1, std::vector<double>& row2, std::size_t offset) {
        try {
            double a = row1[offset];
            double b = row2[offset];

            if (b == 0) {
                return row2;
            }

            for (std::size_t i = 0; i < _column; ++i) {
                row1[i] *= b;
            }

            for (std::size_t i = 0; i < _column; ++i) {
                row2[i] *= a;
            }

            for (std::size_t i = 0; i < _column; ++i) {
                row2[i] -= row1[i];
            }

            for (std::size_t i = 0; i < _column; ++i) {
                row1[i] /= b;
            }

            return row2;
        }
        catch (const std::exception& e) {
            std::cout << e.what();
            return row2; // Return original row2 in case of an exception
        }
    }
    Matrix normalize(const Matrix& matrix) {
        Matrix normalizedMatrix = matrix;

        try {
            std::size_t rows = matrix.get_row_size();
            std::size_t cols = matrix.get_column_size();

            // Find the maximum element in the matrix
            double maxElement = matrix[0][0];
            for (std::size_t i = 0; i < rows; ++i) {
                for (std::size_t j = 0; j < cols; ++j) {
                    if (matrix[i][j] > maxElement) {
                        maxElement = matrix[i][j];
                    }
                }
            }

            // Normalize the matrix by dividing each element by the maximum element
            for (std::size_t i = 0; i < rows; ++i) {
                for (std::size_t j = 0; j < cols; ++j) {
                    normalizedMatrix.value_setter(matrix[i][j] / maxElement, i, j);
                }
            }
        }
        catch (const std::exception& e) {
            std::cout << e.what();
        }

        return normalizedMatrix;
    }
    Matrix upperTriangularForm(Matrix& mat, const Matrix& b) {
        mat.addColumn(b);

        std::size_t current = 0;

        while (current < mat.get_row_size() - 1) {
            //mat = mat.normalize(current);
            std::size_t nextrow = current + 1;

            while (nextrow < mat.get_row_size()) {
                mat[nextrow] = mat.rowTransform(mat[current], mat[nextrow], current);
                ++nextrow;
            }

            ++current;
        }

        // Reverse substitution part
        Matrix result(mat.get_row_size(), 1, 0.0);
        result[mat.get_row_size() - 1][0] = mat[mat.get_row_size() - 1][mat.get_column_size() - 1] / mat[mat.get_row_size() - 1][mat.get_row_size() - 1];
        std::size_t i = mat.get_row_size() - 1;

        while (i != 0) {
            double rowsum = 0;
            std::size_t j = i + 1;

            while (j < mat.get_row_size()) {
                rowsum += (mat[i - 1][j - 1]) * (result[j - 1][0]);
                ++j;
            }

            result[i - 1][0] = (mat[i - 1][mat.get_column_size() - 1] - rowsum) / mat[i - 1][i - 1];
            --i;
        }

        return result;
    }
    */

    /// Inverse of matrix
    /*
    Matrix inverse() const
    {
        try {
            if (_row != _column) {
                throw std::logic_error("Matrix is not square to calculate the inverse\n");
            }
            else {
                // Create an augmented matrix [A|I]
                Matrix augmented(_row, _column * 2);
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
                    double pivotElement = augmented[i][i];
                    for (std::size_t j = i; j < _column * 2; ++j) {
                        augmented[i][j] /= pivotElement;
                    }

                    // Elimination
                    for (std::size_t j = 0; j < _row; ++j) {
                        if (j != i) {
                            double factor = augmented[j][i];
                            for (std::size_t k = i; k < _column * 2; ++k) {
                                augmented[j][k] -= factor * augmented[i][k];
                            }
                        }
                    }
                }

                // Extract the inverted matrix from the augmented matrix
                Matrix inverted(_row, _column);
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
    */

    /// Determinant
    /*
    double determinant()
    {
        try {
            if (_row == _column)
            {
                double determinant = 1;
                Matrix new_matrix(_row, _column);
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
    */

    /// Trace
    double trace() const
    {
        double trace_ = 0;
        try {
            if (_column == _row)
            {
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

    void addColumn(const Matrix& column) {
        try {
            if (_row != column.get_row_size()) {
                throw std::logic_error("Incompatible matrix dimensions to add a column\n");
            }
            else {
                for (std::size_t i = 0; i < _row; ++i) {
                    _matrix[i].push_back(column[i][0]);
                }
                ++_column;
            }
        }
        catch (const std::exception& e) {
            std::cout << e.what();
        }
    }

    /// ---------GETTERS----------

    /// Row
    std::size_t get_row_size() const
    {
        return _row;
    }

    /// Column
    std::size_t get_column_size() const
    {
        return _column;
    }

    /// Returns Matrix[i][j]
    double get_value(std::size_t i, std::size_t j) const
    {
        return _matrix[i][j];
    }

    /// Matrix
    std::vector<std::vector<double>> get_matrix() const
    {
        return _matrix;
    }

    /// ---------SETTERS---------

    /// Value setter
    void value_setter(double num, std::size_t index1, std::size_t index2)
    {
        _matrix[index1][index2] = num;
    }

};

/// operator* to multiply by any number
Matrix operator*(const Matrix& matrix, double value) {
    Matrix newMatrix(matrix.get_row_size(), matrix.get_column_size());
    for (std::size_t i = 0; i < matrix.get_row_size(); ++i) {
        for (std::size_t j = 0; j < matrix.get_column_size(); ++j) {
            newMatrix[i][j] = matrix.get_value(i,j) * value;
        }
    }
    return newMatrix;
}


std::ostream& operator<<(std::ostream& out, const Matrix& matr)
{

    for (std::size_t i = 0; i < matr.get_row_size(); ++i)
    {
        (i == 0) ? std::cout << "⎛" : ((i == matr.get_row_size() - 1) ? std::cout<< "⎝" : std::cout << "│");
        for (std::size_t u = 0; u < matr.get_column_size(); ++u)
        {
            (u == matr.get_column_size() - 1) ? (out << matr[i][u]) : (out << matr[i][u] << " ");
        }
        (i == 0) ? std::cout << "⎞" : ((i == matr.get_row_size() - 1) ? std::cout<< "⎠" : std::cout << "│");
        std::cout << std::endl;
    }
    return out;
}
std::istream& operator>>(std::istream& in, Matrix& matr)
{
    double num;
    for (std::size_t i = 0; i < matr.get_row_size(); ++i)
    {
        for (std::size_t u = 0; u < matr.get_column_size(); ++u)
        {
            in >> num;
            matr.value_setter(num, i, u);
        }
    }
    return in;
}
