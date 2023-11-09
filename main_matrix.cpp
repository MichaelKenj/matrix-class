#include "matrix.h"

int main() { 
    Matrix<double> matr{
        {1,2,3},
        {1,2,3},
        {1,2,3}
    };
    Matrix<int> matr2{
        {1,2,3},
        {1,2,3},
        {1,2,3},
        {1,2,3}
    };
    matr += matr2;
    
   
    
    
    
}


//Matrix<T> upperTriangularForm() {
//    Matrix<T> new_matrix(_m, _n);
//    new_matrix.matr = matr;
//    std::size_t numRows = _m;
//    std::size_t numCols = _n;
//    for (std::size_t r = 0; r < numRows; ++r) {
//        std::size_t maxRow = r;
//        for (std::size_t i = r + 1; i < numRows; ++i) {
//            if (new_matrix[i][r] > new_matrix[maxRow][r])
//                maxRow = i;
//        }
//
//        if (new_matrix[maxRow][r] == 0)
//            continue;
//
//        swapRows(new_matrix, r, maxRow);
//
//        for (std::size_t i = r + 1; i < numRows; ++i) {
//            T factor_r_next = new_matrix[r][r];
//            T factor_r = new_matrix[i][r];
//            for (std::size_t j = r; j < numCols; ++j)
//                new_matrix[i][j] = new_matrix[i][j] * factor_r_next - new_matrix[r][j] * factor_r;
//        }
//
//    }
//    return new_matrix;
//}










