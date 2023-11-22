#include "matrix.h"
#include <chrono>
#include <windows.h>

int main() {
    Matrix matr{
        {1,2,3},
        {1,2,3},
        {1,2,3}
    };
    Matrix matr1 = 0.5 * matr;
    Matrix curr = matr;
    while (curr != matr1)
    {
        std::cout << "No\n";
        curr -= 0.1 * matr;
    }
    std::cout << "Yes";
    return 0;
}
