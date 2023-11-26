#pragma once
#include "matrix.h"
#include "printHelperMethods.h"
#include "PolishNotation.h"
#include <sstream>
#include <cstdlib>

/// Solving System of Linear equations
void solveSOLE()
{
	system("cls");
	SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), { 0,0 });
	
	std::cout << "Enter variables count: ";
	std::size_t row, col;
	std::cin >> col;
	std::cout << "Enter systems count: ";
	std::cin >> row;

	Matrix matr(row, col);
	std::vector<double> freeTerm(row);

	/// decorative moments
	for (std::size_t i = 0; i < matr.size().first; ++i)
	{
		std::cout << "Equation " << i + 1 << ": \n";
		for (std::size_t u = 0; u < matr.size().second; ++u)
		{
			std::cout << "           Coefficent x[" << u + 1 << "]: ";
			std::cin >> matr[i][u];
		}
	}
	for (std::size_t i = 0; i < freeTerm.size(); ++i)
	{
		std::cout << "FreeTerm[" << i+1 << "]: ";
		std::cin >> freeTerm[i];
	}

	// can i indecate solve quantity?

	/// Create an augmented matrix [A | B]
	Matrix augmentedMatrix(matr.get_row_size(), matr.get_column_size() + 1);
	for (size_t i = 0; i < matr.get_row_size(); ++i) {
		for (size_t j = 0; j < matr.get_column_size(); ++j) {
			augmentedMatrix[i][j] = matr[i][j];
		}
		augmentedMatrix[i][matr.get_column_size()] = freeTerm[i];
	}

	/// Apply Gaussian elimination
	for (size_t i = 0; i < matr.get_row_size(); ++i) {
		/// Make the diagonal element 1
		double divisor = augmentedMatrix[i][i];
		for (size_t j = i; j <= matr.get_column_size(); ++j) {
			augmentedMatrix[i][j] /= divisor;
		}

		/// Make the other elements in the column 0
		for (size_t k = 0; k < matr.get_row_size(); ++k) {
			if (k != i) {
				double factor = augmentedMatrix[k][i];
				for (size_t j = i; j <= matr.get_column_size(); ++j) {
					augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
				}
			}
		}
	}

	/// Extract the solution vector from the augmented matrix
	std::vector<double> solution(matr.get_row_size());
	for (size_t i = 0; i < matr.get_row_size(); ++i) {
		solution[i] = augmentedMatrix[i][matr.get_column_size()];
	}

	for (std::size_t i = 0; i < solution.size(); ++i)
	{
		std::cout << "x[" << i << "]: " << solution[i] << '\n';
	}
	

	/// read about lucumneri qanak of SOLE
	std::cout << "\nPress q to continue...";
	char op = _getche();
	if (op == 'q') return;
}

/// Solving Matrix equation A * X * C = D
void solveME()
{
	//stringov
}

/// Solving simple statements such as A - B * C
void matrixCalc()
{
	system("cls");
	SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), { 0,0 });
	std::cout << "Enter your expression: ";

	//our statement
	std::string _statement;
	std::getline(std::cin, _statement);
	
	//converting into Polish notation
	std::string polish_statement = convertToPolish(_statement);
	
	//Variable name: its matrix
	std::map<char, Matrix> variables;
	
	//Variable names
	std::vector<char> vars;
	
	//filling vars vector by names
	for (auto i : polish_statement)
	{
		if (i != '-' &&
			i != '+' &&
			i != '*')
		{
			vars.push_back(i);
		}
	}

	//filling map with Variable names and their matrices
	for (std::size_t i = 0; i < vars.size(); ++i)
	{
		std::size_t row, col;
		std::cout << vars[i] << ": ";
		std::cout << "Input " << vars[i] << "'s row: ";
		std::cin >> row;
		std::cout << "Input " << vars[i] << "'s column: ";
		std::cin >> col;

		Matrix temp(row, col);
		for (std::size_t u = 0; u < row; ++u)
		{
			for (std::size_t j = 0; j < col; ++j)
			{
				std::cout << vars[i] << '[' << u << ']' 
					<< '[' << j << ']' << ": ";
				std::cin >> temp[u][j];
			}
		}
		variables.insert({ vars[i], temp });
	}
	
	Matrix result = evaluatePolishExpression(polish_statement, variables);
	variables.clear();

	/// After evaluating waiting for next instruction
	std::cout << "Your result is: \n" << result;
	std::cout << "Press q to continue...";
	char op = _getche();
	if (op == 'q') return;
}
