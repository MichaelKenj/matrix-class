#pragma once
#include "matrix.h"
#include "printHelperMethods.h"
#include "PolishNotation.h"
#include <sstream>
#include <cstdlib>


void solveSOLE()
{
	//read about lucumneri qanak of SOLE

	//matr = upper
}

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
	
	//converting to Polish notation
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

	//Matrices 
	std::vector<Matrix> matrices(vars.size());

	//filling map with Variable names and their matrices
	for (std::size_t i = 0; i < matrices.size(); ++i)
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
	//for (auto it = variables.begin(); it != variables.end(); ++it)
		//std::cout << it -> first << ":\n " << it -> second;

	/// After evaluating waiting for next instruction
	std::cout << "Your result is: \n" << result;
	char op = _getche();
	std::cout << "Press q to continue...";
	if (op == 'q') return;
}