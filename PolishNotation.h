#include <stack>
#include <cctype>
#include <map>
#include <string>
#include "matrix.h"


int getPrecedence(char c) {
    if (c == '+' || c == '-') {
        return 1;
    }
    else if (c == '*') {
        return 2;
    }
    return 0;  // Lower precedence for other characters
}

bool isOperator(char c) {
    return (c == '+' || c == '-' || c == '*');
}

bool isVariable(char c) {
    return (isalpha(c));
}

/// Converting expression to Polish notation
std::string convertToPolish(const std::string& infix) {
    std::stack<char> operators;
    std::stack<std::string> operands;

    for (char token : infix) {
        if (isVariable(token)) {
            operands.push(std::string(1, token));
        }
        else if (isOperator(token)) {
            while (!operators.empty() && getPrecedence(operators.top()) >= getPrecedence(token)) {
                operands.push(std::string(1, operators.top()));
                operators.pop();
            }
            operators.push(token);
        }
        else if (token == '(') {
            operators.push(token);
        }
        else if (token == ')') {
            while (!operators.empty() && operators.top() != '(') {
                operands.push(std::string(1, operators.top()));
                operators.pop();
            }
            operators.pop(); // Pop '('
        }
    }

    while (!operators.empty()) {
        operands.push(std::string(1, operators.top()));
        operators.pop();
    }

    std::string result;
    while (!operands.empty()) {
        result = operands.top() + result;
        operands.pop();
    }

    return result;
}

/// Evaluating the result of Polish notation
Matrix evaluatePolishExpression(const std::string& expression, const std::map<char, Matrix>& variableMap) {
    std::stack<Matrix> operandStack;

    for (char token : expression) {
        if (isVariable(token)) {
            operandStack.push(variableMap.at(token));
        }
        else {
            Matrix operand2 = operandStack.top();
            operandStack.pop();

            Matrix operand1 = operandStack.top();
            operandStack.pop();

            if (token == '+') {
                // Perform matrix addition (assuming matrices have the same size)
                Matrix result = operand1 + operand2;
                operandStack.push(result);
            }
            else if (token == '-') {
                // Perform matrix subtraction (assuming matrices have the same size)
                Matrix result = operand1 - operand2;
                operandStack.push(result);
            }
            else if (token == '*') {
                // Perform matrix multiplication
                Matrix result = operand1 * operand2;
                operandStack.push(result);
            }
        }
    }

    return operandStack.top();
}

std::stack<Matrix> operandStack;

void performOperation(char op) {
    if (operandStack.size() < 2) {
        std::cerr << "Error: Insufficient operands for operator '" << op << "'." << std::endl;
        exit(EXIT_FAILURE);
    }

    Matrix operand2 = operandStack.top();
    operandStack.pop();

    Matrix operand1 = operandStack.top();
    operandStack.pop();

    Matrix result;

    if (op == '+') {
        // Implement matrix addition
        // This is a simple example, replace it with your actual matrix addition code
        result = operand1 + operand2;
        
    }
    else if (op == '*') {
        // Implement matrix multiplication
        // This is a simple example, replace it with your actual matrix multiplication code
        result = operand1 * operand2;
    }
    else {
        std::cerr << "Error: Unsupported operator '" << op << "'." << std::endl;
        exit(EXIT_FAILURE);
    }

    operandStack.push(result);
}