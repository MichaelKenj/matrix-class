#define NOMINMAX
#include <chrono>
#include <thread>
#include <Windows.h>
#include <conio.h>  
#include <dos.h>   
#include <windows.h>  
#include <time.h>  
#include <limits>
#include "solve.h"

///-----------PRINT FRAMES OF TITLES------------

// _hell is a vector of rows of big title
void printFrame(const std::vector<std::string>& _hell)
{
    std::cout << "\n\n\n\n\n\n\n\n\n\n";
    for (int i = 0; i < _hell.size(); ++i)
    {
        std::cout << _hell[i] << '\n';
    }

}
///-----------FOR HELLO------------

// Sliding Hello
void slideHello()
{
    std::string row1 = "±±±±    ±±±±   ±±±±±±±±±±±±   ±±±±           ±±±±              ±±±±±±";
    std::string row2 = "±±±±    ±±±±   ±±±±           ±±±±           ±±±±            ±±±±  ±±±±";
    std::string row3 = "±±±±    ±±±±   ±±±±           ±±±±           ±±±±           ±±±±    ±±±±";
    std::string row4 = "±±±±±±±±±±±±   ±±±±±±±±±±±±   ±±±±           ±±±±           ±±±±    ±±±±";
    std::string row5 = "±±±±±±±±±±±±   ±±±±±±±±±±±±   ±±±±           ±±±±           ±±±±    ±±±±";
    std::string row6 = "±±±±    ±±±±   ±±±±           ±±±±           ±±±±           ±±±±    ±±±±";
    std::string row7 = "±±±±    ±±±±   ±±±±           ±±±±           ±±±±            ±±±±  ±±±±";
    std::string row8 = "±±±±    ±±±±   ±±±±±±±±±±±±   ±±±±±±±±±±±±   ±±±±±±±±±±±±      ±±±±±±";
    
    std::vector<std::string> _hello{ row1, row2, row3, row4, row5, row6, row7, row8 };
    int spaces = 24;

    std::string sps = "";
    for (int i = 0; i < 50; ++i)
        sps += " ";
    for (std::size_t i = 0; i < _hello.size(); ++i)
    {
        _hello[i] = sps + _hello[i];
    }
    
    while (spaces >= 0) 
    {
        printFrame(_hello);
        std::this_thread::sleep_for(std::chrono::milliseconds(30));
        
        --spaces;
        if (spaces == 0)
            std::this_thread::sleep_for(std::chrono::milliseconds(2000));

        system("cls");
        SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), { 0,0 });
        for (std::size_t i = 0; i < _hello.size(); ++i)
        {
            _hello[i].erase(0, 1);
        }
    }
}
///-----------FOR MENU-------------

// Sliding Matrix
void slideMatrixInc()
{
    std::string row1 = "±±±±       ±±±±       ±±±±±      ±±±±±±±±±±±±±±   ±±±±±±±±±±      ±±±±   ±±±±   ±±±±";
    std::string row2 = "±±±±±     ±±±±±     ±±±± ±±±±    ±±±±±±±±±±±±±±   ±±±±±±±±±±±±    ±±±±   ±±±±   ±±±±";
    std::string row3 = "±±±±±±   ±±±±±±    ±±±±   ±±±±        ±±±±        ±±±±     ±±±±           ±±±± ±±±±";
    std::string row4 = "±±±±±±± ±±±±±±±    ±±±±   ±±±±        ±±±±        ±±±±     ±±±    ±±±±     ±±±±±±±";
    std::string row5 = "±±±± ±±±±± ±±±±   ±±±±     ±±±±       ±±±±        ±±±±   ±±±±±    ±±±±      ±±±±±";
    std::string row6 = "±±±±  ±±±  ±±±±   ±±±±±±±±±±±±±       ±±±±        ±±±±±±±±±       ±±±±     ±±±±±±±";
    std::string row7 = "±±±±       ±±±±   ±±±±     ±±±±       ±±±±        ±±±±   ±±±±±    ±±±±    ±±±± ±±±±";
    std::string row8 = "±±±±       ±±±±   ±±±±     ±±±±       ±±±±        ±±±±    ±±±±±   ±±±±   ±±±±   ±±±±";

    std::vector<std::string> _matrix{ row1, row2, row3, row4, row5, row6, row7, row8 };
    int spaces = 28;

    std::string sps = "";
    for (int i = 0; i < 50; ++i)
        sps += " ";
    for (std::size_t i = 0; i < _matrix.size(); ++i)
    {
        _matrix[i] = sps + _matrix[i];
    }

    while (spaces >= 0)
    {
        printFrame(_matrix);
        std::this_thread::sleep_for(std::chrono::milliseconds(30));

        --spaces;
        if (spaces == 0)
            std::this_thread::sleep_for(std::chrono::milliseconds(2000));

        system("cls");
        SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), { 0,0 });
        for (std::size_t i = 0; i < _matrix.size(); ++i)
        {
            _matrix[i].erase(0, 1);
        }
    }
}
///-----------START WORKING------------



/// Printing the main menu
void printInterface();

/// Printing Solving System of Linear Equations
void printSolveSOLEInterface()
{
    setcursor(0, 0);
    srand((unsigned)time(NULL));

    do {
        system("cls");
        gotoxy(10, 5); std::cout << " --------------------------- ";
        gotoxy(10, 6); std::cout << " |System of linear equation| ";
        gotoxy(10, 7); std::cout << " ---------------------------";
        gotoxy(10, 9); std::cout << "1. Input your system ";
        gotoxy(10, 10); std::cout << "2. Return to main menu ";
        gotoxy(10, 12); std::cout << "Select an option: ";

        char op = _getche();

        if (op == '2')
        {
            printInterface();
        }
        else if (op == '1')
        {
            solveSOLE();
            // Clear the input buffer
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    } while (true);
}

/// Printing the Solving Matrix equation
void printSolveMEInterface()
{
    setcursor(0, 0);
    srand((unsigned)time(NULL));

    do {
        system("cls");
        gotoxy(10, 5); std::cout << " --------------------------- ";
        gotoxy(10, 6); std::cout << " | Solving matrix equation | ";
        gotoxy(10, 7); std::cout << " --------------------------- ";
        gotoxy(10, 9); std::cout << "1. Input your equation ";
        gotoxy(10, 10); std::cout << "2. Return to main menu ";
        gotoxy(10, 12); std::cout << "Select an option: ";

        char op = _getche();

        if (op == '2')
        {
            printInterface();
        }
        else if (op == '1')
        {
            solveME();
            // Clear the input buffer
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    } while (true);
}

/// Printing Matrix calculation
void printMatrixCalcInterface()
{
    setcursor(0, 0);
    srand((unsigned)time(NULL));

    do {
        system("cls");
        gotoxy(10, 5); std::cout << " -------------------------- ";
        gotoxy(10, 6); std::cout << " |    Matrix calculator   | ";
        gotoxy(10, 7); std::cout << " --------------------------";
        gotoxy(10, 9); std::cout << "1. Input your statement ";
        gotoxy(10, 10); std::cout << "2. Return to main menu ";
        gotoxy(10, 12); std::cout << "Select an option: ";
        
        char op = _getche();
        

        if (op == '2')
        {
            printInterface();
        }
        else if(op == '1')
        {
            matrixCalc();
            // Clear the input buffer
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
    } while (true);
}

/// Main interface
void printInterface()
{
    setcursor(0, 0);
    srand((unsigned)time(NULL));

    do {
        system("cls");
        gotoxy(10, 5); std::cout << " -------------------------- ";
        gotoxy(10, 6); std::cout << " |       Matrix Inc.      | ";
        gotoxy(10, 7); std::cout << " --------------------------";
        gotoxy(10, 9); std::cout << "1. Solve a system of linear equations";
        gotoxy(10, 10); std::cout << "2. Solve matrix equation";
        gotoxy(10, 11); std::cout << "3. Matrix calculator";
        gotoxy(10, 12); std::cout << "4. Quit";
        gotoxy(10, 14); std::cout << "Select option: ";
        char op = _getche();
        if (op == '1')
            printSolveSOLEInterface();
        else if (op == '2') 
            printSolveMEInterface();
        else if (op == '3') 
            printMatrixCalcInterface();
        else if (op == '4') 
            exit(0);

    } while (1);
}
