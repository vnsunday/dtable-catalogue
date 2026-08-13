#include <stdio.h>
#include <stdlib.h>
#include "general.hpp"

void ex3_1() {
    Meta M;
    const char* A[][2] = {
        { "Distribution", "A Table"}
    };

    M.Graph("X", A);
    M.Equation("Z", "X^2");
    M.Unknown("Z");
}

void ex3_3() {
    Meta M;
    M.DefineMapping("X");
    M.DefineMapping("Y");
    M.Unknown("X+Y");
}

int main(int argc, char const *argv[]) 
{
    ex3_1();
    ex3_3();
    return 0;
}
