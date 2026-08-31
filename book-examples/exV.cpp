#include <stdio.h>
#include <stdlib.h>
#include "general.hpp"

void ex_1_1() {
    Meta M;
    const char* szM[][2] = {
        { "Null", "" }, 
        { "Alternative", "" }, 
        { "Power", "" } 
    };
    M.Graph("Model", {});
    M.Graph("Substitution", {});
}

void ex_2_1() {
    Meta M;
    const char* szIS[][2] = {
        { "Average", "8500" },
        { "ActualAverage", "8900" },
        { "Distribution", "Normal" },
        { "sigma", "2600" },
        { "SignificanceLevel", "0.05"}
    };
    M.Graph("Model", szIS);
}

void ex_2_5() {
    Meta M;
    // M.Graph("", );
    /*
    const char* szM[][2] = {
        { "", "" }, 
        { "", "" }
    };
    */
}

int main() {
    return 0;
}