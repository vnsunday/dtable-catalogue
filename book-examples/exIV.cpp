#include <stdio.h>
#include <stdlib.h>
#include "general.hpp"

void ex_$1_1_1() {
    int field_id1 = 1; // Height 
    int field_id2 = 2; // Count
}

void ex_2_2() {
    // Distribution Table
    Meta M;
    M.DefineTable("Distribution", "p(x)", "x");
    M.DefineMapping("X1");
    M.DefineMapping("X2");
    M.Unknown("X1+X2");
}

void ex_2_3() {
    Meta M;
    M.DefineTable("WeightDistribution");
    M.Unknown("Features");
}

void ex_3_1() {
    Meta M;
    const char* sz[][2] = {
        { "Distribution", "Normal" }, 
        { "sigma", "c" },
        { "Expctation", "alpha" },
    };
    const char* szU[][2] = {
        { "AverageSample", "X" } 
    };
    const char* szP[][2] = {
        { "Equals", "Expectation; AverageSample" }
    };
    M.DefineMapping("X", "Variables");
    M.Graph("Assumption", sz );
    M.Graph("Proof", szP);
}

void ex_3_3() {
    Meta M;
    const char* szAS[][2] = {
        { "Distribution", "Possion" }, 
        { "Paramter", "Lambda" }
    };
    const char* szUK[][2] = {
        { "ParameterEstimation", "Lambda" }
    };

    M.Graph("Assumption", szAS);
    M.Graph("Unknown", szUK);
}

void ex_4_1() {
    Meta M;
    const char* szPre[][2] = {
        { "Name", "Normal" }, 
        { "Distribution", "Normal" }, 
        { "Sigma", "0.3" }    
    };
    const char* szUK[][2] = {
        { "ConfidenceInterval", "95%" }
    };

    M.Graph("Known", szPre);
    M.Graph("Unknown", szUK);
}

void ex4_5() {
    Meta M;
    const char* szTest[][2] = {
        { "Count", "600" }, 
        { "Significance", "0.95" }
    };
}

int main()
{
    return 0;
}