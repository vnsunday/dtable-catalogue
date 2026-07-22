#include <stdio.h>
#include <stdlib.h>
#include <string>
#include "general.hpp"

void definition_01() {
    int rangeX; // X < x
    int rangeY; // Y < y

    Meta M;
    const char* szX = "X";
    const char* szF = "F";
    M.DefineMapping(szF);
    M.DefineMapping(szX);
}

void definition_02() {
	// Independent
	Meta M;

	const char* szV1 = "X";
	const char* szV2 = "Y";

	M.DefineMapping(szV1, "Real");
	M.DefineMapping(szV2);
	M.DefineMapping(NULL, NULL);
	/*
	M.IfOnlyIf((void*)M.Independent(szV1, szV2),
				(void*)M.Equation("F", "Fx*Fy"));
	*/
}

void definition_03() {
	Meta M;
	M.DefineTable("Distribution", "ValueSet(Y)", "ValueSet(X)");
}

void definition_04() {
    Meta M;
	M.DefineMapping("F", "ValueSet(X)xValueSet(Y)", "Real");
}

void definition_02_01() {
	Meta M;
	M.Equation("M_xy", 0);
	/*
	M.If(
		(void*)M.Equation("M_xy", 0),
		(void*)M.Independent("X", "Y"));
	*/
}

void definition_02_02() {
	const char* szName = "Standard Distribution";
	Meta M;
	Formula F;
	M.DefineMapping("F", "Space(X)xSpace(Y)", "Real");
    
	F.Fraction(0, 0);
	F.Multiply("ExpectationX", "ExpectationY");
	
}

int main(int argc, char **argv) {

	int n = sizeof(long);
	int n1 = sizeof(void*);

	printf("Size(long)=%d; size(Pointer)=%d\r\n", n, n1);
    return 0;
}
