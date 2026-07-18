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

	M.DefineMapping(szV1);
	M.DefineMapping(szV2);

	M.IfOnlyIf((void*)M.Independent(szV1, szV2),
				(void*)M.Equation("F", "Fx*Fy"));
}

void definition_03() {
	Meta M;
	M.DefineTable("Distribution");
}

void definition_04() {

}

int main(int argc, char **argv) {
    return 0;
}
