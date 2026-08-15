#ifndef _DTABLE_CATALOGUE_BOOKEXAMPLE_GENERAL_H_
#define _DTABLE_CATALOGUE_BOOKEXAMPLE_GENERAL_H_

#include <string>

int combination(int k, int n) {
    int ret = 1;
    for (int i=k+1; i<=n;i++) {
        ret = ret*i;
    }
    for (int i=1; i<=k; i++) {
        ret = ret / i;
    }
    return ret;
}

class Meta {
public:
	static void IfOnlyIf(void* Left, void* Right) {}

	void DefineMapping(const char* szName, const char* szFromSpace=NULL, const char* szToSpace=NULL) {}
	void Space(const char* szVariable) {}
	void Independent(const char* Variable1, const char* Variable2) {}
	long Equation(const char* szLeft, const char* szRight) { return 0l; }
	std::string Name(const char* szVariable) { return ""; }

	void DefineTable(const char* szName, const char* szRowSet=NULL, const char* szColumnSet=NULL) {}
	void Integral() {}
    void If(void* Left, void* Right);

    void Graph(const char* szName, const char* param[][2]) {}
    void Unknown(const char* szName) {}
    void Series(const char* szName, const char* szType=NULL, int startIndex=1) {}

    bool Labeled(char* object, const char* szLabel)  { return true; }
};

// Abstration
class Formula {
public:
    int Space(long variable) { return 0; }
    int Fraction(long numerator, long denominator) { return 0; }
    int Integral(long from, long to, long mapping[], long derivative[]) { return 0; }
    char* Multiply(char* L, char* R) { return 0; }
    char* SquareRoot(char* V) { return 0;}
    char* Power(char* x, char* exp) { return 0;}
};

#endif
