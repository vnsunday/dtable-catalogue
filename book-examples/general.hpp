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
	long Equation(const char* szLeft, const char* szRight);
	std::string Name(const char* szVariable) { return ""; }

	void DefineTable(const char* szName, const char* szRowSet=NULL, const char* szColumnSet=NULL) {}
	void Integral() {}
    void If(void* Left, void* Right);
};

class Formula {
public:
    int Space(long variable) { return 0; }
    int Fraction(long numerator, long denominator) { return 0; }
    int Integral() { return 0; }
};

#endif
