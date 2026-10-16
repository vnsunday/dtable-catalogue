#include <stdio.h>
#include <stdlib.h>

double GetE() {
    double dE = 0.0;
    int nLimit = 15;

    double dFactorial = 1.0;
    int i = 0;
    while (i <= nLimit) {
        dE += 1.0/(dFactorial);
        i++;
        dFactorial *= i;
    }

    return dE;
}

double getPI() {
    double dPI = 0.0;
    double i = 1.0;
    int j = 0;
    int nLimit = 10000;

    while ( i <= nLimit) {
        dPI = (j == 0) ? dPI + 4.0 / i : dPI - 4.0 / i;
        
        j = 1 - j;
        i += 2.0;
    }

    return dPI;
}

int main() {

    double EV = GetE();
    double PV = getPI();
    printf("e=%0.5f\r\n", EV);
    printf("pi=%0.5f\r\n", PV);
    return 0;
}