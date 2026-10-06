#include <stdio.h>
#include <stdlib.h>

double NormalStandard() {
    double E = 0;
    double factorial = 1.0;

    int nE = 10000;
    for (int i=0; i < nE; i++) {
        // E += 
        E += 1.0 / factorial;
        i++;
        factorial *= i;
    }

    return E;
}

int main() {
    double E = NormalStandard();
    printf("E=%0.8f\r\n", E);
    return 0;
}