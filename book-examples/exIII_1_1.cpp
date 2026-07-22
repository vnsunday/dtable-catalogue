#include <stdio.h>
#include <stdlib.h>

int main(int argc, char const *argv[]) 
{
    // Table definition 
    double r1[3] = { 0.1, 0.25, 0.1 };
    double r2[3] = { 0.15, 0.05, 0.35 };

    double a[2][3] = {
        { 0.1, 0.25, 0.1 },
        { 0.15, 0.05, 0.35 }
    };

    // 
    const char* szV1 = "X";
    const char* szV2 = "Y";

    double p = 0.0;
    for (int i=0; i<2; i++) {
        for (int j=0; j<3; j++) {
            p += a[i][j];
        }
    }

    //
    printf("P=%0.2f\r\n", p);
    return 0;
}
