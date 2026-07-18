#include <stdlib.h>
#include <stdio.h>
#include "general.hpp"

int main(int argc, char const *argv[])
{
    double d1 = combination(2,3) /  (double)combination(2,9);
    printf("%05.f\r\n", d1)    ;
    return 0;
}
