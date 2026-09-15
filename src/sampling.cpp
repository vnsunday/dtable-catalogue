#include <stdio.h>
#include <stdlib.h>

void tree_edge(const char* tree, const char* newNode, const char* parent = NULL) {
    // Add Tree Here
}

int main() {
    const char* szObsRDS = "Random Sampling";

    tree_edge(szObsRDS, "X", NULL);
    tree_edge(szObsRDS, "Random Variables", "X");
    tree_edge(szObsRDS, "n", NULL);
    tree_edge(szObsRDS, "Number; Integer", "n");
    
    return 0;
}