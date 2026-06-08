#ifndef _DTABLE_CATALOGUE_BOOKEXAMPLE_GENERAL_H_
#define _DTABLE_CATALOGUE_BOOKEXAMPLE_GENERAL_H_

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

#endif