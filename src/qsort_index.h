#include <Rcpp.h>
#include<stdlib.h>
#include<stdio.h>
#include <string.h>
#include<math.h>

int *qsort_index(void *buf, int N, int size, int (*compfunc)(const void *e1, const void *e2));
