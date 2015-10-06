
#include "qsort_index.h"

/*
qsort function which returns the index of items sorted.
Uses Shell sort.
Params: buf - the array to be sorted
N - number of elements
size - element width
compfunc - qsort-style comparison function
Returns: malloced array with index numbers of sorted values in
original array.
*/
//// [[Rcpp::export]]
int *qsort_index(void *buf, int N, int size, int (*compfunc)(const void
*e1, const void *e2))
{
int *answer;
int ciura_intervals[] = {701, 301, 132, 57, 23, 10, 4, 1};
int i, ii, iii;
unsigned char *buff = (unsigned char*)buf;
unsigned char *tbuff;
int interval, res;
int t;
int passes;

/* set up return array. Initally all entries are in index order */
answer = (int*)malloc(N * sizeof(int));
if(!answer)
return 0;
for(i=0;i<N;i++)
answer[i] = i;

/* temporary buffer to make data movement easier */
tbuff = (unsigned char*)malloc(size);
if(!tbuff)
{
free(answer);
return 0;
}

/* how many times to pass over the array ? */
passes = (int) (log(N)/log(2.3));
if(passes < 8)
passes = 8;

/* i goes from negative to + 7. 0 to 7 we use the ciura intervals
*/
for(i= -passes+8;i<8;i++)
{
if(i >= 0)
interval = ciura_intervals[i];
else
{
interval = 701;
/* if we've over 701 items to sort, use an intervals at about
2.3 spacing */
for(ii=0;ii<-i;ii++)
interval = (int) (interval * 2.3);
}
if(interval >= N)
continue;

/* now we've got our interval, pass over with the sort */
for(ii=0;ii<N;ii++)
{
t = answer[ii];
memcpy(tbuff, buff + ii * size, size);
for(iii = ii - interval; iii >= 0; iii -= interval)
{
res = (*compfunc)(buff + iii*size, tbuff);
if(res < 0)
break;
memcpy(buff + (iii + interval) * size, buff + iii * size, size);
answer[iii+interval] = answer[iii];
}
memcpy(buff + (iii + interval) * size, tbuff, size);
answer[iii + interval] = t;
}
}

/* temporary buffer no longer needed */
free(tbuff);

return answer;
}
