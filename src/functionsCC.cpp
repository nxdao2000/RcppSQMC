

#include "functionsCC.hpp"

// [[Rcpp::export]]
int getSeed()
{
    ifstream rand("/dev/urandom");
    char tmp[sizeof(int)];
    rand.read(tmp,sizeof(int));
    rand.close();
    int* number = reinterpret_cast<int*>(tmp);
    return (*number);
}
