
#include"functionsCC.hpp"
#include "DigitalNetsBase2.h"

void Scrambling(Scrambled*);

void getPoints(Scrambled*, int, int, double*);



void Digital_getPoints(DigitalNetGenerator*, int, int, double*);

DigitalNetGenerator *DigitalNetGenerator_Create(int, int);
void DigitalNetGenerator_Destroy(DigitalNetGenerator*);
void DigitalNetGenerator_GetPoints(DigitalNetGenerator*, int, int, double*);

Scrambled*  Scrambled_Create(int, int, int);
void Scrambled_Randomize(Scrambled*);
void Scrambled_GetPoints(Scrambled*,int, int, double*);
void Scrambled_Destroy(Scrambled*);





