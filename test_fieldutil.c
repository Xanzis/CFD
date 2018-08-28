#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "fieldutil_1d.h"

void printBoundary(coordinate coord) {
	boundaryCell b = FLD_getBoundary(coord);
	printf("--\n%f\n%f\n--\n", b.uw, b.ue);
}

int main() {
	FIELDUTIL_setup((int) 5);

	for (int i = -1; i <= 5; i++) {
		int isit = FLD_isBoundary((coordinate) {i});
		printf("%d is a boundary: %d\n", i, isit);
		if (isit) {
			printBoundary((coordinate) {i});
		}
	}
	
	coordinate coord = {5};
	FLD_setU(coord, (float) 7.56);
	printf("\n\nHi %f", FLD_getU(coord));
}