#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "fieldutil_1d.h"

int main() {
	FIELDUTIL_setup((int) 5);

	//printf("-1 is a boundary: %d\n", FLD_isBoundary((coordinate) {-1}));
}