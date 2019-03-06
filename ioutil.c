#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

// This'll be empty for a little while before I start caring about filesystem-loaded parameters for real computation.

void ioutilerror(char *error_text) {
	printf("IO Utility error. Error message follows:\n%s", error_text);
	printf("\nExiting ...\n");
	exit(1);
}

int IO_intrequest(char *prompt_text) {
	printf(prompt_text);
	int res = scanf("%d", &i);
	return i;
}