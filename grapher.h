#ifndef _GRAPHER_H_
#define _GRAPHER_H_

void GRAPH_setup(int canvas_dim, char itemnum, float minvalue, float maxvalue);
void GRAPH_reset();
void GRAPH_addvector(float *invector, int loc);
void GRAPH_addmatrix(float **inmatrix);

void GRAPH_save();

#endif