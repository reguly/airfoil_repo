#ifndef WRITE_VTK
#define WRITE_VTK

void WriteMeshToVTKAscii(const char* filename, float* nodeCoords_data, int nnode, int* cellsToNodes_data, int ncell, float *values_data);
void WriteMeshToVTKBinary(const char* filename, float* nodeCoords_data, int nnode, int* cellsToNodes_data, int ncell, float *values_data);
#endif
