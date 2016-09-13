#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"


int GzNewRender(GzRender **render, GzDisplay *display)
{
/* 
- malloc a renderer struct
- span interpolator needs pointer to display for pixel writes
*/
	if (display != NULL) {
		*render = (GzRender *)malloc(sizeof(GzRender));
		(*render)->display = display;
		return GZ_SUCCESS;
	}
	return GZ_FAILURE;
}


int GzFreeRender(GzRender *render)
{
/* 
-free all renderer resources
*/
	if (render != NULL) {
		free(render);
	}
	return GZ_SUCCESS;
}


int GzBeginRender(GzRender	*render)
{
/* 
- set up for start of each frame - init frame buffer
*/
	int status = GZ_SUCCESS;
	if (render != NULL && render->display != NULL) {
		status |= GzInitDisplay(render->display);
		render->flatcolor[RED] = 0;
		render->flatcolor[GREEN] = 0;
		render->flatcolor[BLUE] = 0;
	}
	return status;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer *valueList) /* void** valuelist */
{
/*
- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
- later set shaders, interpolaters, texture maps, and lights [another assignment I assume]
*/
	if (render != NULL) {
		for (int i = 0; i < numAttributes; ++i) {
			switch (nameList[i]) {
				case GZ_RGB_COLOR:
					// TODO: Do I have to increment through tokens (ints) and use (sizeof) token type
					// to increment the ponter through the value list
					GzColor * color = (GzColor *)valueList[i];
					// clamp color values
					(*color)[RED] = fmaxf(0, fminf(4095, (*color)[RED]));
					(*color)[GREEN] = fmaxf(0, fminf(4095, (*color)[GREEN]));
					(*color)[BLUE] = fmaxf(0, fminf(4095, (*color)[BLUE]));
					render->flatcolor[RED] = (*color)[RED];
					render->flatcolor[GREEN] = (*color)[GREEN];
					render->flatcolor[BLUE] = (*color)[BLUE];
					break;
				// later set shaders, interpolaters, texture maps, and lights
			}
		}
		return GZ_SUCCESS;
	}
	return GZ_FAILURE;
}


int GzPutTriangle(GzRender *render, int	numParts, GzToken *nameList,
	GzPointer *valueList) 
/* numParts - how many names and values */
{
/* 
- pass in a triangle description with tokens and values corresponding to
      GZ_NULL_TOKEN:		do nothing - no values
      GZ_POSITION:		3 vert positions 
- Invoke the scan converter and return an error code
*/
	// TODO
	
	// Edge Classifiers
	const int UNDEFINED_EDGE = -1;
	const int TOP_EDGE = 0;
	const int BOTTOM_EDGE = 1;
	const int LEFT_EDGE = 2;
	const int RIGHT_EDGE = 3;

	if (render != NULL) {
		const int VERTICES_PER_TRIANGLE = 3;
		GzCoord triangleVertices[VERTICES_PER_TRIANGLE];
		for (int i = 0; i < numParts; ++i) {
			switch (nameList[i])
			{
				case GZ_POSITION:
					// Load all coordinates from provided input
					float yValues[VERTICES_PER_TRIANGLE];
					for (int t = 0; t < VERTICES_PER_TRIANGLE; ++t) {
						triangleVertices[t][X] = ((GzCoord *)(valueList[i]))[t][X];
						triangleVertices[t][Y] = ((GzCoord *)(valueList[i]))[t][Y];
						triangleVertices[t][Z] = ((GzCoord *)(valueList[i]))[t][Z];
						yValues[t] = triangleVertices[t][Y];
					}
					
					// Get indices of vertices corresponding to increasing Y
					int * indices;
					sortTriangleVertices(yValues, &indices);
					int i0 = indices[0];
					int i1 = indices[1];
					int i2 = indices[2];

					// Get vertices in increasing Y
					GzCoord Vertices[3];
					memcpy(Vertices[0], triangleVertices[i0], sizeof(float) * 3);
					memcpy(Vertices[1], triangleVertices[i1], sizeof(float) * 3);
					memcpy(Vertices[2], triangleVertices[i2], sizeof(float) * 3);
					
					// Consistently arrange vertices and generate edges
					// also making sure to note what type of edge it is
					// note: E_i starts with V_i at its tail
					int edgeTypes[] = { UNDEFINED_EDGE, UNDEFINED_EDGE, UNDEFINED_EDGE };

					// Check for Horizontal edges first
					if (Vertices[0][Y] == Vertices[1][Y] && Vertices[1][Y] == Vertices[2][Y]) {
						// All points have same Y value (horizontal lines)
						// TODO: what???
						// Assuming we won't have to deal with this for now?
						return GZ_SUCCESS;
					}
					else if (Vertices[0][Y] == Vertices[1][Y]) {
						// Top Edge
						if (Vertices[0][X] < Vertices[1][X]) {
							// Leave alone
						}
						else if (Vertices[0][X] > Vertices[1][X]) {
							// Swap V0 and V1
							GzCoord temp;
							memcpy(temp, Vertices[0], sizeof(GzCoord));
							memcpy(Vertices[0], Vertices[1], sizeof(GzCoord));
							memcpy(Vertices[1], temp, sizeof(GzCoord));
						}
						else {
							// Top edge is a line in the z-axis (triangle appears as vertical line)
							// TODO: what???
							// Assuming we won't have to deal with this for now?
							return GZ_SUCCESS;
						}
						// E0 is a top edge
						// E1 is a right edge
						// E2 is a left edge
						edgeTypes[0] = TOP_EDGE;
						edgeTypes[1] = RIGHT_EDGE;
						edgeTypes[2] = LEFT_EDGE;
					}
					else if (Vertices[1][Y] == Vertices[2][Y]) {
						// Bottom Edge
						if (Vertices[1][X] < Vertices[2][X]) {
							// Swap V1 and V2
							GzCoord temp;
							memcpy(temp, Vertices[1], sizeof(GzCoord));
							memcpy(Vertices[1], Vertices[2], sizeof(GzCoord));
							memcpy(Vertices[2], temp, sizeof(GzCoord));
						}
						else if (Vertices[1][X] > Vertices[2][X]) {
							// Leave alone
						}
						else {
							// Bottom edge is a line in the z-axis (triangle appears as vertical line)
							// TODO: what???
							// Assuming we won't have to deal with this for now?
							return GZ_SUCCESS;
						}
						// E0 is a right edge
						// E1 is a bottom edge
						// E2 is a left edge
						edgeTypes[0] = RIGHT_EDGE;
						edgeTypes[1] = BOTTOM_EDGE;
						edgeTypes[2] = LEFT_EDGE;
					}
					else {
						// Non-Horizontal Edges
						float dY = Vertices[2][Y] - Vertices[0][Y];
						float dX = Vertices[2][X] - Vertices[0][X];
						float oppositeEdgeX = Vertices[0][X] + (dX / dY) * (Vertices[1][Y] - Vertices[0][Y]);
						if (oppositeEdgeX < Vertices[1][X]) {
							// V1 is R-edge
							// Don't need to do anything (leave order of vertices)
							// E1 is a right edge
							edgeTypes[1] = RIGHT_EDGE;
						}
						else if (oppositeEdgeX > Vertices[1][X])
						{
							// V1 is L-edge
							// Want to switch V1 and V2
							GzCoord temp;
							memcpy(temp, Vertices[1], sizeof(GzCoord));
							memcpy(Vertices[1], Vertices[2], sizeof(GzCoord));
							memcpy(Vertices[2], temp, sizeof(GzCoord));
							// E1 is a left edge
							edgeTypes[1] = LEFT_EDGE;
						}
						else {
							// Weird case where all edges have same slope (triangle appears as line)
							// TODO: what???
							// Assuming we won't have to deal with this for now?
							return GZ_SUCCESS;
						}
						// E0 is a right edge
						// E1 is either a L/R edge (see above)
						// E2 is a left edge
						edgeTypes[0] = RIGHT_EDGE;
						edgeTypes[2] = LEFT_EDGE;
					}

					// Compute projection for vertices, compute the E_i
					float A[3], B[3], C[3];
					for (int j = 0; j < 3; ++j) {
						float dX = Vertices[(j + 1) % 3][X] - Vertices[j][X];
						float dY = Vertices[(j + 1) % 3][Y] - Vertices[j][Y];
						A[j] = dY;
						B[j] = -dX;
						C[j] = dX * Vertices[j][Y] - dY * Vertices[j][X];
					}

					// Compute bbox
					int xmin, xmax, ymin, ymax;
					xmin = floor(fminf(fminf(triangleVertices[0][X], triangleVertices[1][X]), triangleVertices[2][X]));
					xmax = ceil(fmaxf(fmaxf(triangleVertices[0][X], triangleVertices[1][X]), triangleVertices[2][X]));
					ymin = floor(fminf(fminf(triangleVertices[0][Y], triangleVertices[1][Y]), triangleVertices[2][Y]));
					ymax = ceil(fmaxf(fmaxf(triangleVertices[0][Y], triangleVertices[1][Y]), triangleVertices[2][Y]));
					
					// Clip bbox to screen limits
					int xres, yres;
					xres = render->display->xres;
					yres = render->display->yres;
					xmin = max(0, min(xres - 1, xmin));
					xmax = max(0, min(xres - 1, xmax));
					ymin = max(0, min(yres - 1, ymin));
					ymax = max(0, min(yres - 1, ymax));

					// Find plane equation for triangle
					float NA, NB, NC, ND;
					getPlane(Vertices, &NA, &NB, &NC, &ND);

					// For all pixels in bbox
					for (int i = xmin; i < xmax; ++i) {
						for (int j = ymin; j < ymax; ++j) {
							// Evaluate edge functions a_i * x + b_i * y + c_i
							float E0 = A[0] * i + B[0] * j + C[0];
							float E1 = A[1] * i + B[1] * j + C[1];
							float E2 = A[2] * i + B[2] * j + C[2];
							int s0 = sign(E0);
							int s1 = sign(E1);
							int s2 = sign(E2);

							if (s0 == 0) {
								if (!(edgeTypes[0] == LEFT_EDGE || edgeTypes[0] == TOP_EDGE)) {
									// On edge 0 but edge 0 is neither a left or top edge
									continue;
								}
							}
							if (s1 == 0) {
								if (!(edgeTypes[1] == LEFT_EDGE || edgeTypes[1] == TOP_EDGE)) {
									// On edge 1 but edge 1 is neither a left or top edge
									continue;
								}
							}
							if (s2 == 0) {
								if (!(edgeTypes[2] == LEFT_EDGE || edgeTypes[2] == TOP_EDGE)) {
									// On edge 2 but edge 2 is neither a left or top edge
									continue;
								}
							}

							if (!((s0 == s1) && (s1 == s2) && (s0 == s2))) {
								// Outside the triangle
								continue;
							}
							
							// Fill in pixel
							GzIntensity r, g, b, a;
							GzDepth z;
							GzGetDisplay(render->display, i, j, &r, &g, &b, &a, &z);
							// Interpolate z-depth
							float interpZ = interpolateZ(NA, NB, NC, ND, i, j);
							if (interpZ < z || z == 0) {
								// closer - update pixel
								GzPutDisplay(render->display, i, j, ctoi(render->flatcolor[RED]), ctoi(render->flatcolor[GREEN]), ctoi(render->flatcolor[BLUE]), a, interpZ);
							}
						}
					}
						

					break;
			}
		}
		return GZ_SUCCESS;
	}
	return GZ_FAILURE;
}

void sortTriangleVertices(float * values, int ** sortedIndices)
{
	// Default order
	*sortedIndices = (int *)malloc(sizeof(int) * 3);
	(*sortedIndices)[0] = 0;
	(*sortedIndices)[1] = 1;
	(*sortedIndices)[2] = 2;

	float v0 = values[0];
	float v1 = values[1];
	float v2 = values[2];
	
	int tempIndex;
	float tempFloat;

	if (v0 > v1) {
		tempIndex = (*sortedIndices)[0];
		(*sortedIndices)[0] = (*sortedIndices)[1];
		(*sortedIndices)[1] = tempIndex;
		tempFloat = v0;
		v0 = v1;
		v1 = tempFloat;
	}
	if (v1 > v2) {
		tempIndex = (*sortedIndices)[1];
		(*sortedIndices)[1] = (*sortedIndices)[2];
		(*sortedIndices)[2] = tempIndex;
		tempFloat = v1;
		v1 = v2;
		v2 = tempFloat;
	}
	if (v0 > v1) {
		tempIndex = (*sortedIndices)[0];
		(*sortedIndices)[0] = (*sortedIndices)[1];
		(*sortedIndices)[1] = tempIndex;
		tempFloat = v0;
		v0 = v1;
		v1 = tempFloat;
	}
}

int sign(float value) {
	if (value == 0) {
		return 0;
	}
	return fabsf(value) / value;
}

void getPlane(GzCoord * triangleVertices, float * A, float * B, float * C, float * D) {
	float X1 = triangleVertices[1][X] - triangleVertices[0][X];
	float Y1 = triangleVertices[1][Y] - triangleVertices[0][Y];
	float Z1 = triangleVertices[1][Z] - triangleVertices[0][Z];
	float X2 = triangleVertices[2][X] - triangleVertices[0][X];
	float Y2 = triangleVertices[2][Y] - triangleVertices[0][Y];
	float Z2 = triangleVertices[2][Z] - triangleVertices[0][Z];
	*A = Y1*Z2 - Z1*Y2;
	*B = Z1*X2 - X1*Z2;
	*C = X1*Y2 - Y1*X2;
	float x = triangleVertices[0][X];
	float y = triangleVertices[0][Y];
	float z = triangleVertices[0][Z];
	*D = -1 * ((*A)*x + (*B)*y + (*C)*z);
}

float interpolateZ(float A, float B, float C, float D, float x, float y) {
	return -1 * (A*x + B*y + D) / C;
}

/* NOT part of API - just for general assistance */

short	ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}

