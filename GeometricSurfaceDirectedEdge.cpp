///////////////////////////////////////////////////
//
//	Hamish Carr
//	January, 2018
//
//	--------------------------------
//	GeometricSurfaceDirectedEdge.cpp
//	--------------------------------
//
//	Fully implemented Directed Edge code
//
///////////////////////////////////////////////////

#include "GeometricSurfaceDirectedEdge.h"
#include <stdio.h>
#include <math.h>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif
#include <cassert>

// use macros for the "previous" and "next" IDs
#define PREVIOUS_EDGE(x) ((x) % 3) ? ((x) - 1) : ((x) + 2)
#define NEXT_EDGE(x) (((x) % 3) == 2) ? ((x) - 2) : ((x) + 1)

// constructor will initialise to safe values
GeometricSurfaceDirectedEdge::GeometricSurfaceDirectedEdge()
	{ // GeometricSurfaceDirectedEdge::GeometricSurfaceDirectedEdge()
	// force the size to nil (should not be necessary, but . . .)
	position.resize(0);
	firstDirectedEdge.resize(0);
	faceVertices.resize(0);
	otherHalf.resize(0);
	normal.resize(0);

	// set this to something reasonable
	boundingSphereSize = 1.0;

	// set the midpoint to the origin
	midPoint = Cartesian3(0.0, 0.0, 0.0);

	} // GeometricSurfaceDirectedEdge::GeometricSurfaceDirectedEdge()

// read routine returns true on success, failure otherwise
bool GeometricSurfaceDirectedEdge::ReadFileDirEdge(char *fileName
                                                   )
	{ // GeometricSurfaceDirectedEdge::ReadFileDirEdge()
	// these are for accumulating a bounding box for the object

    originalName = fileName;

	Cartesian3 minCoords(1000000.0, 1000000.0, 1000000.0);
	Cartesian3 maxCoords(-1000000.0, -1000000.0, -1000000.0);

	// open the input file
	FILE *inFile = fopen(fileName, "r");
	if (inFile == NULL)
		return false;

	// set the number of vertices and faceVertices
	indexType nVertices = 0, nFaces = 0;

	// set the midpoint to the origin
	midPoint = Cartesian3(0.0, 0.0, 0.0);

	// assume that file format is ABSOLUTELY predictable and the lines aren't too long
	char lineBuffer[256];

	// the first four lines will be skipped completely
	fgets(lineBuffer, 256, inFile);
	fgets(lineBuffer, 256, inFile);
	fgets(lineBuffer, 256, inFile);
	fgets(lineBuffer, 256, inFile);

	// read the only header line we care about
	fscanf(inFile, "# %s vertices=%ld faces=%ld \n", lineBuffer, &nVertices, &nFaces);

	// the next line is skipped
	fgets(lineBuffer, 256, inFile);

	// now allocate space for them all
	position.resize(nVertices);
	normal.resize(nVertices);
	firstDirectedEdge.resize(nVertices);
	faceVertices.resize(nFaces*3);
	otherHalf.resize(nFaces*3);	

	// now loop to read the vertices in, and hope nothing goes wrong
	for (indexType vertex = 0; vertex < nVertices; vertex++)
		{ // for each vertex
		indexType vertexID;
		// scan in the vertex ID
		indexType nScanned = fscanf(inFile, "Vertex %ld", &vertexID);
		// if it didn't scan, or was the wrong ID, exit
		if ((nScanned != 1) || (vertexID != vertex))
			{ // scan failed
			printf("Invalid vertex %ld\n", vertex);
			return false;
			} // scan failed
		// now try to read in the vertex (failure => non-zero)
		if (3 != fscanf(inFile, "%f %f %f ", &(position[vertex].x), &(position[vertex].y), &(position[vertex].z)))
			{ // scan failed
			printf("Invalid vertex %ld\n", vertex);
			return false;
			} // scan failed

		// keep running track of midpoint, &c.
		midPoint = midPoint + position[vertex];
		if (position[vertex].x < minCoords.x) minCoords.x = position[vertex].x;
		if (position[vertex].y < minCoords.y) minCoords.y = position[vertex].y;
		if (position[vertex].z < minCoords.z) minCoords.z = position[vertex].z;

		if (position[vertex].x > maxCoords.x) maxCoords.x = position[vertex].x;
		if (position[vertex].y > maxCoords.y) maxCoords.y = position[vertex].y;
		if (position[vertex].z > maxCoords.z) maxCoords.z = position[vertex].z;
		} // for each vertex

	// now sort out the size of a bounding sphere for viewing
	// and also set the midpoint's location
	midPoint = midPoint / position.size();

	// now go back through the vertices, subtracting the mid point
	for (int vertex = 0; vertex < nVertices; vertex++)
		{ // per vertex
        //position[vertex] = position[vertex] - midPoint;
		} // per vertex

	// the bounding sphere radius is just half the distance between these
	boundingSphereSize = sqrt((maxCoords - minCoords).length()) * 0.5;

	// now we read in the vertex normal vectors
	for (indexType vertex = 0; vertex < nVertices; vertex++)
		{ // for each vertex
		indexType vertexID;
		// scan in the vertex ID
		indexType nScanned = fscanf(inFile, "Normal %ld", &vertexID);
		// if it didn't scan, or was the wrong ID, exit
		if ((nScanned != 1) || (vertexID != vertex))
			{ // scan failed
			printf("Invalid vertex ID %ld\n", vertex);
			return false;
			} // scan failed
		// now try to read in the vertex (failure => non-zero)
		if (3 != fscanf(inFile, "%f %f %f ", &(normal[vertex].x), &(normal[vertex].y), &(normal[vertex].z)))
			{ // scan failed
			printf("Invalid vertex normal %ld\n", vertex);
			return false;
			} // scan failed
		} // for each vertex

	// next, we read in the first directed edge array
	for (indexType vertex = 0; vertex < nVertices; vertex++)
		{ // for each vertex
		indexType vertexID;
		// scan in the vertex ID
		indexType nScanned = fscanf(inFile, "FirstDirectedEdge %ld", &vertexID);
		// if it didn't scan, or was the wrong ID, exit
		if ((nScanned != 1) || (vertexID != vertex))
			{ // scan failed
			printf("Invalid vertex ID for first directed edge %ld %ld\n", vertexID, vertex);
			return false;
			} // scan failed
		// now try to read in the first directed edge (failure => non-zero)
		if (1 != fscanf(inFile, "%ld ", &(firstDirectedEdge[vertex])))
			{ // scan failed
			printf("Invalid first directed edge %ld for vertex %ld\n", firstDirectedEdge[vertex], vertex);
			return false;
			} // scan failed
		} // for each vertex

	// similarly, loop to read the faceVertices in
	for (indexType face = 0; face < nFaces; face++)
		{ // for each face
		indexType faceID;

		// scan in the vertex ID
		indexType nScanned = fscanf(inFile, "Face %ld %ld %ld %ld ",
			&faceID, &(faceVertices[3*face]), &(faceVertices[3*face+1]), &(faceVertices[3*face+2]));

		// if it didn't scan, or was the wrong ID, exit
		if ((nScanned != 4) || (faceID != face))
			{ // scan failed
			printf("Invalid face %ld\n", face);
			return false;
			} // scan failed
		} // for each face

	// repeat for other half
	for (indexType dirEdge = 0; dirEdge < nFaces * 3; dirEdge++)
		{ // for each directed edge
		indexType dirEdgeID;

		// scan in the vertex ID
		indexType nScanned = fscanf(inFile, "OtherHalf %ld %ld ", &dirEdgeID, &(otherHalf[dirEdge]));

		// if it didn't scan, or was the wrong ID, exit
		if ((nScanned != 2) || (dirEdgeID != dirEdge))
			{ // scan failed
			printf("Invalid directed edge ID %ld %ld\n", dirEdgeID, dirEdge);
			return false;
			} // scan failed
		} // for each directed edge

    originalPosition = position;
    originalFaceVertices = faceVertices;

	return true;
	} // GeometricSurfaceDirectedEdge::ReadFileDirEdge()

// write routine
bool GeometricSurfaceDirectedEdge::WriteFileDirEdge(const char *fileName,
                                                    const std::vector<Cartesian3> & position,
                                                    const std::vector<indexType> & firstDirectedEdge,
                                                    const std::vector<indexType> & faceVertices,
                                                    const std::vector<indexType> & otherHalf,
                                                    const std::vector<Cartesian3> & normal

                                                    )
	{ // GeometricSurfaceDirectedEdge::WriteFileDirEdge()
	// open the input file
	FILE *outFile = fopen(fileName, "w");
	if (outFile == NULL)
		return false;

	// the first four lines will be skipped completely
	fprintf(outFile, "#\n");
	fprintf(outFile, "# Created for Leeds COMP 5821M Spring 2018\n");
	fprintf(outFile, "#\n");
	fprintf(outFile, "#\n");
	fprintf(outFile, "# Surface vertices=%ld faces=%ld\n", position.size(), faceVertices.size()/3);
	fprintf(outFile, "#\n");

	// loop to write the vertices out
	for (indexType vertex = 0; vertex < (indexType) position.size(); vertex++)
 		fprintf(outFile, "Vertex %1ld %9.6f %9.6f %9.6f\n", vertex, position[vertex].x, position[vertex].y, position[vertex].z);

	// loop to write out the normals
	for (indexType vertex = 0; vertex < (indexType) position.size(); vertex++)
 		fprintf(outFile, "Normal %1ld %9.6f %9.6f %9.6f\n", vertex, normal[vertex].x, normal[vertex].y, normal[vertex].z);

	// loop to write out the first directed edge index
	for (indexType vertex = 0; vertex < (indexType) position.size(); vertex++)
		fprintf(outFile, "FirstDirectedEdge %1ld  %2ld\n", vertex, firstDirectedEdge[vertex]);

	// loop to write the face vertices out
	for (indexType face = 0; face < (indexType) faceVertices.size()/3; face++)
		fprintf(outFile, "Face %2ld  %2ld %2ld %2ld\n", face, faceVertices[3*face], faceVertices[3*face+1], faceVertices[3*face+2]);

	// loop to write out the edge pairs
	for (indexType dirEdge = 0; dirEdge < (indexType) faceVertices.size(); dirEdge++)
        fprintf(outFile, "OtherHalf %1ld  %2ld\n", dirEdge, otherHalf[dirEdge]);

	// close the file & return
	fclose(outFile);
	return true;
	} // GeometricSurfaceDirectedEdge::WriteFileDirEdge()

// routine to render
void GeometricSurfaceDirectedEdge::Render()
	{ // GeometricSurfaceDirectedEdge::Render()
	glBegin(GL_TRIANGLES);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

	// since normals pair with vertices, loop is simpler
	for (unsigned int face = 0; face < faceVertices.size(); face ++)
		{ // per face
        Cartesian3 v0 = (position[faceVertices[face]]) - midPoint;
        Cartesian3 *n0 = &(normal[faceVertices[face]]);

        glNormal3fv((GLfloat *) n0);
        glVertex3fv((GLfloat *) &v0);
		} // per face
	glEnd();
	} // GeometricSurfaceDirectedEdge::Render()

// debug routine to dump arrays
void GeometricSurfaceDirectedEdge::DumpArrays()
	{ // GeometricSurfaceDirectedEdge::DumpArrays()
	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("Dumping Arrays\n");
	// print a header line
	printf("Vertex ID:    |");
	for (indexType vertex = 0; vertex < (indexType) position.size(); vertex++)
		printf("%9ld ", vertex);
	printf("\n");
	printf("--------------+");
	for (indexType vertex = 0; vertex < (indexType) position.size(); vertex++)
		printf("----------");
	printf("\n");
	// print the data
	printf("Position x    |");
	for (indexType vertex = 0; vertex < (indexType) position.size(); vertex++)
		printf("%9.5f ", position[vertex].x);
	printf("\n");
	printf("Position y    |");
	for (indexType vertex = 0; vertex < (indexType) position.size(); vertex++)
		printf("%9.5f ", position[vertex].y);
	printf("\n");
	printf("Position z    |");
	for (indexType vertex = 0; vertex < (indexType) position.size(); vertex++)
		printf("%9.5f ", position[vertex].z);
	printf("\n");
	printf("Normal x      |");
	for (indexType vertex = 0; vertex < (indexType) position.size(); vertex++)
		printf("%9.5f ", normal[vertex].x);
	printf("\n");
	printf("Normal y      |");
	for (indexType vertex = 0; vertex < (indexType) position.size(); vertex++)
		printf("%9.5f ", normal[vertex].y);
	printf("\n");
	printf("Normal z      |");
	for (indexType vertex = 0; vertex < (indexType) position.size(); vertex++)
		printf("%9.5f ", normal[vertex].z);
	printf("\n");
	printf("First DirEdge |");
	for (indexType vertex = 0; vertex < (indexType) position.size(); vertex++)
		printf("%9ld ", firstDirectedEdge[vertex]);
	printf("\n");
	// print a header line
	printf("Vertex ID:    |");
	for (indexType dirEdge = 0; dirEdge < (indexType) faceVertices.size(); dirEdge++)
		printf("%9ld ", dirEdge);
	printf("\n");
	printf("--------------+");
	for (indexType dirEdge = 0; dirEdge < (indexType) faceVertices.size(); dirEdge++)
		printf("----------");
	printf("\n");
	printf("Face Vertices |");
	for (indexType dirEdge = 0; dirEdge < (indexType) faceVertices.size(); dirEdge++)
		printf("%9ld ", faceVertices[dirEdge]);
	printf("\n");
	printf("Other Half    |");
	for (indexType dirEdge = 0; dirEdge < (indexType) faceVertices.size(); dirEdge++)
		printf("%9ld ", otherHalf[dirEdge]);
	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("\n");
	} // GeometricSurfaceDirectedEdge::DumpArrays()


void GeometricSurfaceDirectedEdge::Start(float ratio)
{
    //pass params to simplier
    simplier.position = originalPosition;
    simplier.faceVertices = originalFaceVertices;
    //begin to simplify
    simplier.simplification(ratio);

    faceVertices.clear();
    position = simplier.position;
    //after simplifying, create new mesh.
    simplier.createNewMesh(faceVertices);
}
