///////////////////////////////////////////////////
//
//	Hamish Carr
//	January, 2018
//
//	------------------------
//	GeometricSurfaceDirectedEdge.h
//	------------------------
//	
//	Base code for geometric assignments.
//
//	This is the minimalistic Face-based D/S for storing
//	surfaces, to be used as the basis for fuller versions
//	
//	It will include object load / save code & render code
//	
///////////////////////////////////////////////////

#ifndef _GEOMETRIC_SURFACE_DIRECTED_EDGE_H
#define _GEOMETRIC_SURFACE_DIRECTED_EDGE_H

#include <vector>
#include <cstdint>
#include <queue>
#include "EdgeHeap.h"
#include <set>
#include <map>
#include <unordered_map>

// an index type - signed so we can use negative numbers as flags
typedef signed long indexType;

// define a macro for "not used" flag
#define NO_SUCH_ELEMENT -1

#include "Cartesian3.h"
#include "Matrix4.h"

#include "Simplifier.h"
#include "Tarjan.h"


class GeometricSurfaceDirectedEdge
	{ // class GeometricSurfaceDirectedEdge
	public:
	// the vertex positions
	std::vector<Cartesian3> position;
    std::vector<Cartesian3> originalPosition;
    std::vector<indexType> originalFaceVertices;



	// the "first" directed edge for each vertex
	std::vector<indexType> firstDirectedEdge;
	// the face vertices - doubles as the "to" array for edges
	std::vector<indexType> faceVertices;
	// the other half of the directed edges
	std::vector<indexType> otherHalf;

	// array to hold the normals
	std::vector<Cartesian3> normal;

	// bounding sphere size
	float boundingSphereSize;

	// midpoint of object
	Cartesian3 midPoint;

	// constructor will initialise to safe values
	GeometricSurfaceDirectedEdge();
		
	// read routine returns true on success, failure otherwise
	// does *NOT* check consistency
	bool ReadFileDirEdge(char *fileName);

	// write routine
    bool WriteFileDirEdge(const char *fileName,
                          const std::vector<Cartesian3> & position,
                          const std::vector<indexType> & firstDirectedEdge,
                          const std::vector<indexType> & faceVertices,
                          const std::vector<indexType> & otherHalf,
                          const std::vector<Cartesian3> & normal
                          );
	
	// routine to render
	void Render();
	
	// debug routine to dump arrays
	void DumpArrays();

    void Start(float ratio);

    Simplifier simplier;

    Tarjan tarjan;

    const char * originalName = nullptr;

	}; // class GeometricSurfaceDirectedEdge

#endif
