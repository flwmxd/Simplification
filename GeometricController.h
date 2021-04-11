///////////////////////////////////////////////////
//
//	Hamish Carr
//	February, 2019
//
//	------------------------
//	GeometricController.h
//	------------------------
//	
///////////////////////////////////////////////////

#include "GeometricSurfaceDirectedEdge.h"
#include "GeometricWindow.h"

#ifndef __GEOMETRIC_CONTROLLER_H
#define __GEOMETRIC_CONTROLLER_H 1

class GeometricController: public QObject
	{ // class GeometricController
	Q_OBJECT
	public:
	// the model to be edited
	GeometricSurfaceDirectedEdge *surface;
	// the window to display in
	GeometricWindow *window;

	// constructor
	GeometricController(GeometricSurfaceDirectedEdge *Surface, GeometricWindow *Window);

	public slots:
	// slot for responding to slider changes
	void smoothnessSliderChanged(int value);

	}; // class GeometricController

#endif