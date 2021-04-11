///////////////////////////////////////////////////
//
//	Hamish Carr
//	February, 2019
//
//	------------------------
//	GeometricWindow.h
//	------------------------
//	
//	The main window for our little geometric application
//	
///////////////////////////////////////////////////

#ifndef _GEOMETRIC_WINDOW_H
#define _GEOMETRIC_WINDOW_H

#include <QWidget>
#include <QGridLayout>
#include <QSlider>
#include <QLabel>
#include "GeometricSurfaceDirectedEdge.h"
#include "GeometricWidget.h"

class GeometricWindow : public QWidget										
	{ // class GeometricWindow
	public:	
	// the model - i.e. the surface
	GeometricSurfaceDirectedEdge *surface;

	// the layout widget
	QGridLayout *windowLayout;

	// the main widget
	GeometricWidget *geometricWidget;

	// a label for the slider
	QLabel *smoothnessLabel;

	// the smoothness slider
	QSlider *smoothnessSlider;
	
	// constructor
	GeometricWindow(GeometricSurfaceDirectedEdge *newSurface, QWidget *parent);
	
	// destructor
	~GeometricWindow();
	
	// routine to force interface refresh
	void ResetInterfaceElements();
	}; // class GeometricWindow

#endif
