///////////////////////////////////////////////////
//
//	Hamish Carr
//	February, 2019
//
//	------------------------
//	GeometricWindow.cpp
//	------------------------
//	
//	The main window for our little geometric application
//	
///////////////////////////////////////////////////

#include "GeometricWindow.h"
	
// constructor
GeometricWindow::GeometricWindow(GeometricSurfaceDirectedEdge *newSurface, QWidget *parent)
	: QWidget(parent), surface(newSurface)
	{ // GeometricWindow constructor
	setWindowTitle(QString("COMP5821M Assignment 3"));

	// create a layout
	windowLayout = new QGridLayout(this);

	// a label for the slider
	smoothnessLabel = new QLabel("Simplification Parameter", this);

	// and a slider for controlling the computation
	smoothnessSlider = new QSlider(Qt::Horizontal, this);

	//	create the geometric widget
	geometricWidget = new GeometricWidget(surface, this);

	// add the widgets to the layout
	windowLayout->addWidget(geometricWidget, 0, 0, 1, 1);
	windowLayout->setRowStretch(0, 1);
	windowLayout->addWidget(smoothnessLabel, 1, 0, 1, 1);
	windowLayout->addWidget(smoothnessSlider, 2, 0, 1, 1);
	} // GeometricWindow constructor

// destructor
GeometricWindow::~GeometricWindow()
	{ // GeometricWindow destructor
	delete smoothnessSlider;
	delete smoothnessLabel;
	delete geometricWidget;
	} // GeometricWindow destructor

// routine to force interface refresh
void GeometricWindow::ResetInterfaceElements()
	{ // ResetInterfaceElements()
	geometricWidget->update();
	smoothnessLabel->update();
	smoothnessSlider->update();
	} // ResetInterfaceElements()
