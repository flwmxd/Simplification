///////////////////////////////////////////////////
//
//	Hamish Carr
//	February, 2019
//
//	------------------------
//	GeometricController.cpp
//	------------------------
//	
///////////////////////////////////////////////////

#include "GeometricController.h"

// constructor
GeometricController::GeometricController(GeometricSurfaceDirectedEdge *Surface, GeometricWindow *Window)
	: surface(Surface), window(Window)
	{ // constructor
	QObject::connect(	window->smoothnessSlider,	SIGNAL(valueChanged(int)), 
						this,						SLOT(smoothnessSliderChanged(int)));
	} // constructor

// slot for responding to slider changes
void GeometricController::smoothnessSliderChanged(int value)
{ // smoothnessSliderChanged()
    float percentage = ( 100 - (float)value)/100.0;
    int32_t newSize = percentage * surface->position.size();
    std::cout<<"percentage : "<<percentage<<std::endl;

    surface->Start(percentage);

} // smoothnessSliderChanged()
