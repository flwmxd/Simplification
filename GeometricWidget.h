///////////////////////////////////////////////////
//
//	Hamish Carr
//	January, 2018
//
//	------------------------
//	GeometricWidget.h
//	------------------------
//	
//	The main widget that shows the geometry
//	
///////////////////////////////////////////////////

#ifndef _GEOMETRIC_WIDGET_H
#define _GEOMETRIC_WIDGET_H

#include <QGLWidget>
#include <QMouseEvent>
#include "GeometricSurfaceDirectedEdge.h"
#include "Ball.h"

class GeometricWidget : public QGLWidget										
	{ // class GeometricWidget
	Q_OBJECT
	public:	
	// the model - i.e. the surface
	GeometricSurfaceDirectedEdge *surface;
	
	// arcball for storing light rotation
	BallData lightBall;

	// arcball for storing object rotation
	BallData objectBall;
	
	// translation in window x,y
	GLfloat translate_x, translate_y;
	GLfloat last_x, last_y;
	
	// which button was last pressed
	int whichButton;

	// constructor
	GeometricWidget(GeometricSurfaceDirectedEdge *newSurface, QWidget *parent);
	
	// destructor
	~GeometricWidget();
			
	protected:
	// called when OpenGL context is set up
	void initializeGL();
	// called every time the widget is resized
	void resizeGL(int w, int h);
	// called every time the widget needs painting
	void paintGL();

	// mouse-handling
	virtual void mousePressEvent(QMouseEvent *event);
	virtual void mouseMoveEvent(QMouseEvent *event);
	virtual void mouseReleaseEvent(QMouseEvent *event);

	}; // class GeometricWidget

#endif
