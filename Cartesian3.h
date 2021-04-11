///////////////////////////////////////////////////
//
//	Hamish Carr
//	January, 2018
//
//	------------------------
//	Cartesian3.h
//	------------------------
//	
//	A minimal class for a point in Cartesian space
//	
///////////////////////////////////////////////////

#ifndef CARTESIAN3_H
#define CARTESIAN3_H

// the class - we will rely on POD for sending to GPU
class Cartesian3
	{ // Cartesian3
	public:
	// the coordinates
	float x, y, z;

	// constructors
	Cartesian3();
	Cartesian3(float X, float Y, float Z);
	Cartesian3(const Cartesian3 &other);
	
	// equality operator
    bool operator ==(const Cartesian3 &other) const ;

	// addition operator
    Cartesian3 operator +(const Cartesian3 &other) const ;

	// subtraction operator
    Cartesian3 operator -(const Cartesian3 &other) const ;
	
	// multiplication operator
    Cartesian3 operator *(float factor) const ;

	// division operator
    Cartesian3 operator /(float factor) const ;

	// crossproduct routine
	Cartesian3 cross(const Cartesian3 &other);
	

    // dot routine
    float dot(const Cartesian3 &other) const;


    float operator[](int index) const;
    float &operator [] (const int index);
	// routine to find the length
	float length();
	
	// normalisation routine
    Cartesian3 normalise() const;

	}; // Cartesian3
	
	
#endif
