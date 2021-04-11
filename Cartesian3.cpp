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

#include "Cartesian3.h"
#include "math.h"
#include <stdexcept>

// constructors
Cartesian3::Cartesian3() 
	: x(0.0), y(0.0), z(0.0) 
	{}

Cartesian3::Cartesian3(float X, float Y, float Z)
	: x(X), y(Y), z(Z) 
	{}

Cartesian3::Cartesian3(const Cartesian3 &other)
	: x(other.x), y(other.y), z(other.z) 
	{}
	
// equality operator
bool Cartesian3::operator ==(const Cartesian3 &other) const
	{ // Cartesian3::operator ==()
	return ((x == other.x) && (y == other.y) && (z == other.z));
	} // Cartesian3::operator ==()

// addition operator
Cartesian3 Cartesian3::operator +(const Cartesian3 &other) const
	{ // Cartesian3::operator +()
	Cartesian3 returnVal(x + other.x, y + other.y, z + other.z);
	return returnVal;
	} // Cartesian3::operator +()

// subtraction operator
Cartesian3 Cartesian3::operator -(const Cartesian3 &other) const
	{ // Cartesian3::operator -()
	Cartesian3 returnVal(x - other.x, y - other.y, z - other.z);
	return returnVal;
	} // Cartesian3::operator -()

// multiplication operator
Cartesian3 Cartesian3::operator *(float factor) const
	{ // Cartesian3::operator *()
	Cartesian3 returnVal(x * factor, y * factor, z * factor);
	return returnVal;
	} // Cartesian3::operator *()

// division operator
Cartesian3 Cartesian3::operator /(float factor) const
	{ // Cartesian3::operator /()
	Cartesian3 returnVal(x / factor, y / factor, z / factor);
	return returnVal;
	} // Cartesian3::operator /()

// crossproduct routine
Cartesian3 Cartesian3::cross(const Cartesian3 &other)
	{ // Cartesian3::operator ==()
	Cartesian3 returnVal(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
	return returnVal;
	} // Cartesian3::operator ==()


float Cartesian3::operator[](int index) const
{
    switch (index) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
    }
    throw std::runtime_error("index should be less 3 and greater than -1");
}

float &Cartesian3::operator [] (const int index)
{
    switch (index) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
    }
}
// dot product routine
float Cartesian3::dot(const Cartesian3 &other) const
{ // Cartesian3::dot()
    float returnVal = x * other.x + y * other.y + z * other.z;
    return returnVal;
} // Cartesian3::dot()

// routine to find the length
float Cartesian3::length()
	{ // Cartesian3::length()
	return (x*x + y*y + z*z);	
	} // Cartesian3::length()

// normalisation routine
Cartesian3 Cartesian3::normalise() const
	{ // Cartesian3::normalise()
	float length = sqrt(x*x+y*y+z*z);
	Cartesian3 returnVal(x/length, y/length, z/length);
	return returnVal;
	} // Cartesian3::normalise()
