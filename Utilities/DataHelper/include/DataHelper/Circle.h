// Circle.h: interface for the Circle class.
// Circle class.
// Purpose : Represent the circle object
// Input : 3 different points
// Process : Calcuate the radius and center
// Output : Circle
//           
// This class originally designed for representation of discretized curvature information 
// of sequential pointlist  
// KJIST CAD/CAM     Ryu, Jae Hun ( ryu@geguri.kjist.ac.kr)
// Last update : 1999. 7. 4


#include "CLHEP/Vector/TwoVector.h"

class Circle  
{
public:
	double GetRadius();
	CLHEP::Hep2Vector* GetCenter();
	Circle(CLHEP::Hep2Vector *p1, CLHEP::Hep2Vector *p2, CLHEP::Hep2Vector *p3);	// p1, p2, p3 are co-planar
	Circle();
	virtual ~Circle();

private:
	double CalcCircle(CLHEP::Hep2Vector *pt1, CLHEP::Hep2Vector *pt2, CLHEP::Hep2Vector *pt3);
	bool IsPerpendicular(CLHEP::Hep2Vector *pt1, CLHEP::Hep2Vector *pt2, CLHEP::Hep2Vector *pt3);
	double m_dRadius;
	CLHEP::Hep2Vector m_Center;
};

