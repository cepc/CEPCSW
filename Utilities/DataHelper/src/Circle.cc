// Circle.cpp: implementation of the Circle class.
//
//////////////////////////////////////////////////////////////////////

#include "DataHelper/Circle.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Circle::Circle()
{
	this->m_dRadius=-1;		// error checking 
}

Circle::~Circle()
{

}

Circle::Circle(CLHEP::Hep2Vector *V1, CLHEP::Hep2Vector *V2, CLHEP::Hep2Vector *V3)
{
	this->m_dRadius=-1;		// error checking 

	CLHEP::Hep2Vector *pt1=new CLHEP::Hep2Vector;	CLHEP::Hep2Vector *pt2=new CLHEP::Hep2Vector;	CLHEP::Hep2Vector *pt3=new CLHEP::Hep2Vector;

	*pt1=*V1;			*pt2=*V2;			*pt3=*V3;
	
	if (!this->IsPerpendicular(pt1, pt2, pt3) )			this->CalcCircle(pt1, pt2, pt3);	
	else if (!this->IsPerpendicular(pt1, pt3, pt2) )		this->CalcCircle(pt1, pt3, pt2);	
	else if (!this->IsPerpendicular(pt2, pt1, pt3) )		this->CalcCircle(pt2, pt1, pt3);	
	else if (!this->IsPerpendicular(pt2, pt3, pt1) )		this->CalcCircle(pt2, pt3, pt1);	
	else if (!this->IsPerpendicular(pt3, pt2, pt1) )		this->CalcCircle(pt3, pt2, pt1);	
	else if (!this->IsPerpendicular(pt3, pt1, pt2) )		this->CalcCircle(pt3, pt1, pt2);	
	else { 
	  std::cout << "The three pts are perpendicular to axis" << std::endl;
	  //		pt1->trace();			pt2->trace();			pt3->trace();
	  delete pt1;				delete pt2;				delete pt3;
	  this->m_dRadius=-1;
	  return ;
	}
	delete pt1;				delete pt2;				delete pt3;
	
}

bool Circle::IsPerpendicular(CLHEP::Hep2Vector *pt1, CLHEP::Hep2Vector *pt2, CLHEP::Hep2Vector *pt3)
// Check the given point are perpendicular to x or y axis 
{
	double yDelta_a= pt2->y() - pt1->y();
	double xDelta_a= pt2->x() - pt1->x();
	double yDelta_b= pt3->y() - pt2->y();
	double xDelta_b= pt3->x() - pt2->x();
	

//	TRACE(" yDelta_a: %f xDelta_a: %f \n",yDelta_a,xDelta_a);
//	TRACE(" yDelta_b: %f xDelta_b: %f \n",yDelta_b,xDelta_b);

	// checking whether the line of the two pts are vertical
	if (fabs(xDelta_a) <= 0.000000001 && fabs(yDelta_b) <= 0.000000001){
	  std::cout << "The points are pependicular and parallel to x-y axis" << std::endl;
	  return false;
	}

	if (fabs(yDelta_a) <= 0.0000001){
//		TRACE(" A line of two point are perpendicular to x-axis 1\n");
		return true;
	}
	else if (fabs(yDelta_b) <= 0.0000001){
//		TRACE(" A line of two point are perpendicular to x-axis 2\n");
		return true;
	}
	else if (fabs(xDelta_a)<= 0.000000001){
//		TRACE(" A line of two point are perpendicular to y-axis 1\n");
		return true;
	}
	else if (fabs(xDelta_b)<= 0.000000001){
//		TRACE(" A line of two point are perpendicular to y-axis 2\n");
		return true;
	}
	else return false ;
}

double Circle::CalcCircle(CLHEP::Hep2Vector *pt1, CLHEP::Hep2Vector *pt2, CLHEP::Hep2Vector *pt3)
{
  double yDelta_a= pt2->y() - pt1->y();
  double xDelta_a= pt2->x() - pt1->x();
  double yDelta_b= pt3->y() - pt2->y();
  double xDelta_b= pt3->x() - pt2->x();
  
  if (fabs(xDelta_a) <= 0.000000001 && fabs(yDelta_b) <= 0.000000001){
    // TRACE("Calc cirlce \n");
    this->m_Center.setX(0.5*(pt2->x() + pt3->x()));
    this->m_Center.setY(0.5*(pt1->y() + pt2->y()));
    //    this->m_dRadius = m_Center.howNear(*pt1);		// calc. radius
    this->m_dRadius = sqrt((pt1->x()-m_Center.x())*(pt1->x()-m_Center.x()) + (pt1->y()-m_Center.y())*(pt1->y()-m_Center.y()));
    //		TRACE(" Center: %f %f %f\n", m_Center.x(), m_Center.y(), m_Center.z());
    //		TRACE(" radius: %f %f %f\n", length(&m_Center,pt1), length(&m_Center,pt2),length(&m_Center,pt3));
    
    return this->m_dRadius;
  }
  
  // IsPerpendicular() assure that xDelta(s) are not zero
  double aSlope=yDelta_a/xDelta_a; // 
  double bSlope=yDelta_b/xDelta_b;
  if (fabs(aSlope-bSlope) <= 0.000000001){	// checking whether the given points are colinear. 	
    std::cout << "The three pts are colinear" << std::endl;
		return -1;
	}

  // calc center
  this->m_Center.setX((aSlope*bSlope*(pt1->y() - pt3->y()) + bSlope*(pt1->x() + pt2->x())
		       - aSlope*(pt2->x()+pt3->x()) )/(2* (bSlope-aSlope) ));
  this->m_Center.setY(-1*(m_Center.x() - (pt1->x()+pt2->x())/2)/aSlope +  (pt1->y()+pt2->y())/2);
  
  //this->m_dRadius= m_Center.howNear(*pt1);			// calc. radius
  this->m_dRadius = sqrt((pt1->x()-m_Center.x())*(pt1->x()-m_Center.x()) + (pt1->y()-m_Center.y())*(pt1->y()-m_Center.y()));
  //	TRACE(" Center: %f %f %f\n", m_Center.x(), m_Center.y(), m_Center.z());
  //	TRACE(" radius: %f %f %f\n", length(&m_Center,pt1), length(&m_Center,pt2),length(&m_Center,pt3));
  return this->m_dRadius;
}

CLHEP::Hep2Vector* Circle::GetCenter()
{
  return &this->m_Center;
  
}

double Circle::GetRadius()
{
  return this->m_dRadius;
}
