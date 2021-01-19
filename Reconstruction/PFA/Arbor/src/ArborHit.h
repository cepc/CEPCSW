#ifndef _ARBORHIT_H_
#define _ARBORHIT_H_

#include "TVector3.h"
#include <iostream>

class ArborHit
{
	float m_hitTime; 
	float m_depth; 
	int m_hitLayer; 
	TVector3 m_hitPos; 
	int m_subD; 
	int m_stave; 

	public:
	
	// ArborHit();
	ArborHit( TVector3 hitPos, int hitLayer, float hitTime, float depth, int stave, int subD );
	
	void setHit( TVector3 hitPos, int hitLayer, float hitTime, float depth, int stave, int subD );

	float GetTime()
	{
		return m_hitTime;
	} 
	int GetLayer()
	{
		return m_hitLayer;
	} 
	TVector3 GetPosition()
	{ 
		return m_hitPos;
	}
	int GetSubD()
	{
		return m_subD; 
	}
	int GetStave()
	{
		return m_stave; 
	}
	float GetDepth()
	{
		return m_depth; 
	}
};

#endif
