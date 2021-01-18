#include "ArborHit.h"
#include <iostream>
#include <vector>

ArborHit::ArborHit( TVector3 hitPos, int hitLayer, float hitTime, float depth, int stave, int subD )
{
	setHit( hitPos, hitLayer, hitTime, depth, stave, subD );
}

void ArborHit::setHit(TVector3 hitPos, int hitLayer, float hitTime, float depth, int stave, int subD )
{
	m_hitTime = hitTime;
	m_hitLayer = hitLayer; 
	m_hitPos = hitPos; 
	m_subD = subD; 
	m_stave = stave; 
	m_depth = depth;
}
