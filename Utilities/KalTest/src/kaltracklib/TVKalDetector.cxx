//*************************************************************************
//* =====================
//*  TVKalDetector Class
//* =====================
//*
//* (Description)
//*   Base class to hold information of detector system
//*   used in Kalman filter classes.
//* (Requires)
//* 	TObjArray
//* (Provides)
//* 	class TVKalDetector
//* (Update Recored)
//*   2003/09/30  K.Fujii	Original version.
//*   2005/08/14  K.Fujii  Moved GetEnergyLoss() and CalcQms() to
//*                        TVMeasLayer.
//*
//*************************************************************************

#include "TVKalDetector.h"

ClassImp(TVKalDetector)
