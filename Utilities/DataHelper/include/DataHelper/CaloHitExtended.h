#ifndef CaloHitExtended_H
#define CaloHitExtended_H 1

//#include "lcio.h"
//#include "EVENT/LCIO.h"
#include "edm4hep/CalorimeterHit.h"
#include "ClusterExtended.h"

using namespace edm4hep ;

class ClusterExtended;

/**
 * Class extending native LCIO class CalorimeterHit. <br>
 * Class CaloHitExtended is used in TrackwiseClustering <br>
 * and Wolf processors. <br>
 * @author A. Raspereza (DESY)<br>
 * @version $Id: CaloHitExtended.h,v 1.4 2006-02-22 14:53:27 owendt Exp $<br>
 */
class CaloHitExtended {

    public:

  CaloHitExtended(const edm4hep::CalorimeterHit& calhit, int type);
    
    ~CaloHitExtended();

    CalorimeterHit * getCalorimeterHit();
    CaloHitExtended * getCaloHitTo();
    CaloHitExtended * getCaloHitFrom();
    ClusterExtended * getClusterExtended();
    int    getIndex();
    int    getType();
    const float* getDirVec();
    float  getYresTo();
    float  getYresFrom();
    float  getGenericDistance();

    //void setCalorimeterHit(CalorimeterHit* calhit);
    void setCaloHitTo(CaloHitExtended* calhitto);
    void setCaloHitFrom(CaloHitExtended* calohitfrom);
    void setClusterExtended(ClusterExtended* cluster);
    void setIndex(int index);
    void setType(int type);
    void setDirVec(float* dirVec);
    void setYresTo(float yresto);
    void setYresFrom(float yresfrom);
    void setGenericDistance(float genericDistance);
    void setDistanceToCalo(float distanceToCalo);
    float getDistanceToCalo();
    void setDistanceToNearestHit(float distanceToNearest);
    float getDistanceToNearestHit();
      

    private:
    
    edm4hep::CalorimeterHit _calohit;
    CaloHitExtended * _calohitTo;
    CaloHitExtended * _calohitFrom;
    ClusterExtended * _clusterAR;
    int _index;
    int _type;
    float _dirVec[3];
    float _yresTo;
    float _yresFrom;
    float _genericDistance;
    float _caloDistance;
    float _distanceToNearestHit;

};

#endif
