#ifndef ClusterExtended_H
#define ClusterExtended_H 1

//#include "CaloHitExtended.h"
//#include "TrackExtended.h"
#include "edm4hep/Cluster.h"
#include "HelixClass.h"

using namespace edm4hep;

class TrackExtended;
typedef std::vector<TrackExtended*> TrackExtendedVec;

class CaloHitExtended;
typedef std::vector<CaloHitExtended*> CaloHitExtendedVec;

//class ClusterExtended;
//typedef std::vector<ClusterExtended*> ClusterExtendedVec;

/**
 * Class extending native LCIO class Cluster. <br>
 * Class ClusterExtended is used in TrackwiseClustering <br>
 * and Wolf processors. <br>
 * @author A. Raspereza (DESY)<br>
 * @version $Id: ClusterExtended.h,v 1.4 2006-02-22 14:41:41 owendt Exp $<br>
 */

class ClusterExtended {

 public:

    ClusterExtended();
    ClusterExtended( Cluster * cluster );
    ClusterExtended(CaloHitExtended * calohit);
    ClusterExtended(TrackExtended * track);

    ~ClusterExtended();
    
    CaloHitExtendedVec & getCaloHitExtendedVec();
    TrackExtendedVec & getTrackExtendedVec();
    const float* getStartingPoint();
    const float* getDirection();
    void setStartingPoint(float* sPoint);
    void setDirection(float* direct);
    void addCaloHitExtended(CaloHitExtended * calohit);
    void addTrackExtended(TrackExtended * track);
    void setType( int type );
    int getType();

    void Clear();

    void setCluster(Cluster * cluster);
    Cluster * getCluster();

    void setAxis(float * axis);
    float * getAxis();

    void setEccentricity( float eccentricity);
    float getEccentricity();

    void setHelix(HelixClass helix);
    HelixClass & getHelix();

    void setHelixChi2R(float helixChi2);
    float getHelixChi2R();

    void setHelixChi2Z(float helixChi2);
    float getHelixChi2Z();

    void setPosition(float * position);
    float * getPosition();

    void setLowEdge(float * lowEdge);
    float * getLowEdge();
    void setUpEdge(float * upEdge);
    float * getUpEdge();


 private:

    TrackExtendedVec _trackVector;
    CaloHitExtendedVec _hitVector;
    float _startingPoint[3];
    float _direction[3];

    int _type;
    Cluster * _cluster;

    float _axis[3];
    float _position[3];
    float _eccentricity;

    HelixClass _helix;
    float _helixChi2R;
    float _helixChi2Z;
    
    float _lowEdge[3];
    float _upEdge[3];


};

typedef std::vector<ClusterExtended*> ClusterExtendedVec;

#endif
