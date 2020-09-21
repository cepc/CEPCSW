#ifndef TRACKERHITAR_H
#define TRACKERHITAR_H 1

//#include "lcio.h"
//#include "EVENT/LCIO.h"
#include "edm4hep/TrackerHit.h"
//#include "TrackExtended.h"
#include <vector>


//using namespace edm4hep;

class TrackExtended;
typedef std::vector<TrackExtended*> TrackExtendedVec;

/**
 * Class extending native LCIO class TrackerHit. <br>
 * Class TrackerHitExtended is used in TrackwiseClustering <br>
 * and Wolf processors. <br>
 * @author A. Raspereza (DESY)<br>
 * @version $Id: TrackerHitExtended.h,v 1.7 2007-09-05 09:39:49 rasp Exp $<br>
 */
class TrackerHitExtended {
  
 public:
  
  TrackerHitExtended(const edm4hep::ConstTrackerHit trackerhit);
  ~TrackerHitExtended();
  void setTrackExtended(TrackExtended * trackAR);
    void addTrackExtended(TrackExtended * trackAR);
    void setTrackerHitTo(TrackerHitExtended * hitTo);
    void setTrackerHitFrom(TrackerHitExtended * hitFrom);
    void setGenericDistance(float genericDistance);
    //void setTrackerHit(edm4hep::ConstTrackerHit hit);
    void setYresTo(float yresTo);
    void setYresFrom(float yresFrom);
    void setDirVec(float * dirVec);
    void clearTrackVec();
    void setResolutionRPhi(float rphiReso);
    void setResolutionZ(float zReso);
    void setType(int type);
    void setDet(int idet);
    void setUsedInFit(bool usedInFit);

    edm4hep::ConstTrackerHit getTrackerHit();
    TrackExtended * getTrackExtended();
    TrackExtendedVec & getTrackExtendedVec();
    TrackerHitExtended * getTrackerHitFrom();
    TrackerHitExtended * getTrackerHitTo();
    float getGenericDistance();
    float getYresTo();
    float getYresFrom();
    float * getDirVec();
    float getResolutionRPhi();
    float getResolutionZ();
    int getType();
    int getDet();
    bool getUsedInFit();

 private:

    edm4hep::ConstTrackerHit _trackerHit;
    TrackExtended * _trackAR;
    TrackerHitExtended * _hitTo;
    TrackerHitExtended * _hitFrom;
    TrackExtendedVec _trackVecAR;

    float _rphiReso;
    float _zReso;
    float _yresTo;
    float _yresFrom;
    float _genericDistance;
    float _dirVec[3];

    int _type;
    int _idet;

    bool _usedInFit;
	
};

typedef std::vector<TrackerHitExtended*> TrackerHitExtendedVec; 

#endif
