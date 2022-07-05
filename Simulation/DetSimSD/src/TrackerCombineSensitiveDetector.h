// *********************************************************
//
// $Id: TrackerCombineSensitiveDetector.hh,v 1.0 2022/03/27

#ifndef TrackerCombineSensitiveDetector_h
#define TrackerCombineSensitiveDetector_h

#include "DetSimSD/DDG4SensitiveDetector.h"
#include "DDG4/Defs.h"

namespace dd4hep {
  namespace sim {
    struct TrackerCombine {
      Geant4TrackerHit  pre;
      Geant4TrackerHit  post;
      G4Track*          track;
      double            e_cut;
      int               current;
      long long int     cellID;
      TrackerCombine() : pre(), post(), track(0), e_cut(0.0), current(-1), cellID(0)  {}
      void start(long long int cell, G4Step* step, G4StepPoint* point)   {
	pre.storePoint(step,point);
	current = pre.truth.trackID;
	track   = step->GetTrack();
	cellID  = cell;
	post    = pre;
	//std::cout << "start: " << cellID << " " << pre.position << " current = " << current << std::endl;
      }
      void update(G4Step* step) {
	post.storePoint(step,step->GetPostStepPoint());
	pre.truth.deposit += post.truth.deposit;
	//std::cout << "update: " << cellID << " " << post.position << " current = " << current << std::endl;
      }
      void clear()   {
	pre.truth.clear();
	current = -1;
	track = 0;
      }
      Geant4TrackerHit* extractHit(DDG4SensitiveDetector::HitCollection* c) {
	//std::cout << "create Geant4TrackerHit: " << cellID << " current = " << current << " track = " << track << " de = " << pre.truth.deposit << std::endl;  
	if ( current == -1 || !track ) {
	  return 0;
	}
	else if ( pre.truth.deposit <= e_cut && !Geant4Hit::isGeantino(track) ) {
	  clear();
	  return 0;
	}
	double rho1 = pre.position.Rho();
	double rho2 = post.position.Rho();
	double rho = 0.5*(rho1+rho2);
	Position pos = 0.5 * (pre.position + post.position);
	double z = pos.z();
	double r = sqrt(rho*rho+z*z);
	Position path = post.position - pre.position;
	double angle_O_pre_post = acos(-pre.position.Unit().Dot(path.Unit()));
	double angle_O_post_pre = acos(post.position.Unit().Dot(path.Unit()));
	double angle_O_P_pre = asin(pre.position.R()*sin(angle_O_pre_post)/r);
	if(angle_O_pre_post>dd4hep::halfpi||angle_O_post_pre>dd4hep::halfpi){
	  bool backward = angle_O_post_pre>angle_O_pre_post;
	  double angle_O_P_pre = backward ? dd4hep::pi - asin(pre.position.R()*sin(angle_O_pre_post)/r) : asin(pre.position.R()*sin(angle_O_pre_post)/r);
	  double pre2P = r/sin(angle_O_pre_post)*sin(angle_O_pre_post+angle_O_P_pre);
	  pos = pre.position + pre2P*path.Unit();
	}
	Momentum mom = 0.5 * (pre.momentum + post.momentum);
	Geant4TrackerHit* hit = new Geant4TrackerHit(pre.truth.trackID,
						     pre.truth.pdgID,
						     pre.truth.deposit,
						     pre.truth.time);
	hit->cellID   = cellID;
	hit->position = pos;
	hit->momentum = mom;
	hit->length   = path.R();;
	clear();
	c->insert(hit);
	return hit;
      }
    };
  }
}

class TrackerCombineSensitiveDetector: public DDG4SensitiveDetector {
 public:
  TrackerCombineSensitiveDetector(const std::string& name, dd4hep::Detector& description);
  
  void Initialize(G4HCofThisEvent* HCE);
  G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
  void EndOfEvent(G4HCofThisEvent* HCE);
  
 protected:
  HitCollection* m_hc = nullptr;

 private:
  dd4hep::sim::TrackerCombine userData;
};
#endif
