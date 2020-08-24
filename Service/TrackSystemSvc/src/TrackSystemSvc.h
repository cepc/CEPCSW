#ifndef TrackSystemSvc_h
#define TrackSystemSvc_h

#include "TrackSystemSvc/ITrackSystemSvc.h"
#include <GaudiKernel/Service.h>

class TrackSystemSvc : public extends<Service, ITrackSystemSvc>{
 public:
  TrackSystemSvc(const std::string& name, ISvcLocator* svc);
  ~TrackSystemSvc();

  MarlinTrk::IMarlinTrkSystem* getTrackSystem() override;
  void removeTrackSystem() override;

  StatusCode initialize() override;
  StatusCode finalize() override;

 private:
  MarlinTrk::IMarlinTrkSystem* m_trackSystem;
};

#endif
