#ifndef TrackSystemSvc_h
#define TrackSystemSvc_h

#include "TrackSystemSvc/ITrackSystemSvc.h"
#include <GaudiKernel/Service.h>

class TrackSystemSvc : public extends<Service, ITrackSystemSvc>{
 public:
  TrackSystemSvc(const std::string& name, ISvcLocator* svc);
  ~TrackSystemSvc();

  MarlinTrk::IMarlinTrkSystem* getTrackSystem(void* address=0) override;
  void removeTrackSystem(void* address=0) override;

  StatusCode initialize() override;
  StatusCode finalize() override;

 private:
  std::map<void*, MarlinTrk::IMarlinTrkSystem*> m_trackSystems;
};

#endif
