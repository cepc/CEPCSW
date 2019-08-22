#ifndef FWCORE_CEPCDATASVC_H
#define FWCORE_CEPCDATASVC_H

#include "FWCore/PodioDataSvc.h"

class CEPCDataSvc : public PodioDataSvc
{
  friend class SvcFactory<CEPCDataSvc>;

public:
  /// Standard Constructor
  CEPCDataSvc(const std::string& name, ISvcLocator* svc);

  /// Standard Destructor
  virtual ~CEPCDataSvc();
};

#endif  // FWCORE_CEPCDATASVC_H
