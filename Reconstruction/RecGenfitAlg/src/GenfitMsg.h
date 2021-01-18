//////////////////////////////////////////////////////////////////////
///
/// This is a class for genfit interfaces classes to access Gaudi log
///   system
///
/// Authors:
///   Yao ZHANG(zhangyao@ihep.ac.cn)
///
//////////////////////////////////////////////////////////////////////
#ifndef RECGENFITALG_GENFITMSG_H
#define RECGENFITALG_GENFITMSG_H

#include "GaudiKernel/MsgStream.h"

class GenfitMsg {
  public:
    static MsgStream get();
  private:
    GenfitMsg(){};
    virtual ~GenfitMsg();
    static MsgStream* m_log;

};

#endif

