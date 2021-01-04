#include "GenfitMsg.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/IMessageSvc.h"

MsgStream* GenfitMsg::m_log = nullptr;

MsgStream GenfitMsg::get(){
  if(nullptr == m_log){
    SmartIF<IMessageSvc> msgSvc(
        Gaudi::svcLocator()->service<IMessageSvc>("MessageSvc"));
    m_log = new MsgStream(msgSvc, "GenfitMsg");
    (*m_log) << MSG::DEBUG << "initialize GenfitMsg" << endmsg;
  }
  return (*m_log);
}

GenfitMsg::~GenfitMsg(){
  delete m_log;
}
