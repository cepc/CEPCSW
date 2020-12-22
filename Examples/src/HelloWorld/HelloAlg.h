#ifndef HelloAlg_h
#define HelloAlg_h

#include <GaudiKernel/Algorithm.h>
#include <Gaudi/Property.h>

class HelloAlg: public Algorithm {
public:
    HelloAlg(const std::string& name, ISvcLocator* pSvcLocator);

    StatusCode initialize() override;
    StatusCode execute() override;
    StatusCode finalize() override;

private:

    Gaudi::Property<int> m_int{this, "MyInt", 42};
};


#endif
