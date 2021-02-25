#ifndef GenEvent_h 
#define GenEvent_h 1

#include "edm4hep/MCParticleCollection.h"//plico

namespace MyHepMC {

class GenEvent{
    public: 
        //GenEvent();
        GenEvent(edm4hep::MCParticleCollection& mcCol);
        ~GenEvent();
        void SetEventHeader(long event_id_, long run_id_, float time_, const std::string& det_name_);
        //void SetMCCollection(edm4hep::MCParticleCollection vec_);
        long getID();
        long getRun();
        long getTime();
        void ReSet();
        std::string getName();
        edm4hep::MCParticleCollection& getMCVec();
        edm4hep::MCParticleCollection& m_mc_vec;
        //edm4hep::MCParticleCollection m_mc_vec;
    private:

        long m_event_id;
        long m_run_id;
        float m_time;
        std::string m_det_name;
};

}
#endif
