#ifndef GenEvent_h 
#define GenEvent_h 1

#include "plcio/MCParticleCollection.h"//plico

namespace MyHepMC {

class GenEvent{
    public: 
        //GenEvent();
        GenEvent(plcio::MCParticleCollection& mcCol);
        ~GenEvent();
        void SetEventHeader(long event_id_, long run_id_, float time_, std::string det_name_);
        //void SetMCCollection(plcio::MCParticleCollection vec_);
        long getID();
        long getRun();
        long getTime();
        void ReSet();
        std::string getName();
        plcio::MCParticleCollection getMCVec();
        plcio::MCParticleCollection& m_mc_vec;
        //plcio::MCParticleCollection m_mc_vec;
    private:
        long m_event_id;
        long m_run_id;
        float m_time;
        std::string m_det_name;
};

}
#endif
