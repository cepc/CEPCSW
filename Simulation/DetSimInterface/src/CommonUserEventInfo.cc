#include "DetSimInterface/CommonUserEventInfo.hh"
#include <iostream>

CommonUserEventInfo::CommonUserEventInfo() {

}

CommonUserEventInfo::~CommonUserEventInfo() {

}

void
CommonUserEventInfo::Print() const {

}

bool
CommonUserEventInfo::setIdxG4Track2Edm4hep(int idxG4, int idxEdm4hep) {
    auto it = m_g4track_to_edm4hep.find(idxG4);

    // if already exists, return false
    if (it != m_g4track_to_edm4hep.end()) {
        return false;
    }

    m_g4track_to_edm4hep[idxG4] = idxEdm4hep;

    return true;
}

int
CommonUserEventInfo::idxG4Track2Edm4hep(int idxG4) const {
    int ret = -1;

    auto it = m_g4track_to_edm4hep.find(idxG4);

    // if found
    if (it != m_g4track_to_edm4hep.end()) {
        ret = it->second;
    }

    return ret;
}

void
CommonUserEventInfo::dumpIdxG4Track2Edm4hep() const {
    std::cout << "---- Dumping IdxG4Track2Edm4hep: "
              << " total size: " << m_g4track_to_edm4hep.size()
              << std::endl;
    for (auto [idxG4, idxE4]: m_g4track_to_edm4hep) {
        std::cout << " - " << idxG4 << " -> " << idxE4 << std::endl;
    }
}
