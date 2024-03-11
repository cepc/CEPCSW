#include "DetSimInterface/CommonUserTrackInfo.hh"
#include <iostream>

CommonUserTrackInfo::CommonUserTrackInfo() {

}

CommonUserTrackInfo::~CommonUserTrackInfo() {

}

void CommonUserTrackInfo::Print() const {

}

bool CommonUserTrackInfo::setIdxEdm4hep(int idxEdm4hep) {
    m_idxEdm4hep = idxEdm4hep;
}

int CommonUserTrackInfo::idxEdm4hep() const {
    return m_idxEdm4hep;
}
