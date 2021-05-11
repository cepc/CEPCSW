#ifndef CLUPATRA_TRACKINGHELPER
#define CLUPATRA_TRACKINGHELPER

#include "edm4hep/TrackState.h"
#include "edm4hep/Track.h"
#include "UTIL/BitField64.h"
#include "UTIL/CellIDDecoder.h"
#include "UTIL/ILDConf.h"
#include "lcio.h"
#include <array>

inline bool hasTrackStateAt(edm4hep::ConstTrack track, int location) {
    for (auto it = track.trackStates_begin(); it != track.trackStates_end(); it++) {
        if (it->location == location) {
            return true;
        }
    }
    return false;
}

inline edm4hep::TrackState getTrackStateAt(edm4hep::ConstTrack track, int location) {
    for (auto it = track.trackStates_begin(); it != track.trackStates_end(); it++) {
        if (it->location == location) {
	  return *it;
        }
    }
    return edm4hep::TrackState();
}

inline std::array<float,15> getCovMatrix(const edm4hep::ConstTrack &track) {
    return track.getTrackStates(0).covMatrix;
}
inline float getTanLambda(const edm4hep::ConstTrack &track) {
    return track.getTrackStates(0).tanLambda;
}
inline float getOmega(const edm4hep::ConstTrack &track) {
    return track.getTrackStates(0).omega;
}
inline float getD0(const edm4hep::ConstTrack &track) {
    return track.getTrackStates(0).D0;
}
inline float getZ0(const edm4hep::ConstTrack &track) {
    return track.getTrackStates(0).Z0;
}
inline float getPhi(const edm4hep::ConstTrack &track) {
    return track.getTrackStates(0).phi;
}


inline int getLayer(const edm4hep::TrackerHit hit) {
    UTIL::BitField64* _encoder = new UTIL::BitField64(lcio::ILDCellID0::encoder_string);
    _encoder->setValue(hit.getCellID());
    int layer = (*_encoder)[lcio::ILDCellID0::layer];
    delete _encoder;
    return layer;
}

inline int getLayer(const edm4hep::ConstTrackerHit hit) {
    UTIL::BitField64* _encoder = new UTIL::BitField64(lcio::ILDCellID0::encoder_string);
    _encoder->setValue(hit.getCellID());
    int layer = (*_encoder)[lcio::ILDCellID0::layer];
    delete _encoder;
    return layer;
}

#endif
