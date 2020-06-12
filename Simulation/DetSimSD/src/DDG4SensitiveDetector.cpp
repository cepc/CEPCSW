#include "DetSimSD/DDG4SensitiveDetector.h"

void
DDG4SensitiveDetector::Initialize(G4HCofThisEvent* HCE) {

}

G4bool
DDG4SensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) {
    
    return true;
}

void
DDG4SensitiveDetector::EndOfEvent(G4HCofThisEvent* HCE) {

}

long long
DDG4SensitiveDetector::getVolumeID(G4Step* step) {
    long long vid = 0;

    return vid;
}

long long
DDG4SensitiveDetector::getCellID(G4Step* step) {
    long long vid = 0;

    return vid;
}
