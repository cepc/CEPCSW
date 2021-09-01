#include "TrackHeedSimTool.h"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include <G4VProcess.hh>
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include <G4VTouchable.hh>
#include "DDG4/Geant4SensitiveDetector.h"
#include "DDG4/Geant4Converter.h"
#include "DDG4/Geant4Hits.h"
#include "DetSegmentation/GridDriftChamber.h"


#include <math.h>
#include <cmath>
#include <iostream>
#include <time.h>


DECLARE_COMPONENT(TrackHeedSimTool)

double TrackHeedSimTool::dedx(const edm4hep::MCParticle& mc) {return 0;}
double TrackHeedSimTool::dndx(double betagamma) {return 0;}

double TrackHeedSimTool::dedx(const G4Step* Step)
{

    clock_t t0 = clock();
    double de = 0;
    float cm_to_mm = 10; 
    G4Step* aStep = const_cast<G4Step*>(Step);
    G4Track* g4Track  =  aStep->GetTrack();//FIXME
    int pdg_code        = g4Track->GetParticleDefinition()->GetPDGEncoding();
    G4double pdg_mass   = g4Track->GetParticleDefinition()->GetPDGMass();
    G4double pdg_charge = g4Track->GetParticleDefinition()->GetPDGCharge();
    const G4VProcess* creatorProcess = g4Track->GetCreatorProcess();
    const G4String tmp_str_pro = (creatorProcess !=0) ? creatorProcess->GetProcessName() : "normal";
    G4double gammabeta=aStep->GetPreStepPoint()->GetBeta() * aStep->GetPreStepPoint()->GetGamma();
    if(gammabeta<0.01)return 1e-6;//too low momentum
    if(m_only_primary.value() && g4Track->GetParentID() != 0) return 1e-6;//give very small
    if(g4Track->GetParticleDefinition()->GetPDGCharge() ==0) return 0;//skip neutral particle 
    if(g4Track->GetKineticEnergy() <=0) return 0;//skip 
    if(pdg_code == 11 && (tmp_str_pro=="phot" || tmp_str_pro=="hIoni" || tmp_str_pro=="eIoni" || tmp_str_pro=="muIoni" || tmp_str_pro=="ionIoni" ) ) return 1e-6;//skip the electron produced by Ioni, because it is already simulated by TrackHeed
    if(m_particle_map.find(pdg_code) == m_particle_map.end() ) return 0;
    edm4hep::SimTrackerHitCollection* SimHitCol = nullptr;
    edm4hep::MCParticleCollection* mcCol = nullptr;
    try{
        SimHitCol =  const_cast<edm4hep::SimTrackerHitCollection*>(m_handle.get());
        mcCol = const_cast<edm4hep::MCParticleCollection*>(m_mc_handle.get());
    }
    catch(...){
        G4cout<<"can't find "<<m_handle.objKey() <<"in event, please check it have been createAndPut() in Begin of event"<<", or MCParicleCollection?"<<G4endl;
        throw "stop here!";
    }

    G4double track_KE   = aStep->GetPreStepPoint()->GetKineticEnergy();
    G4double track_time = aStep->GetPreStepPoint()->GetGlobalTime();
    G4double track_dx     = aStep->GetPreStepPoint()->GetMomentumDirection ().x();
    G4double track_dy     = aStep->GetPreStepPoint()->GetMomentumDirection ().y();
    G4double track_dz     = aStep->GetPreStepPoint()->GetMomentumDirection ().z();
    G4double track_length = aStep->GetStepLength();
    G4double position_x   = aStep->GetPreStepPoint()->GetPosition().x();
    G4double position_y   = aStep->GetPreStepPoint()->GetPosition().y();
    G4double position_z   = aStep->GetPreStepPoint()->GetPosition().z();
    int track_ID = g4Track->GetTrackID();
    int Parent_ID = g4Track->GetParentID();
    bool update_ke = true;
    if(m_use_max_step.value()){
        bool do_sim = false;
        if(m_isFirst){
            m_pre_x  = aStep->GetPreStepPoint()->GetPosition().x();
            m_pre_y  = aStep->GetPreStepPoint()->GetPosition().y();
            m_pre_z  = aStep->GetPreStepPoint()->GetPosition().z();
            m_pre_dx = aStep->GetPreStepPoint()->GetMomentumDirection().x();
            m_pre_dy = aStep->GetPreStepPoint()->GetMomentumDirection().y();
            m_pre_dz = aStep->GetPreStepPoint()->GetMomentumDirection().z();
            m_pre_t  = aStep->GetPreStepPoint()->GetGlobalTime();
            m_post_point = aStep->GetPostStepPoint();
            m_total_range += track_length;
            m_current_track_ID = g4Track->GetTrackID();
            m_current_Parent_ID = g4Track->GetParentID();
            m_pdg_code = g4Track->GetParticleDefinition()->GetPDGEncoding(); 
            m_isFirst = false;    
            m_pa_KE =  aStep->GetPreStepPoint()->GetKineticEnergy();
        }
        else{
        
            if(g4Track->GetTrackID() != m_current_track_ID){
                do_sim = true;
                m_change_track = true;
                update_ke = false;
            }
            else{
                m_post_point = aStep->GetPostStepPoint();
                m_total_range += track_length;
            }
        }
        if(m_total_range/CLHEP::mm >= m_max_step.value()){
            do_sim = true;
        }
        if(do_sim){
            track_KE = m_pa_KE;
            pdg_code = m_pdg_code;
            track_length = m_total_range;
            track_time = m_pre_t;
            track_dx = m_pre_dx;
            track_dy = m_pre_dy;
            track_dz = m_pre_dz;
            position_x = m_pre_x;
            position_y = m_pre_y;
            position_z = m_pre_z;
            track_ID = m_current_track_ID;
            Parent_ID = m_current_Parent_ID;
            if(m_change_track){
                m_pre_x  = aStep->GetPreStepPoint()->GetPosition().x();
                m_pre_y  = aStep->GetPreStepPoint()->GetPosition().y();
                m_pre_z  = aStep->GetPreStepPoint()->GetPosition().z();
                m_pre_dx = aStep->GetPreStepPoint()->GetMomentumDirection().x();
                m_pre_dy = aStep->GetPreStepPoint()->GetMomentumDirection().y();
                m_pre_dz = aStep->GetPreStepPoint()->GetMomentumDirection().z();
                m_pre_t  = aStep->GetPreStepPoint()->GetGlobalTime();
                m_post_point = aStep->GetPostStepPoint(); 
                m_total_range = aStep->GetStepLength();
                m_current_track_ID = g4Track->GetTrackID();
                m_current_Parent_ID = g4Track->GetParentID();
                m_pdg_code = g4Track->GetParticleDefinition()->GetPDGEncoding(); 
                m_change_track = false;
            }
            else{
                m_pre_x  = aStep->GetPostStepPoint()->GetPosition().x();
                m_pre_y  = aStep->GetPostStepPoint()->GetPosition().y();
                m_pre_z  = aStep->GetPostStepPoint()->GetPosition().z();
                m_pre_dx = aStep->GetPostStepPoint()->GetMomentumDirection().x();
                m_pre_dy = aStep->GetPostStepPoint()->GetMomentumDirection().y();
                m_pre_dz = aStep->GetPostStepPoint()->GetMomentumDirection().z();
                m_pre_t  = aStep->GetPostStepPoint()->GetGlobalTime();
                m_total_range = 0;
            }
        }
        else return 1e-6;
    }

    float init_x = 10;//cm
    float init_y = -10;//cm
    if(pdg_code == 11 && track_KE/CLHEP::keV < m_delta_threshold.value()){
        int nc = 0, ni=0;
        m_track->TransportDeltaElectron(init_x, init_y, 0, track_time/CLHEP::ns, track_KE/CLHEP::eV, track_dx, track_dy, track_dz, nc, ni);
        for (int j = 0; j < nc; ++j) {
            double xe = 0., ye = 0., ze = 0., te = 0., ee = 0.;
            double dx = 0., dy = 0., dz = 0.;
            m_track->GetElectron(j, xe, ye, ze, te, ee, dx, dy, dz);
            auto ehit = SimHitCol->create();
            ehit.setTime(te);
            double epos[3] = { cm_to_mm*( (xe - init_x)+position_x/CLHEP::cm) , cm_to_mm*((ye - init_y)+position_y/CLHEP::cm), cm_to_mm*(ze + position_z/CLHEP::cm)};
            ehit.setPosition(edm4hep::Vector3d(epos));
            /*// no sense of conductor electron
            float emom[3] = {0,0,0};
            getMom(ee, dx, dy, dz, emom);
            ehit.setMomentum(edm4hep::Vector3f(emom));
            */
            ehit.setQuality(3);
        }
        g4Track->SetTrackStatus(fStopAndKill);
        return 0;
    }
    clock_t t01 = clock();
    //cmp.SetMagneticField(0., 0., -3.);
    m_track->SetParticle(m_particle_map[pdg_code]);
    
    bool change_KE = false;
    if( abs(m_previous_KE -(track_KE/CLHEP::eV) )/(track_KE/CLHEP::eV) > m_change_threshold.value()) {
        change_KE = true;
        m_previous_KE = track_KE/CLHEP::eV;
    }
    
    bool change_ID = false;
    if(m_previous_track_ID != track_ID ){
        change_ID = true;
        m_previous_track_ID = track_ID; 
        m_previous_KE = track_KE/CLHEP::eV;
    }
    if(change_ID || change_KE ){
        m_track->SetKineticEnergy(track_KE/CLHEP::eV);
    }
    
  
    m_track->EnableOneStepFly(true);
    m_track->SetSteppingLimits( track_length/CLHEP::cm, 1000, 0.1, 0.2);
    clock_t t012 = clock();
    m_track->NewTrack(init_x, init_y, 0, track_time/CLHEP::ns, track_dx, track_dy, track_dz);//cm
    double xc = 0., yc = 0., zc = 0., tc = 0., ec = 0., extra = 0.;
    int nc = 0;
    int ic = 0;
    int first=true;
    clock_t t02 = clock();
    while (m_track->GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
        auto chit = SimHitCol->create();
        chit.setTime(tc);
        double cpos[3] = { cm_to_mm*( (xc - init_x)+position_x/CLHEP::cm) , cm_to_mm*((yc - init_y)+position_y/CLHEP::cm), cm_to_mm*(zc + position_z/CLHEP::cm)};
        chit.setPosition(edm4hep::Vector3d(cpos));
        float cmom[3]  = {0,0,0};
        getMom(ec, 1, 0, 0, cmom);//FIXME direction is not important?
        chit.setMomentum(edm4hep::Vector3f(cmom));
        chit.setQuality(1);
        if(m_save_cellID) chit.setCellID( getCellID(cpos[0], cpos[1], cpos[2]) );
        if(m_save_mc && Parent_ID == 0 && track_ID <= mcCol->size() ) chit.setMCParticle(  mcCol->at(track_ID-1) );
        if(first){
        chit.setPathLength( track_length/CLHEP::mm ); 
        first = false;
        }
        de += ec;
        for (int j = 0; j < nc; ++j) {
            double xe = 0., ye = 0., ze = 0., te = 0., ee = 0.;
            double dx = 0., dy = 0., dz = 0.;
            m_track->GetElectron(j, xe, ye, ze, te, ee, dx, dy, dz);
            auto ehit = SimHitCol->create();
            ehit.setTime(te);
            double epos[3] = { cm_to_mm*( (xe - init_x)+position_x/CLHEP::cm) , cm_to_mm*((ye - init_y)+position_y/CLHEP::cm), cm_to_mm*(ze + position_z/CLHEP::cm)};
            ehit.setPosition(edm4hep::Vector3d(epos));
            if(m_save_mc && Parent_ID == 0 && track_ID <= mcCol->size() ) chit.setMCParticle(  mcCol->at(track_ID-1) );
            /* //no sense of conductor electron
            float emom[3] = {0,0,0};
            getMom(ee, dx, dy, dz, emom); 
            ehit.setMomentum(edm4hep::Vector3f(emom));
            */
            if(m_save_cellID) ehit.setCellID( getCellID(epos[0], epos[1], epos[2]) );
            ehit.setQuality(2);
        }
    }
    double Dedx = (de*1e-6/(track_length/CLHEP::cm) ) ;//MeV/cm
    double new_KE = track_KE/CLHEP::MeV - de*1e-6; 
    if( update_ke ){
        g4Track->SetKineticEnergy ( new_KE*CLHEP::MeV );
        aStep->GetPostStepPoint()->SetKineticEnergy ( new_KE*CLHEP::MeV );
        m_pa_KE = new_KE;
    }
    else{
        m_pa_KE = aStep->GetPreStepPoint()->GetKineticEnergy(); 
    }
    m_tot_edep += de;
    return Dedx;
}


long long TrackHeedSimTool::getCellID(float x, float y, float z)
{
    float MM_2_CM = 0.1;
    G4Navigator* gNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    G4ThreeVector global(x,y,z);
    dd4hep::sim::Geant4VolumeManager volMgr = dd4hep::sim::Geant4Mapping::instance().volumeManager();
    G4VPhysicalVolume* pv = gNavigator->LocateGlobalPointAndSetup( global, 0, true);
    if(!pv) return 0;
    G4TouchableHistory *hist = gNavigator->CreateTouchableHistory();
    dd4hep::VolumeID volID  = volMgr.volumeID(hist);
    const G4AffineTransform & affine = gNavigator->GetGlobalToLocalTransform();
    G4ThreeVector local = affine.TransformPoint(global);
    dd4hep::Position loc(local.x()*MM_2_CM, local.y()*MM_2_CM, local.z()*MM_2_CM);
    dd4hep::Position glob(global.x()*MM_2_CM, global.y()*MM_2_CM, global.z()*MM_2_CM);
    dd4hep::VolumeID cID = m_segmentation->cellID(loc,glob,volID);
    
    if(m_debug){
        TVector3 Wstart(0,0,0);
        TVector3 Wend  (0,0,0);
        m_segmentation->cellposition(cID, Wstart, Wend);
        std::cout<<"Name="<<pv->GetName()<<",CopyNo="<<pv->GetCopyNo()<<",cID="<<cID<<",volID="<<volID<<",glob="<<glob<<",loc="<<loc<<",ws_x="<<Wstart.X()<<",y="<<Wstart.Y()<<",z="<<Wstart.Z()<<",we_x="<<Wend.X()<<",y="<<Wend.Y()<<",z="<<Wend.Z()<<std::endl;
    }
    delete hist;
    return cID;
}

void TrackHeedSimTool::getMom(float ee, float dx, float dy,float dz, float mom[3])
{
    double tot_E = 0.511*1e6 + ee;//eV
    double Mom = sqrt(tot_E*tot_E - pow(0.511*1e6,2) );
    double mom_direction =  sqrt(dx*dx + dy*dy + dz*dz);
    if (mom_direction == 0){
        mom[0] = 0;
        mom[1] = 0;
        mom[2] = 0;
    }
    else{
        double scale = Mom/mom_direction;
        mom[0] = scale*dx/1e9;
        mom[1] = scale*dy/1e9;
        mom[2] = scale*dz/1e9;
    }
}

StatusCode TrackHeedSimTool::initialize()
{

  m_geosvc = service<IGeomSvc>("GeomSvc");
  if ( !m_geosvc )  throw "TrackHeedSimTool :Failed to find GeomSvc ...";
  m_dd4hep = m_geosvc->lcdd();
  if ( !m_dd4hep )  throw "TrackHeedSimTool :Failed to get dd4hep::Detector ...";
  m_readout = new dd4hep::Readout( m_dd4hep->readout(m_readout_name) );
  if ( !m_readout )  throw "TrackHeedSimTool :Failed to get readout ...";
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::GridDriftChamber*>(m_readout->segmentation().segmentation());
  if ( !m_segmentation )  throw "TrackHeedSimTool :Failed to get segmentation ...";

  m_particle_map[ 11] = "e-";
  m_particle_map[-11] = "e+";
  m_particle_map[ 13] = "mu-";
  m_particle_map[-13] = "mu+";
  m_particle_map[ 211] = "pi+";
  m_particle_map[-211] = "pi-";
  m_particle_map[ 321] = "K+";
  m_particle_map[-321] = "K-";
  m_particle_map[2212] = "p";
  m_particle_map[-2212] = "pbar";
  m_particle_map[700201] = "d";
  m_particle_map[700202] = "alpha";

  m_gas.SetComposition("he", m_he,"isobutane", m_isob);
  m_gas.SetTemperature(293.15);
  m_gas.SetPressure(760.0);
  m_gas.SetMaxElectronEnergy(200.);
  m_gas.EnablePenningTransfer(0.44, 0.0, "He");
  m_gas.LoadGasFile(m_gas_file.value());
  m_gas.LoadIonMobility(m_IonMobility.value());
  //std::this_thread::sleep_for(std::chrono::milliseconds(m_delay_time));
  //m_gas.LoadGasFile("/junofs/users/wxfang/MyGit/tmp/check_G4FastSim_20210121/CEPCSW/Digitisers/DigiGarfield/He_50_isobutane_50.gas");
  //m_gas.LoadIonMobility("/junofs/users/wxfang/MyGit/tmp/check_G4FastSim_20210121/CEPCSW/Digitisers/DigiGarfield/IonMobility_He+_He.txt");
  /*
  m_gas.SetComposition("he", 90.,"isobutane", 10.);  // cepc gas
  m_gas.SetPressure(760.0);
  m_gas.SetTemperature(293.15);
  m_gas.SetFieldGrid(100., 100000., 20, true);
  m_gas.GenerateGasTable(10);
  m_gas.WriteGasFile("he_90_isobutane_10.gas");
  */

  cmp.SetMedium(&m_gas);
  // Field Wire radius [cm]
  const double rFWire = 110.e-4;
  // Signa Wire radius [cm]
  const double rSWire = 25.e-4;
  // Cell radius [cm]
  float rCell = 50;
  // Voltages
  const double vSWire = 2000.;
  const double vFWire = 0.;
  // Add the signal wire in the centre.
  cmp.AddWire(0, 0, 2 * rSWire, vSWire, "s");
  // Add the field wire around the signal wire.
  cmp.AddWire(-rCell, -rCell, 2 * rFWire, vFWire, "f");
  cmp.AddWire(    0., -rCell, 2 * rFWire, vFWire, "f");
  cmp.AddWire( rCell, -rCell, 2 * rFWire, vFWire, "f");
  cmp.AddWire(-rCell,     0., 2 * rFWire, vFWire, "f");
  cmp.AddWire( rCell,     0., 2 * rFWire, vFWire, "f");
  cmp.AddWire(-rCell,  rCell, 2 * rFWire, vFWire, "f");
  cmp.AddWire(    0.,  rCell, 2 * rFWire, vFWire, "f");
  cmp.AddWire( rCell,  rCell, 2 * rFWire, vFWire, "f");
  if(m_BField !=0 ) cmp.SetMagneticField(0., 0., m_BField);
  cmp.AddReadout("s");

  
  ///
  /// Make a sensor.
  ///
  m_sensor = new Sensor(); 
  m_sensor->AddComponent(&cmp);
  m_sensor->AddElectrode(&cmp, "s");
  // Set the signal time window. [ns]
  const double tstep = 0.5;
  const double tmin = -0.5 * 0.5;
  const unsigned int nbins = 1000;
  m_sensor->SetTimeWindow(tmin, tstep, nbins);
  m_sensor->ClearSignal();



  m_track = new Garfield::TrackHeed();
  //track->EnableDebugging();
  m_track->SetSensor(m_sensor);
  m_track->EnableDeltaElectronTransport();
  //track->DisableDeltaElectronTransport();
  m_track->EnableMagneticField();
  m_track->EnableElectricField();//almost no effect here

  m_current_track_ID = 0;
  m_previous_track_ID =0;
  m_previous_KE = 0;
  m_current_Parent_ID = -1;
  m_change_track = false;
  m_total_range = 0;
  m_isFirst = false;
  m_tot_edep = 0;
  m_pa_KE =0;
  m_pdg_code = 0;
  m_pre_x  = 0;
  m_pre_y  = 0;
  m_pre_z  = 0;
  m_pre_dx = 0;
  m_pre_dy = 0;
  m_pre_dz = 0;
  m_pre_t  = 0;
    return StatusCode::SUCCESS;
}

StatusCode TrackHeedSimTool::finalize()
{
    if(m_debug)std::cout << "m_tot_edep="<<m_tot_edep<<" eV"<<std::endl;
    return StatusCode::SUCCESS;
}
