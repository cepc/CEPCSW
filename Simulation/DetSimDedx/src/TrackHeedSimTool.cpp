#include "TrackHeedSimTool.h"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include <G4VProcess.hh>
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include <G4VTouchable.hh>
#include "DDG4/Geant4Converter.h"
#include "DetSegmentation/GridDriftChamber.h"


#include <math.h>
#include <cmath>
#include <iostream>
#include <time.h>
#include "CLHEP/Random/RandGauss.h"


DECLARE_COMPONENT(TrackHeedSimTool)

double TrackHeedSimTool::dedx(const edm4hep::MCParticle& mc) {return 0;}
double TrackHeedSimTool::dndx(double betagamma) {return 0;}

double TrackHeedSimTool::dedx(const G4Step* Step)
{

    clock_t t0 = clock();
    double de = 0;
    float cm_to_mm = 10; 
    G4Step* aStep = const_cast<G4Step*>(Step);
    G4Track* g4Track  =  aStep->GetTrack();
    int pdg_code        = g4Track->GetParticleDefinition()->GetPDGEncoding();
    G4double pdg_mass   = g4Track->GetParticleDefinition()->GetPDGMass();
    G4double pdg_charge = g4Track->GetParticleDefinition()->GetPDGCharge();
    const G4VProcess* creatorProcess = g4Track->GetCreatorProcess();
    const G4String tmp_str_pro = (creatorProcess !=0) ? creatorProcess->GetProcessName() : "normal";
    G4double gammabeta=aStep->GetPreStepPoint()->GetBeta() * aStep->GetPreStepPoint()->GetGamma();
    if(g4Track->GetParticleDefinition()->GetPDGCharge() ==0) return 0;//skip neutral particle 
    if(gammabeta<0.01)return m_eps;//too low momentum
    if(m_only_primary.value() && g4Track->GetParentID() != 0) return m_eps;
    if(g4Track->GetKineticEnergy() <=0) return 0;
    if(pdg_code == 11 && (tmp_str_pro=="phot" || tmp_str_pro=="hIoni" || tmp_str_pro=="eIoni" || tmp_str_pro=="muIoni" || tmp_str_pro=="ionIoni" ) ) return m_eps;//skip the electron produced by Ioni, because it is already simulated by TrackHeed
    if(m_particle_map.find(pdg_code) == m_particle_map.end() ) return m_eps;
    edm4hep::SimPrimaryIonizationClusterCollection* SimPrimaryIonizationCol = nullptr;
    edm4hep::MCParticleCollection* mcCol = nullptr;
    try{
        SimPrimaryIonizationCol =  const_cast<edm4hep::SimPrimaryIonizationClusterCollection*>(m_SimPrimaryIonizationCol.get());
        mcCol = const_cast<edm4hep::MCParticleCollection*>(m_mc_handle.get());
    }
    catch(...){
        G4cout<<"Error! Can't find collection in event, please check it have been createAndPut() in Begin of event"<<G4endl;
        G4cout<<"SimPrimaryIonizationCol="<<SimPrimaryIonizationCol<<",mcCol="<<mcCol<<G4endl;
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
        else return m_eps;
    }

    float init_x = 10;//cm
    float init_y = -10;//cm
    /*
    if(pdg_code == 11 && track_KE/CLHEP::keV < m_delta_threshold.value()){
        int nc = 0, ni=0;
        m_track->TransportDeltaElectron(init_x, init_y, 0, track_time/CLHEP::ns, track_KE/CLHEP::eV, track_dx, track_dy, track_dz, nc, ni);
        for (int j = 0; j < nc; ++j) {
            double xe = 0., ye = 0., ze = 0., te = 0., ee = 0.;
            double dx = 0., dy = 0., dz = 0.;
            m_track->GetElectron(j, xe, ye, ze, te, ee, dx, dy, dz);
            auto ehit  = SimIonizationCol->create();
            ehit.setTime(te);
            double epos[3] = { cm_to_mm*( (xe - init_x)+position_x/CLHEP::cm) , cm_to_mm*((ye - init_y)+position_y/CLHEP::cm), cm_to_mm*(ze + position_z/CLHEP::cm)};
            ehit.setPosition(edm4hep::Vector3d(epos));
            ehit.setType(11);
        }
        g4Track->SetTrackStatus(fStopAndKill);
        return 0;
    }
    */
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
        //auto chit = SimHitCol->create();
        auto chit = SimPrimaryIonizationCol->create();
        chit.setTime(tc);
        double cpos[3] = { cm_to_mm*( (xc - init_x)+position_x/CLHEP::cm) , cm_to_mm*((yc - init_y)+position_y/CLHEP::cm), cm_to_mm*(zc + position_z/CLHEP::cm)};
        chit.setPosition(edm4hep::Vector3d(cpos));
        //float cmom[3]  = {0,0,0};
        //getMom(ec, 1, 0, 0, cmom);//FIXME direction is not important?
        chit.setType(0);//default
        if(m_save_cellID) chit.setCellID( getCellID(cpos[0], cpos[1], cpos[2]) );
        if(m_save_mc && Parent_ID == 0 && track_ID <= mcCol->size() && mcCol ){ 
            chit.setMCParticle(  mcCol->at(track_ID-1) );
            //std::cout<<"mc obj index="<<mcCol->at(track_ID-1).getObjectID().index<<std::endl;
            //std::cout<<"mc obj index1="<<chit.getMCParticle().getObjectID().index<<std::endl;
        }
        de += ec;
        for (int j = 0; j < nc; ++j) {
            double xe = 0., ye = 0., ze = 0., te = 0., ee = 0.;
            double dx = 0., dy = 0., dz = 0.;
            m_track->GetElectron(j, xe, ye, ze, te, ee, dx, dy, dz);
            //auto ehit = SimHitCol->create();
            //auto ehit = SimIonizationCol->create();
            //ehit.setPrimaryIonization(chit);
            chit.addToElectronTime(te);
            //ehit.setTime(te);
            double epos[3] = { cm_to_mm*( (xe - init_x)+position_x/CLHEP::cm) , cm_to_mm*((ye - init_y)+position_y/CLHEP::cm), cm_to_mm*(ze + position_z/CLHEP::cm)};
            //ehit.setPosition(edm4hep::Vector3d(epos));
            //ehit.setPosition(edm4hep::Vector3d(epos));
            chit.addToElectronPosition(edm4hep::Vector3d(epos));
            //if(m_save_mc && Parent_ID == 0 && track_ID <= mcCol->size() && mcCol ){ 
            //    ehit.setMCParticle(  mcCol->at(track_ID-1) );
                //ehit.setMcParticleObjID( mcCol->at(track_ID-1).id() );
                //ehit.setMcParticleColID( mcCol->at(track_ID-1).getObjectID().collectionID );
            //}
            /* //no sense of conductor electron
            float emom[3] = {0,0,0};
            getMom(ee, dx, dy, dz, emom); 
            ehit.setMomentum(edm4hep::Vector3f(emom));
            */
            //if(m_save_cellID) ehit.setCellID( getCellID(epos[0], epos[1], epos[2]) );
            if(m_save_cellID) chit.addToElectronCellID( getCellID(epos[0], epos[1], epos[2]) );
            //ehit.setQuality(2);
            //ehit.setType(0);//default
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
    m_tot_length += track_length;
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
  float rCell = 50;//As the ionization process is almost not effected by cell geometry and wire voltage. Here the radius is to make sure the ionization process is completed.
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
  m_tot_length = 0;
  m_pa_KE =0;
  m_pdg_code = 0;
  m_pre_x  = 0;
  m_pre_y  = 0;
  m_pre_z  = 0;
  m_pre_dx = 0;
  m_pre_dy = 0;
  m_pre_dz = 0;
  m_pre_t  = 0;


   // for NN pulse simulation//
   m_env = std::make_shared<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "ENV");
   m_seesion_options = std::make_shared<Ort::SessionOptions>();
   m_seesion_options->SetIntraOpNumThreads(m_intra_op_nthreads);
   m_seesion_options->SetInterOpNumThreads(m_inter_op_nthreads);
   if(m_debug) std::cout << "before load model " << m_model_file.value() << std::endl;
   m_session = std::make_shared<Ort::Session>(*m_env, m_model_file.value().c_str(), *m_seesion_options);
   if(m_debug) std::cout << "after load model " << m_model_file.value() << std::endl;
   // lambda function to print the dims.
   auto dims_str = [&](const auto& dims) {
      return std::accumulate(dims.begin(), dims.end(), std::to_string(dims[0]),
                             [](const std::string& a, int64_t b){
                                 return a + "x" + std::to_string(b);
                             });
   };
   // prepare the input
   auto num_input_nodes = m_session->GetInputCount();
   if(m_debug) std::cout << "num_input_nodes: " << num_input_nodes << std::endl;
   for (size_t i = 0; i < num_input_nodes; ++i) {
      auto name = m_session->GetInputNameAllocated(i, m_allocator);
      m_inputNodeNameAllocatedStrings.push_back(std::move(name));
      m_input_node_names.push_back(m_inputNodeNameAllocatedStrings.back().get());

      Ort::TypeInfo type_info  = m_session->GetInputTypeInfo(i);
      auto tensor_info = type_info.GetTensorTypeAndShapeInfo();
      auto dims = tensor_info.GetShape();
      //dims[0] = 1; //wxfang, FIXME, if it is -1 (dynamic axis), need overwrite it manually
      dims[0] = 10; //wxfang, FIXME, if it is -1 (dynamic axis), need overwrite it manually
      m_input_node_dims.push_back(dims);


      if(m_debug) std::cout<< "[" << i << "]"
              << " input_name: " << m_inputNodeNameAllocatedStrings.back().get()
              << " ndims: " << dims.size()
              << " dims: " << dims_str(dims)
              << std::endl;
   }
   // prepare the output
   size_t num_output_nodes = m_session->GetOutputCount();
   for(std::size_t i = 0; i < num_output_nodes; i++) {
       auto output_name = m_session->GetOutputNameAllocated(i, m_allocator);
       m_outputNodeNameAllocatedStrings.push_back(std::move(output_name));
       m_output_node_names.push_back(m_outputNodeNameAllocatedStrings.back().get());
       Ort::TypeInfo type_info        = m_session->GetOutputTypeInfo(i);
       auto tensor_info               = type_info.GetTensorTypeAndShapeInfo();
       ONNXTensorElementDataType type = tensor_info.GetElementType();
       m_output_node_dims               = tensor_info.GetShape();
       if(m_debug) std::cout << "[" << i << "]"
               << " output_name: " << m_outputNodeNameAllocatedStrings.back().get()
               << " ndims: " << m_output_node_dims.size()
               << " dims: " << dims_str(m_output_node_dims)
               << std::endl;

   }

  return StatusCode::SUCCESS;
}

void TrackHeedSimTool::wire_xy(float x1, float y1, float z1, float x2, float y2, float z2, float z, float &x, float &y){
    //linear function: 
    //(x-x1)/(x2-x1)=(y-y1)/(y2-y1)=(z-z1)/(z2-z1)
    x = x1+(x2-x1)*(z-z1)/(z2-z1);
    y = y1+(y2-y1)*(z-z1)/(z2-z1);
}

float TrackHeedSimTool::xy2phi(float x, float y){
    float phi = acos(x/sqrt(x*x+y*y));
    if(y < 0) phi = 2*M_PI-phi;
    return phi; 
}
void TrackHeedSimTool::getLocal(float x1, float y1, float x2, float y2, float& dx, float& dy){
    /*  .    
        .     *(x2,y2)
        .   .
        . .    
        *(x1, y1)        
        .  
        .  
        . 
        o   
    */
    float mo1 = sqrt(x1*x1+y1*y1);
    float mo2 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1) );
    float costheta = (x1*(x2-x1)+y1*(y2-y1))/(mo1*mo2);
    dy = mo2*costheta;
    dx = xy2phi(x2,y2)>xy2phi(x1,y1) ? mo2*sqrt(1-costheta*costheta) : -mo2*sqrt(1-costheta*costheta) ;
}

float* TrackHeedSimTool::NNPred(std::vector<float>& inputs)
{


    std::vector<Ort::Value> input_tensors;
    auto& dims = m_input_node_dims[0];
    //std::cout << "inputs.size()="<<inputs.size() << std::endl;
    dims[0] = int(inputs.size()/3);
    Ort::MemoryInfo info("Cpu", OrtDeviceAllocator, 0, OrtMemTypeDefault);
    // prepare a dummy input for the model
    
    
    auto input_tensor = Ort::Value::CreateTensor(info,
                                                 inputs.data(),
                                                 inputs.size(),
                                                 dims.data(),
                                                 dims.size());
    
    input_tensors.push_back(std::move(input_tensor));
    auto output_tensors = m_session->Run(Ort::RunOptions{ nullptr }, m_input_node_names.data(), input_tensors.data(), input_tensors.size(), m_output_node_names.data(), m_output_node_names.size());
    //const auto& output_tensor = output_tensors[0];
    auto& output_tensor = output_tensors[0];
    int num_elements = output_tensor.GetTensorTypeAndShapeInfo().GetElementCount();
    //std::cout << "output_tensor num_elements=" << num_elements<< std::endl;

    float* vec2 = new float[num_elements];
    std::memcpy(vec2, output_tensor.GetTensorMutableData<float>(), num_elements * sizeof(float));
    /*
    for (int k=0;k<num_elements;k++){
       std::cout<<"k="<<k<< ",v=" <<vec2[k]<< std::endl;
    }
    */ 
    return vec2;
}


void TrackHeedSimTool::endOfEvent() {
    if(m_sim_pulse){
        edm4hep::SimPrimaryIonizationClusterCollection* SimPrimaryIonizationCol = nullptr;
        try{
            SimPrimaryIonizationCol =  const_cast<edm4hep::SimPrimaryIonizationClusterCollection*>(m_SimPrimaryIonizationCol.get());
        }
        catch(...){
            G4cout<<"Error! Can't find collection in event, please check it have been createAndPut() in Begin of event"<<G4endl;
            G4cout<<"SimPrimaryIonizationCol="<<SimPrimaryIonizationCol<<G4endl;
            throw "stop here!";
        }
        if(m_debug) G4cout<<"SimPrimaryIonizationCol size="<<SimPrimaryIonizationCol->size()<<G4endl;
        clock_t t01 = clock();
        std::vector<float> inputs;
        std::vector<unsigned long> indexs_c;
        std::vector<unsigned long> indexs_i;
        std::map<unsigned long long, std::vector<std::pair<float, float> > > id_pulse_map;
        for (unsigned long z=0; z<SimPrimaryIonizationCol->size(); z++) {
            for (unsigned long k=0; k<SimPrimaryIonizationCol->at(z).electronCellID_size(); k++) {
                //auto simIon = SimIonizationCol->at(k);
                auto position = SimPrimaryIonizationCol->at(z).getElectronPosition(k);//mm
                auto cellId = SimPrimaryIonizationCol->at(z).getElectronCellID(k);
                TVector3 Wstart(0,0,0);
                TVector3 Wend  (0,0,0);
                m_segmentation->cellposition(cellId, Wstart, Wend);
                float dd4hep_mm = dd4hep::mm;
                Wstart =(1/dd4hep_mm)* Wstart;// from DD4HEP cm to mm
                Wend   =(1/dd4hep_mm)* Wend  ;
                //std::cout<<"cellid="<<cellId<<",s_x="<<Wstart.x()<<",s_y="<<Wstart.y()<<",s_z="<<Wstart.z()<<",E_x="<<Wend.x()<<",E_y="<<Wend.y()<<",E_z="<<Wend.z()<<std::endl;
                float wire_x = 0;
                float wire_y = 0;
                double pos_z = position[2];
                wire_xy(Wend.x(), Wend.y(), Wend.z(), Wstart.x(), Wstart.y(), Wstart.z(), pos_z, wire_x, wire_y);
                float local_x = 0;
                float local_y = 0;
                getLocal(wire_x, wire_y, position[0], position[1], local_x, local_y);
                //std::cout<<"pos_z="<<pos_z<<",wire_x="<<wire_x<<",wire_y="<<wire_y<<",position[0]="<<position[0]<<",position[1]="<<position[1]<<",local_x="<<local_x<<",local_y="<<local_y<<",dr="<<sqrt(local_x*local_x+local_y*local_y)<<",dr1="<<sqrt( (wire_x-position[0])*(wire_x-position[0])+(wire_y-position[1])*(wire_y-position[1]) )<<std::endl;
                float m_x_scale = 1;
                float m_y_scale = 1;
                local_x = local_x/m_x_scale;//FIXME, default is 18mm x 18mm, the real cell size maybe a bit different. need konw size of cell for each layer, and do normalization
                local_y = local_y/m_y_scale;
                float noise = CLHEP::RandGauss::shoot(0,1);
                inputs.push_back(local_x);
                inputs.push_back(local_y);
                inputs.push_back(noise  );
                indexs_c.push_back(z);
                indexs_i.push_back(k);
                if(indexs_c.size()==m_batchsize){
                    float* res = NNPred(inputs);
                    for(unsigned int i=0; i<m_batchsize; i++){
                        float tmp_time = res[i*2  ]*m_time_scale + m_time_shift;// in ns
                        float tmp_amp  = res[i*2+1]*m_amp_scale  + m_amp_shift ;
                        //unsigned long tmp_index = indexs.at(i);
                        //tmp_pluse.setCellID(SimIonizationCol->at(tmp_index).getCellID());
                        //tmp_pluse.setTime(tmp_time + SimIonizationCol->at(tmp_index).getTime());//ns
                        //tmp_pluse.setValue(tmp_amp);
                        //tmp_pluse.setType(SimIonizationCol->at(tmp_index).getType());
                        //tmp_pluse.setSimIonization(SimIonizationCol->at(tmp_index));
                        auto ion_time = SimPrimaryIonizationCol->at(indexs_c.at(i)).getElectronTime(indexs_i.at(i));
                        id_pulse_map[indexs_c.at(i)].push_back(std::make_pair(tmp_time+ion_time,tmp_amp) );
                    }
                    inputs.clear();
                    indexs_c.clear();
                    indexs_i.clear();
                    delete [] res;
                }
            } //end of k
        }//end of z
        if(indexs_c.size()!=0){
            float* res = NNPred(inputs);
            for(unsigned int i=0; i<indexs_c.size(); i++){
                float tmp_time = res[i*2  ]*m_time_scale + m_time_shift;
                float tmp_amp  = res[i*2+1]*m_amp_scale  + m_amp_shift ;
                //tmp_pluse.setCellID(SimIonizationCol->at(tmp_index).getCellID());
                //tmp_pluse.setTime(tmp_time + SimIonizationCol->at(tmp_index).getTime());//ns
                //tmp_pluse.setValue(tmp_amp);
                //tmp_pluse.setType(SimIonizationCol->at(tmp_index).getType());
                //tmp_pluse.setSimIonization(SimIonizationCol->at(tmp_index));
                //id_pulse_map[SimIonizationCol->at(tmp_index).getCellID()].push_back(std::make_pair(tmp_pluse.getTime(), tmp_pluse.getValue() ) );
                auto ion_time = SimPrimaryIonizationCol->at(indexs_c.at(i)).getElectronTime(indexs_i.at(i));
                id_pulse_map[indexs_c.at(i)].push_back(std::make_pair(tmp_time+ion_time,tmp_amp) );
            }
            inputs.clear();
            indexs_c.clear();
            indexs_i.clear();
            delete [] res;
        }
        for(auto iter = id_pulse_map.begin(); iter != id_pulse_map.end(); iter++){
            edm4hep::MutableSimPrimaryIonizationCluster dcIonCls = SimPrimaryIonizationCol->at(iter->first); 
            for(unsigned int i=0; i< iter->second.size(); i++){
                auto tmp_time = iter->second.at(i).first ;
                auto tmp_amp  = iter->second.at(i).second;
                dcIonCls.addToPulseTime(tmp_time);  
                dcIonCls.addToPulseAmplitude(tmp_amp);  
            }
            if(dcIonCls.electronPosition_size() != dcIonCls.pulseTime_size()){
                G4cout<<"Error ion size != pulse size"<<G4endl;
                throw "stop here!";
            }
        }
        clock_t t02 = clock();
        if(m_debug) std::cout<<"time for Pulse Simulation=" << (double)(t02 - t01) / CLOCKS_PER_SEC <<" seconds"<< std::endl;
    }
}

StatusCode TrackHeedSimTool::finalize()
{
    //if(m_debug)std::cout << "m_tot_edep="<<m_tot_edep<<" eV"<<std::endl;
    //std::cout << "m_tot_length="<<m_tot_length<<" mm"<<std::endl;
    return StatusCode::SUCCESS;
}
