/**
 * 
 *  @brief  Implementation of the track creator class.
 * 
 *  $Log: $
 */


#include "UTIL/ILDConf.h"

#include "edm4hep/Vertex.h"
#include "edm4hep/ReconstructedParticle.h"
#if __has_include("edm4hep/EDM4hepVersion.h")
#include "edm4hep/EDM4hepVersion.h"
#else
// Copy the necessary parts from  the header above to make whatever we need to work here
#define EDM4HEP_VERSION(major, minor, patch) ((UINT64_C(major) << 32) | (UINT64_C(minor) << 16) | (UINT64_C(patch)))
// v00-09 is the last version without the capitalization change of the track vector members
#define EDM4HEP_BUILD_VERSION EDM4HEP_VERSION(0, 9, 0)
#endif

#include "gear/BField.h"
#include "gear/CalorimeterParameters.h"
#include "gear/PadRowLayout2D.h"
#include "gear/TPCParameters.h"
#include "gear/FTDParameters.h"
#include "gear/FTDLayerLayout.h"

#include "GaudiKernel/IService.h"
#include "GearSvc/IGearSvc.h"
#include "PandoraMatrixAlg.h"

#include "TrackCreator.h"
#include "Pandora/PdgTable.h"

#include <algorithm>
#include <cmath>
#include <limits>

TrackCreator::TrackCreator(const Settings &settings, const pandora::Pandora *const pPandora, ISvcLocator* svcloc) :
    m_settings(settings),
    m_pPandora(pPandora)
{

    IGearSvc*  iSvc = 0;
    StatusCode sc = svcloc->service("GearSvc", iSvc, false);
    if ( !sc ) 
    {
        throw "Failed to find GearSvc ...";
    }
    _GEAR = iSvc->getGearMgr();


    m_bField                  = (_GEAR->getBField().at(gear::Vector3D(0., 0., 0.)).z());
    m_tpcInnerR               = (_GEAR->getTPCParameters().getPadLayout().getPlaneExtent()[0]);
    m_tpcOuterR               = (_GEAR->getTPCParameters().getPadLayout().getPlaneExtent()[1]);
    m_tpcMaxRow               = (_GEAR->getTPCParameters().getPadLayout().getNRows());
    m_tpcZmax                 = (_GEAR->getTPCParameters().getMaxDriftLength());
    m_eCalBarrelInnerSymmetry = (_GEAR->getEcalBarrelParameters().getSymmetryOrder());
    m_eCalBarrelInnerPhi0     = (_GEAR->getEcalBarrelParameters().getPhi0());
    m_eCalBarrelInnerR        = (_GEAR->getEcalBarrelParameters().getExtent()[0]);
    m_eCalEndCapInnerZ        = (_GEAR->getEcalEndcapParameters().getExtent()[2]);
    // fg: FTD description in GEAR has changed ...
    try
    {
        m_ftdInnerRadii = _GEAR->getGearParameters("FTD").getDoubleVals("FTDInnerRadius");
        m_ftdOuterRadii = _GEAR->getGearParameters("FTD").getDoubleVals("FTDOuterRadius");
        m_ftdZPositions = _GEAR->getGearParameters("FTD").getDoubleVals("FTDZCoordinate");
        m_nFtdLayers = m_ftdZPositions.size();
    }
    catch (gear::UnknownParameterException &)
    {
        const gear::FTDLayerLayout &ftdLayerLayout(_GEAR->getFTDParameters().getFTDLayerLayout());
        std::cout << " Filling FTD parameters from gear::FTDParameters - n layers: " << ftdLayerLayout.getNLayers() << std::endl;

        for(unsigned int i = 0, N = ftdLayerLayout.getNLayers(); i < N; ++i)
        {
            // Create a disk to represent even number petals front side
            m_ftdInnerRadii.push_back(ftdLayerLayout.getSensitiveRinner(i));
            m_ftdOuterRadii.push_back(ftdLayerLayout.getMaxRadius(i));

            // Take the mean z position of the staggered petals
            const double zpos(ftdLayerLayout.getZposition(i));
            m_ftdZPositions.push_back(zpos);
            std::cout << "     layer " << i << " - mean z position = " << zpos << std::endl;
        }

        m_nFtdLayers = m_ftdZPositions.size() ;
    }

    // Check tpc parameters
    if ((std::fabs(m_tpcZmax) < std::numeric_limits<float>::epsilon()) || (std::fabs(m_tpcInnerR) < std::numeric_limits<float>::epsilon())
        || (std::fabs(m_tpcOuterR - m_tpcInnerR) < std::numeric_limits<float>::epsilon()))
    {
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    m_cosTpc = m_tpcZmax / std::sqrt(m_tpcZmax * m_tpcZmax + m_tpcInnerR * m_tpcInnerR);

    // Check ftd parameters
    if ((0 == m_nFtdLayers) || (m_nFtdLayers != m_ftdInnerRadii.size()) || (m_nFtdLayers != m_ftdOuterRadii.size()))
    {
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }

    for (unsigned int iFtdLayer = 0; iFtdLayer < m_nFtdLayers; ++iFtdLayer)
    {
        if ((std::fabs(m_ftdOuterRadii[iFtdLayer]) < std::numeric_limits<float>::epsilon()) ||
            (std::fabs(m_ftdInnerRadii[iFtdLayer]) < std::numeric_limits<float>::epsilon()))
        {
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
        }
    }

    m_tanLambdaFtd = m_ftdZPositions[0] / m_ftdOuterRadii[0];

    // Calculate etd and set parameters
    // fg: make SET and ETD optional - as they might not be in the model ...
    try
    {
        const DoubleVector &etdZPositions(_GEAR->getGearParameters("ETD").getDoubleVals("ETDLayerZ"));
        const DoubleVector &setInnerRadii(_GEAR->getGearParameters("SET").getDoubleVals("SETLayerRadius"));

        if (etdZPositions.empty() || setInnerRadii.empty())
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

        m_minEtdZPosition = *(std::min_element(etdZPositions.begin(), etdZPositions.end()));
        m_minSetRadius = *(std::min_element(setInnerRadii.begin(), setInnerRadii.end()));
    }
    catch(gear::UnknownParameterException &)
    {
        std::cout << "Warnning, ETDLayerZ or SETLayerRadius parameters missing from GEAR parameters!" << std::endl
                               << "     -> both will be set to " << std::numeric_limits<float>::quiet_NaN() << std::endl;

        //fg: Set them to NAN, so that they cannot be used to set   trackParameters.m_reachesCalorimeter = true;
        m_minEtdZPosition = std::numeric_limits<float>::quiet_NaN();
        m_minSetRadius = std::numeric_limits<float>::quiet_NaN();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackCreator::~TrackCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode TrackCreator::CreateTrackAssociations(const CollectionMaps& collectionMaps)
{
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractKinks(collectionMaps));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractProngsAndSplits(collectionMaps));
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->ExtractV0s(collectionMaps));

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
pandora::StatusCode TrackCreator::ExtractKinks(const CollectionMaps& collectionMaps)
{
    std::cout<<"start TrackCreator::ExtractKinks:"<<std::endl;
    for (StringVector::const_iterator iter = m_settings.m_kinkVertexCollections.begin(), iterEnd = m_settings.m_kinkVertexCollections.end();
        iter != iterEnd; ++iter)
    {
        if(collectionMaps.collectionMap_Vertex.find(*iter) == collectionMaps.collectionMap_Vertex.end()) { std::cout<<"not find "<<(*iter)<<std::endl; continue;}
        try
        {
            auto pKinkCollection = (collectionMaps.collectionMap_Vertex.find(*iter))->second;

            for (int i = 0, iMax = pKinkCollection.size(); i < iMax; ++i)
            {
                try
                {
                    auto  pVertex0 = pKinkCollection.at(i);
                    auto  pVertex  = &(pVertex0);

                    if (NULL == pVertex) throw ("Collection type mismatch");

                    auto pReconstructedParticle = pVertex->getAssociatedParticle();
                    if (this->IsConflictingRelationship(pReconstructedParticle))continue;

                    const int vertexPdgCode(pReconstructedParticle.getType());

                    // Extract the kink vertex information
                    for (unsigned int iTrack = 0, nTracks = pReconstructedParticle.tracks_size(); iTrack < nTracks; ++iTrack)
                    {
                        auto pTrack = pReconstructedParticle.getTracks(iTrack);
                        (0 == iTrack) ? m_parentTrackList.insert(pTrack.id().index) : m_daughterTrackList.insert(pTrack.id().index);

                        int trackPdgCode = pandora::UNKNOWN_PARTICLE_TYPE;

                        if (0 == iTrack)
                        {
                            trackPdgCode = vertexPdgCode;
                        }
                        else
                        {
                            switch (vertexPdgCode)
                            {
                            case pandora::PI_PLUS :
                            case pandora::K_PLUS :
                                trackPdgCode = pandora::MU_PLUS;
                                break;
                            case pandora::PI_MINUS :
                            case pandora::K_MINUS :
                                trackPdgCode = pandora::MU_MINUS;
                                break;
                            case pandora::HYPERON_MINUS_BAR :
                            case pandora::SIGMA_PLUS :
                                trackPdgCode = pandora::PI_PLUS;
                                break;
                            case pandora::SIGMA_MINUS :
                            case pandora::HYPERON_MINUS :
                                trackPdgCode = pandora::PI_PLUS;
                                break;
                            default :
                                (pTrack.getTrackStates(0).omega > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PI_MINUS;
                                break;
                            }
                        }

                        m_trackToPidMap.insert(TrackToPidMap::value_type(pTrack, trackPdgCode));

                        if (0 == m_settings.m_shouldFormTrackRelationships)
                            continue;

                        // Make track parent-daughter relationships
                        if (0 == iTrack)
                        {
                            for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack)
                            {
                                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackParentDaughterRelationship(*m_pPandora, GetTrackAddress(collectionMaps, pTrack), GetTrackAddress(collectionMaps, pReconstructedParticle.getTracks(jTrack) ) ) );
                            }
                        }

                        // Make track sibling relationships
                        else
                        {
                            for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack)
                            {
                                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackSiblingRelationship(*m_pPandora, GetTrackAddress(collectionMaps, pTrack), GetTrackAddress(collectionMaps, pReconstructedParticle.getTracks(jTrack) ) ) );
                            }
                        }
                    }
                }
                catch (...)
                {
                    std::cout << "Failed to extract kink vertex: " <<  std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout << "Failed to extract kink vertex collection: " << *iter <<  std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode TrackCreator::ExtractProngsAndSplits(const CollectionMaps& collectionMaps)
{
    std::cout<<"start TrackCreator::ExtractProngsAndSplits:"<<std::endl;
    for (StringVector::const_iterator iter = m_settings.m_prongSplitVertexCollections.begin(), iterEnd = m_settings.m_prongSplitVertexCollections.end(); iter != iterEnd; ++iter)
    {
        if(collectionMaps.collectionMap_Vertex.find(*iter) == collectionMaps.collectionMap_Vertex.end()) { std::cout<<"not find "<<(*iter)<<std::endl; continue;}
        try
        {
            auto pProngOrSplitCollection = (collectionMaps.collectionMap_Vertex.find(*iter))->second;

            for (int i = 0, iMax = pProngOrSplitCollection.size(); i < iMax; ++i)
            {
                try
                {
                    auto  pVertex0 = pProngOrSplitCollection.at(i);
                    auto  pVertex  = &(pVertex0);

                    if (NULL == pVertex) throw ("Collection type mismatch");
                    auto pReconstructedParticle = pVertex->getAssociatedParticle();

                    if (this->IsConflictingRelationship(pReconstructedParticle))continue;

                    // Extract the prong/split vertex information
                    for (unsigned int iTrack = 0, nTracks = pReconstructedParticle.tracks_size(); iTrack < nTracks; ++iTrack)
                    {
                        auto pTrack = pReconstructedParticle.getTracks(iTrack);
                        (0 == iTrack) ? m_parentTrackList.insert(pTrack.id().index) : m_daughterTrackList.insert(pTrack.id().index);

                        if (0 == m_settings.m_shouldFormTrackRelationships) continue;

                        // Make track parent-daughter relationships
                        if (0 == iTrack)
                        {
                            for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack)
                            {
                                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackParentDaughterRelationship(*m_pPandora, GetTrackAddress(collectionMaps, pTrack), GetTrackAddress(collectionMaps, pReconstructedParticle.getTracks(jTrack) ) ) );
                            }
                        }

                        // Make track sibling relationships
                        else
                        {
                            for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack)
                            {
                                PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackSiblingRelationship(*m_pPandora, GetTrackAddress(collectionMaps, pTrack), GetTrackAddress(collectionMaps, pReconstructedParticle.getTracks(jTrack) ) ) );
                            }
                        }
                    }
                }
                catch (...)
                {
                    std::cout << "Failed to extract prong/split vertex" <<std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout<< "Failed to extract prong/split vertex collection: " << *iter << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode TrackCreator::ExtractV0s(const CollectionMaps& collectionMaps)
{
    std::cout<<"start TrackCreator::ExtractV0s:"<<std::endl;
    for (StringVector::const_iterator iter = m_settings.m_v0VertexCollections.begin(), iterEnd = m_settings.m_v0VertexCollections.end(); iter != iterEnd; ++iter)
    {
        if(collectionMaps.collectionMap_Vertex.find(*iter) == collectionMaps.collectionMap_Vertex.end()) { std::cout<<"not find "<<(*iter)<<std::endl; continue;}
        try
        {
            auto pV0Collection = (collectionMaps.collectionMap_Vertex.find(*iter))->second;

            for (int i = 0, iMax = pV0Collection.size(); i < iMax; ++i)
            {
                try
                {
                    auto  pVertex0 = pV0Collection.at(i);
                    auto  pVertex  = &(pVertex0);

                    if (NULL == pVertex) throw ("Collection type mismatch");

                    auto pReconstructedParticle = pVertex->getAssociatedParticle();

                    if (this->IsConflictingRelationship(pReconstructedParticle))continue;

                    // Extract the v0 vertex information
                    const int vertexPdgCode(pReconstructedParticle.getType());

                    for (unsigned int iTrack = 0, nTracks = pReconstructedParticle.tracks_size(); iTrack < nTracks; ++iTrack)
                    {
                        auto pTrack = pReconstructedParticle.getTracks(iTrack);
                        m_v0TrackList.insert(pTrack.id().index);

                        int trackPdgCode = pandora::UNKNOWN_PARTICLE_TYPE;

                        switch (vertexPdgCode)
                        {
                        case pandora::PHOTON :
                            (pTrack.getTrackStates(0).omega > 0) ? trackPdgCode = pandora::E_PLUS : trackPdgCode = pandora::E_MINUS;
                            break;
                        case pandora::LAMBDA :
                            (pTrack.getTrackStates(0).omega > 0) ? trackPdgCode = pandora::PROTON : trackPdgCode = pandora::PI_MINUS;
                            break;
                        case pandora::LAMBDA_BAR :
                            (pTrack.getTrackStates(0).omega > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PROTON_BAR;
                            break;
                        case pandora::K_SHORT :
                            (pTrack.getTrackStates(0).omega > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PI_MINUS;
                            break;
                        default :
                            (pTrack.getTrackStates(0).omega > 0) ? trackPdgCode = pandora::PI_PLUS : trackPdgCode = pandora::PI_MINUS;
                            break;
                        }

                        m_trackToPidMap.insert(TrackToPidMap::value_type(pTrack, trackPdgCode));

                        if (0 == m_settings.m_shouldFormTrackRelationships) continue;

                        // Make track sibling relationships
                        for (unsigned int jTrack = iTrack + 1; jTrack < nTracks; ++jTrack)
                        {
                            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackSiblingRelationship(*m_pPandora, GetTrackAddress(collectionMaps, pTrack), GetTrackAddress(collectionMaps, pReconstructedParticle.getTracks(jTrack) ) ) );
                        }
                    }
                }
                catch (...)
                {
                    std::cout<< "Failed to extract v0 vertex" << std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout << "Failed to extract v0 vertex collection: " << *iter << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackCreator::IsConflictingRelationship(const edm4hep::ReconstructedParticle &Particle) const
{
    for (unsigned int iTrack = 0, nTracks = Particle.tracks_size(); iTrack < nTracks; ++iTrack)
    {
        edm4hep::Track pTrack = Particle.getTracks(iTrack) ;
        unsigned int pTrack_id = pTrack.id().index ;

        if (this->IsDaughter(pTrack_id) || this->IsParent(pTrack_id) || this->IsV0(pTrack_id))
            return true;
    }

    return false;
}

edm4hep::Track* TrackCreator::GetTrackAddress(const CollectionMaps& collectionMaps, const edm4hep::Track& pTrack )
{
    for (StringVector::const_iterator iter = m_settings.m_trackCollections.begin(), iterEnd = m_settings.m_trackCollections.end(); iter != iterEnd; ++iter)
    {
        if(collectionMaps.collectionMap_Track.find(*iter) == collectionMaps.collectionMap_Track.end()) { std::cout<<"not find "<<(*iter)<<std::endl; continue;}
        auto pTrackCollection = (collectionMaps.collectionMap_Track.find(*iter))->second;
        for (int i = 0, iMax = pTrackCollection.size(); i < iMax; ++i)
        {
            auto pTrack0 = pTrackCollection.at(i);
            if (pTrack.id() == pTrack0.id()) return (&pTrack0);
        }
    }
    return NULL;
}
//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode TrackCreator::CreateTracks(const CollectionMaps& collectionMaps)
{
    std::cout<<"start TrackCreator::CreateTracks:"<<std::endl;
    for (StringVector::const_iterator iter = m_settings.m_trackCollections.begin(), iterEnd = m_settings.m_trackCollections.end();
        iter != iterEnd; ++iter)
    {
        if(collectionMaps.collectionMap_Track.find(*iter) == collectionMaps.collectionMap_Track.end()) { std::cout<<"not find "<<(*iter)<<std::endl; continue;}
        try
        {
            auto pTrackCollection = (collectionMaps.collectionMap_Track.find(*iter))->second;
            std::cout<<(*iter)<<" size="<<pTrackCollection.size()<<std::endl;

            for (int i = 0, iMax = pTrackCollection.size(); i < iMax; ++i)
            {
                try
                {
                    auto pTrack0 = pTrackCollection.at(i);
                    auto pTrack  = (&pTrack0);

                    if (NULL == pTrack) throw ("Collection type mismatch");

                    int minTrackHits = m_settings.m_minTrackHits;
                    const float tanLambda(std::fabs(pTrack->getTrackStates(0).tanLambda));

                    if (tanLambda > m_tanLambdaFtd)
                    {
                        int expectedFtdHits(0);

                        for (unsigned int iFtdLayer = 0; iFtdLayer < m_nFtdLayers; ++iFtdLayer)
                        {
                            if ((tanLambda > m_ftdZPositions[iFtdLayer] / m_ftdOuterRadii[iFtdLayer]) &&
                                (tanLambda < m_ftdZPositions[iFtdLayer] / m_ftdInnerRadii[iFtdLayer]))
                            {
                                expectedFtdHits++;
                            }
                        }

                        minTrackHits = std::max(m_settings.m_minFtdTrackHits, expectedFtdHits);
                    }

                    const int nTrackHits(static_cast<int>(pTrack->trackerHits_size()));

                    if ((nTrackHits < minTrackHits) || (nTrackHits > m_settings.m_maxTrackHits))
                        continue;

                    // Proceed to create the pandora track
                    PandoraApi::Track::Parameters trackParameters;
                    trackParameters.m_d0 = pTrack->getTrackStates(0).D0;
                    trackParameters.m_z0 = pTrack->getTrackStates(0).Z0;
                    trackParameters.m_pParentAddress = pTrack;
                    // By default, assume tracks are charged pions
                    const float signedCurvature(pTrack->getTrackStates(0).omega);
                    trackParameters.m_particleId = (signedCurvature > 0) ? pandora::PI_PLUS : pandora::PI_MINUS;
                    trackParameters.m_mass = pandora::PdgTable::GetParticleMass(pandora::PI_PLUS);

                    // Use particle id information from V0 and Kink finders
                    TrackToPidMap::const_iterator iter_t = m_trackToPidMap.find(*pTrack);

                    if(iter_t != m_trackToPidMap.end())
                    {
                        trackParameters.m_particleId = (*iter_t).second;
                        trackParameters.m_mass = pandora::PdgTable::GetParticleMass((*iter_t).second);
                    }

                    if (std::numeric_limits<float>::epsilon() < std::fabs(signedCurvature))
                        trackParameters.m_charge = static_cast<int>(signedCurvature / std::fabs(signedCurvature));

                    this->GetTrackStates(pTrack, trackParameters);
                    this->TrackReachesECAL(pTrack, trackParameters);
                    this->DefineTrackPfoUsage(pTrack, trackParameters);

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Track::Create(*m_pPandora, trackParameters));
                    m_trackVector.push_back(pTrack);
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    std::cout<<"ERROR Failed to extract a track: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    std::cout << "WARNNING Failed to extract a track "<< std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout<<"WARNING Failed to extract track collection: " << *iter << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackCreator::GetTrackStates(edm4hep::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    edm4hep::TrackState pTrackState = pTrack->getTrackStates(1); // ref  /cvmfs/cepcsw.ihep.ac.cn/prototype/LCIO/include/EVENT/TrackState.h 


    const double pt(m_bField * 2.99792e-4 / std::fabs(pTrackState.omega));
    trackParameters.m_momentumAtDca = pandora::CartesianVector(std::cos(pTrackState.phi), std::sin(pTrackState.phi), pTrackState.tanLambda) * pt;
    this->CopyTrackState(pTrack->getTrackStates(2), trackParameters.m_trackStateAtStart);

    auto pEndTrack = (pTrack->tracks_size() ==0 ) ?  *pTrack  :  pTrack->getTracks(pTrack->tracks_size()-1);

    
    
    this->CopyTrackState(pEndTrack.getTrackStates(3), trackParameters.m_trackStateAtEnd);
    //FIXME ? LCIO input only has 4 states, so 4 can't be used.
    if( pEndTrack.trackStates_size()<5) this->CopyTrackState(pEndTrack.getTrackStates(3), trackParameters.m_trackStateAtCalorimeter);
    else  this->CopyTrackState(pEndTrack.getTrackStates(4), trackParameters.m_trackStateAtCalorimeter);
    
    
    trackParameters.m_isProjectedToEndCap = ((std::fabs(trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetZ()) < m_eCalEndCapInnerZ) ? false : true);

    // Convert generic time (length from reference point to intersection, divided by momentum) into nanoseconds
    const float minGenericTime(this->CalculateTrackTimeAtCalorimeter(pTrack));
    const float particleMass(trackParameters.m_mass.Get());
    const float particleEnergy(std::sqrt(particleMass * particleMass + trackParameters.m_momentumAtDca.Get().GetMagnitudeSquared()));
    trackParameters.m_timeAtCalorimeter = minGenericTime * particleEnergy / 299.792f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TrackCreator::CalculateTrackTimeAtCalorimeter(edm4hep::Track *const pTrack) const
{
    const pandora::Helix helix(pTrack->getTrackStates(0).phi, pTrack->getTrackStates(0).D0, pTrack->getTrackStates(0).Z0, pTrack->getTrackStates(0).omega, pTrack->getTrackStates(0).tanLambda, m_bField);
    const pandora::CartesianVector &referencePoint(helix.GetReferencePoint());

    // First project to endcap
    float minGenericTime(std::numeric_limits<float>::max());

    pandora::CartesianVector bestECalProjection(0.f, 0.f, 0.f);
    const int signPz((helix.GetMomentum().GetZ() > 0.f) ? 1 : -1);
    (void) helix.GetPointInZ(static_cast<float>(signPz) * m_eCalEndCapInnerZ, referencePoint, bestECalProjection, minGenericTime);

    // Then project to barrel surface(s)
    pandora::CartesianVector barrelProjection(0.f, 0.f, 0.f);
    if (m_eCalBarrelInnerSymmetry > 0)
    {
        // Polygon
        float twopi_n = 2. * M_PI / (static_cast<float>(m_eCalBarrelInnerSymmetry));

        for (int i = 0; i < m_eCalBarrelInnerSymmetry; ++i)
        {
            float genericTime(std::numeric_limits<float>::max());
            const float phi(twopi_n * static_cast<float>(i) + m_eCalBarrelInnerPhi0);

            const pandora::StatusCode statusCode(helix.GetPointInXY(m_eCalBarrelInnerR * std::cos(phi), m_eCalBarrelInnerR * std::sin(phi),
                std::cos(phi + 0.5 * M_PI), std::sin(phi + 0.5 * M_PI), referencePoint, barrelProjection, genericTime));

            if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime))
            {
                minGenericTime = genericTime;
                bestECalProjection = barrelProjection;
            }
        }
    }
    else
    {
        // Cylinder
        float genericTime(std::numeric_limits<float>::max());
        const pandora::StatusCode statusCode(helix.GetPointOnCircle(m_eCalBarrelInnerR, referencePoint, barrelProjection, genericTime));

        if ((pandora::STATUS_CODE_SUCCESS == statusCode) && (genericTime < minGenericTime))
        {
            minGenericTime = genericTime;
            bestECalProjection = barrelProjection;
        }
    }

    if (bestECalProjection.GetMagnitudeSquared() < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_NOT_INITIALIZED);

    return minGenericTime;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackCreator::CopyTrackState(const edm4hep::TrackState & pTrackState, pandora::InputTrackState &inputTrackState) const
{
    const double pt(m_bField * 2.99792e-4 / std::fabs(pTrackState.omega));

    const double px(pt * std::cos(pTrackState.phi));
    const double py(pt * std::sin(pTrackState.phi));
    const double pz(pt * pTrackState.tanLambda);

    const double xs(pTrackState.referencePoint[0]);
    const double ys(pTrackState.referencePoint[1]);
    const double zs(pTrackState.referencePoint[2]);

    inputTrackState = pandora::TrackState(xs, ys, zs, px, py, pz);
}


//------------------------------------------------------------------------------------------------------------------------------------------

void TrackCreator::TrackReachesECAL(edm4hep::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    
    // Calculate hit position information
    float hitZMin(std::numeric_limits<float>::max());
    float hitZMax(-std::numeric_limits<float>::max());
    float hitOuterR(-std::numeric_limits<float>::max());

    int maxOccupiedFtdLayer=0;
    
    
    const unsigned int nTrackHits(pTrack->trackerHits_size());
    for (unsigned int i = 0; i < nTrackHits; ++i)
    {
        
        const edm4hep::TrackerHit Hit ( pTrack->getTrackerHits(i) );
        const edm4hep::Vector3d pos = Hit.getPosition();
        
        float x = float(pos[0]);
        float y = float(pos[1]);
        float z = float(pos[2]);
        float r = std::sqrt(x * x + y * y);

        if (z > hitZMax) hitZMax = z;

        if (z < hitZMin) hitZMin = z;

        if (r > hitOuterR) hitOuterR = r;

        if ((r > m_tpcInnerR) && (r < m_tpcOuterR) && (std::fabs(z) <= m_tpcZmax))  continue;

        for (unsigned int j = 0; j < m_nFtdLayers; ++j)
        {
            if ((r > m_ftdInnerRadii[j]) && (r < m_ftdOuterRadii[j]) &&
                (std::fabs(z) - m_settings.m_reachesECalFtdZMaxDistance < m_ftdZPositions[j]) &&
                (std::fabs(z) + m_settings.m_reachesECalFtdZMaxDistance > m_ftdZPositions[j]))
            {
                if ( j > maxOccupiedFtdLayer) maxOccupiedFtdLayer = j;
                break;
            }
        }
    }
    const int nTpcHits(this->GetNTpcHits(pTrack));
    const int nFtdHits(this->GetNFtdHits(pTrack));

    // Look to see if there are hits in etd or set, implying track has reached edge of ecal
    if ((hitOuterR > m_minSetRadius) || (hitZMax > m_minEtdZPosition))
    {
        trackParameters.m_reachesCalorimeter = true;
        return;
    }

    // Require sufficient hits in tpc or ftd, then compare extremal hit positions with tracker dimensions
    if ((nTpcHits >= m_settings.m_reachesECalNTpcHits) || (nFtdHits >= m_settings.m_reachesECalNFtdHits))
    {
        if ((hitOuterR - m_tpcOuterR > m_settings.m_reachesECalTpcOuterDistance) ||
            (std::fabs(hitZMax) - m_tpcZmax > m_settings.m_reachesECalTpcZMaxDistance) ||
            (std::fabs(hitZMin) - m_tpcZmax > m_settings.m_reachesECalTpcZMaxDistance) ||
            (maxOccupiedFtdLayer >= m_settings.m_reachesECalMinFtdLayer))
        {
            trackParameters.m_reachesCalorimeter = true;
            return;
        }
    }

    // If track is lowpt, it may curl up and end inside tpc inner radius
    const pandora::CartesianVector &momentumAtDca(trackParameters.m_momentumAtDca.Get());
    const float cosAngleAtDca(std::fabs(momentumAtDca.GetZ()) / momentumAtDca.GetMagnitude());
    const float pX(momentumAtDca.GetX()), pY(momentumAtDca.GetY());
    const float pT(std::sqrt(pX * pX + pY * pY));

    if ((cosAngleAtDca > m_cosTpc) || (pT < m_settings.m_curvatureToMomentumFactor * m_bField * m_tpcOuterR))
    {
        trackParameters.m_reachesCalorimeter = true;
        return;
    }
    
    trackParameters.m_reachesCalorimeter = false;
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackCreator::DefineTrackPfoUsage(edm4hep::Track *const pTrack, PandoraApi::Track::Parameters &trackParameters) const
{
    bool canFormPfo(false);
    bool canFormClusterlessPfo(false);

    if (trackParameters.m_reachesCalorimeter.Get() && !this->IsParent(pTrack->id().index))
    {
        const float d0(std::fabs(pTrack->getTrackStates(0).D0)), z0(std::fabs(pTrack->getTrackStates(0).Z0));

        float rInner(std::numeric_limits<float>::max()), zMin(std::numeric_limits<float>::max());

        for (auto iter = pTrack->trackerHits_begin(), iterEnd = pTrack->trackerHits_end(); iter != iterEnd; ++iter)
        {
            const edm4hep::Vector3d pPosition = (*iter).getPosition(); 
            const float x(pPosition[0]), y(pPosition[1]), absoluteZ(std::fabs(pPosition[2]));
            const float r(std::sqrt(x * x + y * y));

            if (r < rInner)
                rInner = r;

            if (absoluteZ < zMin)
                zMin = absoluteZ;
        }

        if (this->PassesQualityCuts(pTrack, trackParameters))
        {
            const pandora::CartesianVector &momentumAtDca(trackParameters.m_momentumAtDca.Get());
            const float pX(momentumAtDca.GetX()), pY(momentumAtDca.GetY()), pZ(momentumAtDca.GetZ());
            const float pT(std::sqrt(pX * pX + pY * pY));

            const float zCutForNonVertexTracks(m_tpcInnerR * std::fabs(pZ / pT) + m_settings.m_zCutForNonVertexTracks);
            const bool passRzQualityCuts((zMin < zCutForNonVertexTracks) && (rInner < m_tpcInnerR + m_settings.m_maxTpcInnerRDistance));

            const bool isV0(this->IsV0(pTrack->id().index));
            const bool isDaughter(this->IsDaughter(pTrack->id().index));

            // Decide whether track can be associated with a pandora cluster and used to form a charged PFO
            if ((d0 < m_settings.m_d0TrackCut) && (z0 < m_settings.m_z0TrackCut) && (rInner < m_tpcInnerR + m_settings.m_maxTpcInnerRDistance))
            {
                canFormPfo = true;
            }
            else if (passRzQualityCuts && (0 != m_settings.m_usingNonVertexTracks))
            {
                canFormPfo = true;
            }
            else if (isV0 || isDaughter)
            {
                canFormPfo = true;
            }

            // Decide whether track can be used to form a charged PFO, even if track fails to be associated with a pandora cluster
            const float particleMass(trackParameters.m_mass.Get());
            const float trackEnergy(std::sqrt(momentumAtDca.GetMagnitudeSquared() + particleMass * particleMass));

            if ((0 != m_settings.m_usingUnmatchedVertexTracks) && (trackEnergy < m_settings.m_unmatchedVertexTrackMaxEnergy))
            {
                if ((d0 < m_settings.m_d0UnmatchedVertexTrackCut) && (z0 < m_settings.m_z0UnmatchedVertexTrackCut) &&
                    (rInner < m_tpcInnerR + m_settings.m_maxTpcInnerRDistance))
                {
                    canFormClusterlessPfo = true;
                }
                else if (passRzQualityCuts && (0 != m_settings.m_usingNonVertexTracks) && (0 != m_settings.m_usingUnmatchedNonVertexTracks))
                {
                    canFormClusterlessPfo = true;
                }
                else if (isV0 || isDaughter)
                {
                    canFormClusterlessPfo = true;
                }
            }
        }
        else if (this->IsDaughter(pTrack->id().index) || this->IsV0(pTrack->id().index))
        {
            std::cout<<"WARNING Recovering daughter or v0 track " << trackParameters.m_momentumAtDca.Get().GetMagnitude() << std::endl;
            canFormPfo = true;
        }
    }

    trackParameters.m_canFormPfo = canFormPfo;
    trackParameters.m_canFormClusterlessPfo = canFormClusterlessPfo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackCreator::PassesQualityCuts(edm4hep::Track *const pTrack, const PandoraApi::Track::Parameters &trackParameters) const
{
    // First simple sanity checks
    if (trackParameters.m_trackStateAtCalorimeter.Get().GetPosition().GetMagnitude() < m_settings.m_minTrackECalDistanceFromIp)
        return false;

    if (std::fabs(pTrack->getTrackStates(0).omega) < std::numeric_limits<float>::epsilon())
    {
        std::cout<<"ERROR Track has Omega = 0 " << std::endl;
        return false;
    }

    // Check momentum uncertainty is reasonable to use track
    const pandora::CartesianVector &momentumAtDca(trackParameters.m_momentumAtDca.Get());
    const float sigmaPOverP(std::sqrt(pTrack->getTrackStates(0).covMatrix[5]) / std::fabs(pTrack->getTrackStates(0).omega));

    if (sigmaPOverP > m_settings.m_maxTrackSigmaPOverP)
    {
        std::cout<<"WARNING Dropping track : " << momentumAtDca.GetMagnitude() << "+-" << sigmaPOverP * (momentumAtDca.GetMagnitude())
                               << " chi2 = " <<  pTrack->getChi2() << " " << pTrack->getNdf()
                               << " from " << pTrack->trackerHits_size() << std::endl;
        return false;
    }

    // Require reasonable number of TPC hits 
    if (momentumAtDca.GetMagnitude() > m_settings.m_minMomentumForTrackHitChecks)
    {
        const float pX(fabs(momentumAtDca.GetX()));
        const float pY(fabs(momentumAtDca.GetY()));
        const float pZ(fabs(momentumAtDca.GetZ()));
        const float pT(std::sqrt(pX * pX + pY * pY));
        const float rInnermostHit(pTrack->getRadiusOfInnermostHit());

        if ((std::numeric_limits<float>::epsilon() > std::fabs(pT)) || (std::numeric_limits<float>::epsilon() > std::fabs(pZ)) || (rInnermostHit == m_tpcOuterR))
        {
            std::cout<<"ERROR Invalid track parameter, pT " << pT << ", pZ " << pZ << ", rInnermostHit " << rInnermostHit << std::endl;
            return false;
        }

        float nExpectedTpcHits(0.);

        if (pZ < m_tpcZmax / m_tpcOuterR * pT)
        {
            const float innerExpectedHitRadius(std::max(m_tpcInnerR, rInnermostHit));
            const float frac((m_tpcOuterR - innerExpectedHitRadius) / (m_tpcOuterR - m_tpcInnerR));
            nExpectedTpcHits = m_tpcMaxRow * frac;
        }

        if ((pZ <= m_tpcZmax / m_tpcInnerR * pT) && (pZ >= m_tpcZmax / m_tpcOuterR * pT))
        {
            const float innerExpectedHitRadius(std::max(m_tpcInnerR, rInnermostHit));
            const float frac((m_tpcZmax * pT / pZ - innerExpectedHitRadius) / (m_tpcOuterR - innerExpectedHitRadius));
            nExpectedTpcHits = frac * m_tpcMaxRow;
        }

        // TODO Get TPC membrane information from GEAR when available
        if (std::fabs(pZ) / momentumAtDca.GetMagnitude() < m_settings.m_tpcMembraneMaxZ / m_tpcInnerR)
            nExpectedTpcHits = 0;

        const int nTpcHits(this->GetNTpcHits(pTrack));
        const int nFtdHits(this->GetNFtdHits(pTrack));

        const int minTpcHits = static_cast<int>(nExpectedTpcHits * m_settings.m_minTpcHitFractionOfExpected);

        if ((nTpcHits < minTpcHits) && (nFtdHits < m_settings.m_minFtdHitsForTpcHitFraction))
        {
            std::cout<<"WARNING Dropping track : " << momentumAtDca.GetMagnitude() << " Number of TPC hits = " << nTpcHits
                                   << " < " << minTpcHits << " nftd = " << nFtdHits  << std::endl;
            return false;
        }
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int TrackCreator::GetNTpcHits(edm4hep::Track *const pTrack) const
{
#if EDM4HEP_BUILD_VERSION > EDM4HEP_VERSION(0, 9, 0)
    return pTrack->getSubdetectorHitNumbers(2 * lcio::ILDDetID::TPC - 1);// still use LCIO code now
#else
    return pTrack->getSubDetectorHitNumbers(2 * lcio::ILDDetID::TPC - 1);// still use LCIO code now
#endif
}

//------------------------------------------------------------------------------------------------------------------------------------------

int TrackCreator::GetNFtdHits(edm4hep::Track *const pTrack) const
{
#if EDM4HEP_BUILD_VERSION > EDM4HEP_VERSION(0, 9, 0)
    return pTrack->getSubdetectorHitNumbers( 2 * lcio::ILDDetID::FTD - 1 );
#else
    return pTrack->getSubDetectorHitNumbers( 2 * lcio::ILDDetID::FTD - 1 );
#endif
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TrackCreator::Settings::Settings() :
    m_shouldFormTrackRelationships(1),
    m_minTrackHits(5),
    m_minFtdTrackHits(0),
    m_maxTrackHits(5000.f),
    m_d0TrackCut(50.f),
    m_z0TrackCut(50.f),
    m_usingNonVertexTracks(1),
    m_usingUnmatchedNonVertexTracks(0),
    m_usingUnmatchedVertexTracks(1),
    m_unmatchedVertexTrackMaxEnergy(5.f),
    m_d0UnmatchedVertexTrackCut(5.f),
    m_z0UnmatchedVertexTrackCut(5.f),
    m_zCutForNonVertexTracks(250.f),
    m_reachesECalNTpcHits(11),
    m_reachesECalNFtdHits(4),
    m_reachesECalTpcOuterDistance(-100.f),
    m_reachesECalMinFtdLayer(9),
    m_reachesECalTpcZMaxDistance(-50.f),
    m_reachesECalFtdZMaxDistance(1.f),
    m_curvatureToMomentumFactor(0.3f / 2000.f),
    m_minTrackECalDistanceFromIp(100.f),
    m_maxTrackSigmaPOverP(0.15f),
    m_minMomentumForTrackHitChecks(1.f),
    m_tpcMembraneMaxZ(10.f),
    m_maxTpcInnerRDistance(50.f),
    m_minTpcHitFractionOfExpected(0.2f),
    m_minFtdHitsForTpcHitFraction(2)
{
}
