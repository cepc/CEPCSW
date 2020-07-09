/**
 *  @brief  Implementation of the pfo creator class.
 * 
 *  $Log: $
 */


#include "Api/PandoraApi.h"

#include "Objects/Cluster.h"
#include "Objects/ParticleFlowObject.h"
#include "Objects/Track.h"

#include "Pandora/PdgTable.h"
#include "PfoCreator.h"
#include "PandoraPFAlg.h"

#include <algorithm>
#include <cmath>

PfoCreator::PfoCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pPandora(pPandora)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoCreator::~PfoCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode PfoCreator::CreateParticleFlowObjects(CollectionMaps& collectionMaps, DataHandle<edm4hep::ClusterCollection>& _pClusterCollection, DataHandle<edm4hep::ReconstructedParticleCollection>& _pReconstructedParticleCollection, DataHandle<edm4hep::VertexCollection>& _pStartVertexCollection)
{
    m_collectionMaps = &collectionMaps;
    edm4hep::ClusterCollection* pClusterCollection                              = _pClusterCollection.createAndPut();
    edm4hep::ReconstructedParticleCollection* pReconstructedParticleCollection  = _pReconstructedParticleCollection.createAndPut();
    edm4hep::VertexCollection* pStartVertexCollection                           = _pStartVertexCollection.createAndPut();
 
    const pandora::PfoList *pPandoraPfoList = NULL;
    PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::GetCurrentPfoList(*m_pPandora, pPandoraPfoList));

    pandora::StringVector subDetectorNames;
    this->InitialiseSubDetectorNames(subDetectorNames);

    std::cout<<"pPandoraPfoList size="<<pPandoraPfoList->size()<<std::endl;
    for (pandora::PfoList::const_iterator pIter = pPandoraPfoList->begin(), pIterEnd = pPandoraPfoList->end(); pIter != pIterEnd; ++pIter)
    {
        const pandora::ParticleFlowObject *const pPandoraPfo(*pIter);
        edm4hep::ReconstructedParticle pReconstructedParticle0 = pReconstructedParticleCollection->create();
        edm4hep::ReconstructedParticle* pReconstructedParticle = &pReconstructedParticle0;

        const bool hasTrack(!pPandoraPfo->GetTrackList().empty());
        const pandora::ClusterList &clusterList(pPandoraPfo->GetClusterList());

        float clustersTotalEnergy(0.f);
        pandora::CartesianVector referencePoint(0.f, 0.f, 0.f), clustersWeightedPosition(0.f, 0.f, 0.f);
        for (pandora::ClusterList::const_iterator cIter = clusterList.begin(), cIterEnd = clusterList.end(); cIter != cIterEnd; ++cIter)
        {
            const pandora::Cluster *const pPandoraCluster(*cIter);
            pandora::CaloHitList pandoraCaloHitList;
            pPandoraCluster->GetOrderedCaloHitList().FillCaloHitList(pandoraCaloHitList);
            pandoraCaloHitList.insert(pandoraCaloHitList.end(), pPandoraCluster->GetIsolatedCaloHitList().begin(), pPandoraCluster->GetIsolatedCaloHitList().end());

            pandora::FloatVector hitE, hitX, hitY, hitZ;
            edm4hep::Cluster p_Cluster0 = pClusterCollection->create();
            edm4hep::Cluster* p_Cluster = &p_Cluster0;
            this->SetClusterSubDetectorEnergies(subDetectorNames, p_Cluster, pandoraCaloHitList, hitE, hitX, hitY, hitZ);

            float clusterCorrectEnergy(0.f);
            this->SetClusterEnergyAndError(pPandoraPfo, pPandoraCluster, p_Cluster, clusterCorrectEnergy);

            pandora::CartesianVector clusterPosition(0.f, 0.f, 0.f);
            const unsigned int nHitsInCluster(pandoraCaloHitList.size());
            this->SetClusterPositionAndError(nHitsInCluster, hitE, hitX, hitY, hitZ, p_Cluster, clusterPosition);

            if (!hasTrack)
            {
                clustersWeightedPosition += clusterPosition * clusterCorrectEnergy;
                clustersTotalEnergy += clusterCorrectEnergy;
            }

            edm4hep::ConstCluster p_ClusterCon = *p_Cluster;
            pReconstructedParticle->addToClusters(p_ClusterCon);
        }

        if (!hasTrack)
        {
            if (clustersTotalEnergy < std::numeric_limits<float>::epsilon())
            {
                std::cout<<"WARNING PfoCreator::CreateParticleFlowObjects: invalid cluster energy " << clustersTotalEnergy << std::endl;
                throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
            }
            else
            {
                referencePoint = clustersWeightedPosition * (1.f / clustersTotalEnergy);
            }
        }
        else
        {
            PANDORA_RETURN_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, this->CalculateTrackBasedReferencePoint(pPandoraPfo, referencePoint));
        }

        this->SetRecoParticleReferencePoint(referencePoint, pReconstructedParticle);
        this->AddTracksToRecoParticle(pPandoraPfo, pReconstructedParticle);
        this->SetRecoParticlePropertiesFromPFO(pPandoraPfo, pReconstructedParticle);

        edm4hep::Vertex pStartVertex0 = pStartVertexCollection->create();
        edm4hep::Vertex* pStartVertex = &pStartVertex0;
        pStartVertex->setAlgorithmType(0);
        const float ref_value[3] = {referencePoint.GetX(),referencePoint.GetY(),referencePoint.GetZ()};
        pStartVertex->setPosition(edm4hep::Vector3f(ref_value));
        pStartVertex->setAssociatedParticle(*pReconstructedParticle);
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoCreator::InitialiseSubDetectorNames(pandora::StringVector &subDetectorNames) const
{
    subDetectorNames.push_back("ecal");
    subDetectorNames.push_back("hcal");
    subDetectorNames.push_back("yoke");
    subDetectorNames.push_back("lcal");
    subDetectorNames.push_back("lhcal");
    subDetectorNames.push_back("bcal");
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoCreator::SetClusterSubDetectorEnergies(const pandora::StringVector &subDetectorNames, edm4hep::Cluster *const p_Cluster,
    const pandora::CaloHitList &pandoraCaloHitList, pandora::FloatVector &hitE, pandora::FloatVector &hitX, pandora::FloatVector &hitY,
    pandora::FloatVector &hitZ) const
{
    for (pandora::CaloHitList::const_iterator hIter = pandoraCaloHitList.begin(), hIterEnd = pandoraCaloHitList.end(); hIter != hIterEnd; ++hIter)
    {
        const pandora::CaloHit *const pPandoraCaloHit(*hIter);
        edm4hep::CalorimeterHit *const pCalorimeterHit0 = (edm4hep::CalorimeterHit*)(pPandoraCaloHit->GetParentAddress());
        const edm4hep::CalorimeterHit pCalorimeterHit = *pCalorimeterHit0;
        
        p_Cluster->addToHits(pCalorimeterHit);

        const float caloHitEnergy(pCalorimeterHit.getEnergy());
        hitE.push_back(caloHitEnergy);
        hitX.push_back(pCalorimeterHit.getPosition()[0]);
        hitY.push_back(pCalorimeterHit.getPosition()[1]);
        hitZ.push_back(pCalorimeterHit.getPosition()[2]);
        /*Can be added later 
        std::vector<float> &subDetectorEnergies = p_Cluster->subdetectorEnergies();
        subDetectorEnergies.resize(subDetectorNames.size());

        switch (CHT(pCalorimeterHit->getType()).caloID())
        {
            case CHT::ecal:  subDetectorEnergies[ECAL_INDEX ] += caloHitEnergy; break;
            case CHT::hcal:  subDetectorEnergies[HCAL_INDEX ] += caloHitEnergy; break;
            case CHT::yoke:  subDetectorEnergies[YOKE_INDEX ] += caloHitEnergy; break;
            case CHT::lcal:  subDetectorEnergies[LCAL_INDEX ] += caloHitEnergy; break;
            case CHT::lhcal: subDetectorEnergies[LHCAL_INDEX] += caloHitEnergy; break;
            case CHT::bcal:  subDetectorEnergies[BCAL_INDEX ] += caloHitEnergy; break;
            default: streamlog_out(WARNING) << "PfoCreator::SetClusterSubDetectorEnergies: no subdetector found for hit with type: " << pCalorimeterHit->getType() << std::endl;
        }
        */
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoCreator::SetClusterEnergyAndError(const pandora::ParticleFlowObject *const pPandoraPfo, const pandora::Cluster *const pPandoraCluster, 
    edm4hep::Cluster *const p_Cluster, float &clusterCorrectEnergy) const
{
    const bool isEmShower((pandora::PHOTON == pPandoraPfo->GetParticleId()) || (pandora::E_MINUS == std::abs(pPandoraPfo->GetParticleId())));
    clusterCorrectEnergy = (isEmShower ? pPandoraCluster->GetCorrectedElectromagneticEnergy(*m_pPandora) : pPandoraCluster->GetCorrectedHadronicEnergy(*m_pPandora));

    if (clusterCorrectEnergy < std::numeric_limits<float>::epsilon())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);

    const float stochasticTerm(isEmShower ? m_settings.m_emStochasticTerm : m_settings.m_hadStochasticTerm); 
    const float constantTerm(isEmShower ? m_settings.m_emConstantTerm : m_settings.m_hadConstantTerm);
    const float energyError(std::sqrt(stochasticTerm * stochasticTerm / clusterCorrectEnergy + constantTerm * constantTerm) * clusterCorrectEnergy);

    p_Cluster->setEnergy(clusterCorrectEnergy);
    p_Cluster->setEnergyError(energyError);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoCreator::SetClusterPositionAndError(const unsigned int nHitsInCluster, pandora::FloatVector &hitE, pandora::FloatVector &hitX, 
    pandora::FloatVector &hitY, pandora::FloatVector &hitZ, edm4hep::Cluster *const p_Cluster, pandora::CartesianVector &clusterPositionVec) const
{
    ClusterShapes *const pClusterShapes(new ClusterShapes(nHitsInCluster, hitE.data(), hitX.data(), hitY.data(), hitZ.data()));//this need GSL/1.14 

    try
    {
        p_Cluster->setPhi(std::atan2(pClusterShapes->getEigenVecInertia()[1], pClusterShapes->getEigenVecInertia()[0]));
        p_Cluster->setITheta(std::acos(pClusterShapes->getEigenVecInertia()[2]));
        p_Cluster->setPosition(pClusterShapes->getCentreOfGravity());
        clusterPositionVec.SetValues(pClusterShapes->getCentreOfGravity()[0], pClusterShapes->getCentreOfGravity()[1], pClusterShapes->getCentreOfGravity()[2]);
    }
    catch (...)
    {
        std::cout<<"WARNING PfoCreator::SetClusterPositionAndError: unidentified exception caught." << std::endl;
    }

    delete pClusterShapes;
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode PfoCreator::CalculateTrackBasedReferencePoint(const pandora::ParticleFlowObject *const pPandoraPfo, pandora::CartesianVector &referencePoint) const
{
    const pandora::TrackList &trackList(pPandoraPfo->GetTrackList());

    float totalTrackMomentumAtDca(0.f), totalTrackMomentumAtStart(0.f);
    pandora::CartesianVector referencePointAtDCAWeighted(0.f, 0.f, 0.f), referencePointAtStartWeighted(0.f, 0.f, 0.f);

    bool hasSiblings(false);
    for (pandora::TrackList::const_iterator tIter = trackList.begin(), tIterEnd = trackList.end(); tIter != tIterEnd; ++tIter)
    {
        const pandora::Track *const pPandoraTrack(*tIter);

        if (!this->IsValidParentTrack(pPandoraTrack, trackList))
            continue;

        if (this->HasValidSiblingTrack(pPandoraTrack, trackList))
        {
            // Presence of sibling tracks typically represents a conversion
            const pandora::CartesianVector &trackStartPoint((pPandoraTrack->GetTrackStateAtStart()).GetPosition());
            const float trackStartMomentum(((pPandoraTrack->GetTrackStateAtStart()).GetMomentum()).GetMagnitude());
            referencePointAtStartWeighted += trackStartPoint * trackStartMomentum;
            totalTrackMomentumAtStart += trackStartMomentum;
            hasSiblings = true;
        }
        else
        {
            const edm4hep::Track *const pLcioTrack0 = (edm4hep::Track*)(pPandoraTrack->GetParentAddress());
            const edm4hep::Track pLcioTrack = *pLcioTrack0;

            const float z0(pPandoraTrack->GetZ0());
            pandora::CartesianVector intersectionPoint(0.f, 0.f, 0.f);

            if(pLcioTrack.trackStates_size()==0) throw "zero trackStates size find";
            intersectionPoint.SetValues(pLcioTrack.getTrackStates(0).D0 * std::cos(pLcioTrack.getTrackStates(0).phi), pLcioTrack.getTrackStates(0).D0 * std::sin(pLcioTrack.getTrackStates(0).phi), z0);
            const float trackMomentumAtDca((pPandoraTrack->GetMomentumAtDca()).GetMagnitude());
            referencePointAtDCAWeighted += intersectionPoint * trackMomentumAtDca;
            totalTrackMomentumAtDca += trackMomentumAtDca;
        }
    }

    if (hasSiblings)
    {
        if (totalTrackMomentumAtStart < std::numeric_limits<float>::epsilon())
        {
            std::cout<<" WARNING PfoCreator::CalculateTrackBasedReferencePoint: invalid track momentum " << totalTrackMomentumAtStart << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }
        else
        {
            referencePoint = referencePointAtStartWeighted * (1.f / totalTrackMomentumAtStart);
        }
    }
    else
    {
        if (totalTrackMomentumAtDca < std::numeric_limits<float>::epsilon())
        {
            std::cout<<"WARNING PfoCreator::CalculateTrackBasedReferencePoint: invalid track momentum " << totalTrackMomentumAtDca << std::endl;
            throw pandora::StatusCodeException(pandora::STATUS_CODE_FAILURE);
        }
        else
        {
            referencePoint = referencePointAtDCAWeighted * (1.f / totalTrackMomentumAtDca);
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PfoCreator::IsValidParentTrack(const pandora::Track *const pPandoraTrack, const pandora::TrackList &allTrackList) const
{
    const pandora::TrackList &parentTrackList(pPandoraTrack->GetParentList());

    for (pandora::TrackList::const_iterator iter = parentTrackList.begin(), iterEnd = parentTrackList.end(); iter != iterEnd; ++iter)
    {
        if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
            continue;

        // ATTN This track must have a parent not in the all track list; still use it if it is the closest to the ip
        std::cout<<"WARNING PfoCreator::IsValidParentTrack: mismatch in track relationship information, use information as available " << std::endl;

        if (this->IsClosestTrackToIP(pPandoraTrack, allTrackList))
            return true;

        return false;
    }

    // Ideal case: All parents are associated to same pfo
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PfoCreator::HasValidSiblingTrack(const pandora::Track *const pPandoraTrack, const pandora::TrackList &allTrackList) const
{
    const pandora::TrackList &siblingTrackList(pPandoraTrack->GetSiblingList());

    for (pandora::TrackList::const_iterator iter = siblingTrackList.begin(), iterEnd = siblingTrackList.end(); iter != iterEnd; ++iter)
    {
        if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
            continue;

        // ATTN This track must have a sibling not in the all track list; still use it if it has a second sibling that is in the list
        std::cout<<"WARNING PfoCreator::HasValidSiblingTrack: mismatch in track relationship information, use information as available " << std::endl;

        if (this->AreAnyOtherSiblingsInList(pPandoraTrack, allTrackList))
            return true;

        return false;
    }

    // Ideal case: All siblings associated to same pfo
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PfoCreator::IsClosestTrackToIP(const pandora::Track *const pPandoraTrack, const pandora::TrackList &allTrackList) const 
{
    const pandora::Track *pClosestTrack(NULL);
    float closestTrackDisplacement(std::numeric_limits<float>::max()); 

    for (pandora::TrackList::const_iterator iter = allTrackList.begin(), iterEnd = allTrackList.end(); iter != iterEnd; ++iter)
    {
        const pandora::Track *const pTrack(*iter);
        const float trialTrackDisplacement(pTrack->GetTrackStateAtStart().GetPosition().GetMagnitude());

        if (trialTrackDisplacement < closestTrackDisplacement)
        {
            closestTrackDisplacement = trialTrackDisplacement;
            pClosestTrack = pTrack;
        }
    }

    return (pPandoraTrack == pClosestTrack);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PfoCreator::AreAnyOtherSiblingsInList(const pandora::Track *const pPandoraTrack, const pandora::TrackList &allTrackList) const
{
    const pandora::TrackList &siblingTrackList(pPandoraTrack->GetSiblingList());

    for (pandora::TrackList::const_iterator iter = siblingTrackList.begin(), iterEnd = siblingTrackList.end(); iter != iterEnd; ++iter)
    {
        if (allTrackList.end() != std::find(allTrackList.begin(), allTrackList.end(), *iter))
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoCreator::SetRecoParticleReferencePoint(const pandora::CartesianVector &referencePoint, edm4hep::ReconstructedParticle *const pReconstructedParticle) const
{
    const float referencePointArray[3] = {referencePoint.GetX(), referencePoint.GetY(), referencePoint.GetZ()};
    pReconstructedParticle->setReferencePoint(referencePointArray);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoCreator::AddTracksToRecoParticle(const pandora::ParticleFlowObject *const pPandoraPfo, edm4hep::ReconstructedParticle *const pReconstructedParticle) const
{
    const pandora::TrackList &trackList(pPandoraPfo->GetTrackList());

    for (pandora::TrackList::const_iterator tIter = trackList.begin(), tIterEnd = trackList.end(); tIter != tIterEnd; ++tIter)
    {
        const pandora::Track *const pTrack(*tIter);
        const edm4hep::Track *const pLcioTrack0 = (edm4hep::Track*)(pTrack->GetParentAddress());
        const edm4hep::Track pLcioTrack = *pLcioTrack0;
        pReconstructedParticle->addToTracks(pLcioTrack);

    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PfoCreator::SetRecoParticlePropertiesFromPFO(const pandora::ParticleFlowObject *const pPandoraPfo, edm4hep::ReconstructedParticle *const pReconstructedParticle) const
{
    const float momentum[3] = {pPandoraPfo->GetMomentum().GetX(), pPandoraPfo->GetMomentum().GetY(), pPandoraPfo->GetMomentum().GetZ()};
    pReconstructedParticle->setMomentum(momentum);
    pReconstructedParticle->setEnergy(pPandoraPfo->GetEnergy());
    pReconstructedParticle->setMass(pPandoraPfo->GetMass());
    pReconstructedParticle->setCharge(pPandoraPfo->GetCharge());
    pReconstructedParticle->setType(pPandoraPfo->GetParticleId());
}

//------------------------------------------------------------------------------------------------------------------------------------------

PfoCreator::Settings::Settings():
    m_emStochasticTerm(0.17f),
    m_hadStochasticTerm(0.6f),
    m_emConstantTerm(0.01f),
    m_hadConstantTerm(0.03f)
{
}
