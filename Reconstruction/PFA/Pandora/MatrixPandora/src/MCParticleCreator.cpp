/**
 * 
 *  @brief  Implementation of the mc particle creator class.
 * 
 *  $Log: $
 */


#include "edm4hep/MCParticleConst.h"
#include "edm4hep/MCParticle.h" 
#include "edm4hep/MCRecoCaloAssociation.h" 
#include "edm4hep/SimCalorimeterHitConst.h" 
#include "edm4hep/CaloHitContributionConst.h" 
#include "edm4hep/Track.h" 
#include "edm4hep/MCRecoTrackerAssociation.h" 
#include "edm4hep/SimTrackerHitConst.h" 




#include "PandoraMatrixAlg.h"
#include "MCParticleCreator.h"

#include <cmath>
#include <limits>
#include <assert.h>

MCParticleCreator::MCParticleCreator(const Settings &settings, const pandora::Pandora *const pPandora) :
    m_settings(settings),
    m_pPandora(pPandora),
    m_bField(settings.m_bField)
{
m_id_pMC_map = new std::map<unsigned int, const edm4hep::MCParticle*>;
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCParticleCreator::~MCParticleCreator()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

pandora::StatusCode MCParticleCreator::CreateMCParticles(const CollectionMaps& collectionMaps ) const
{
    for (StringVector::const_iterator iter = m_settings.m_mcParticleCollections.begin(), iterEnd = m_settings.m_mcParticleCollections.end();
        iter != iterEnd; ++iter)
    {
        if(collectionMaps.collectionMap_MC.find(*iter) == collectionMaps.collectionMap_MC.end()) continue;
        try
        {
            const std::vector<edm4hep::MCParticle>& pMCParticleCollection = (collectionMaps.collectionMap_MC.find(*iter))->second;
            std::cout<<"Do CreateMCParticles, collection:"<<(*iter)<<", size="<<pMCParticleCollection.size()<<std::endl;
            for (int im = 0; im < pMCParticleCollection.size(); im++)
            {
                try
                {
                    const edm4hep::MCParticle& pMcParticle = pMCParticleCollection.at(im);
                    PandoraApi::MCParticle::Parameters mcParticleParameters;
                    mcParticleParameters.m_energy =   sqrt(pMcParticle.getMomentum()[0] * pMcParticle.getMomentum()[0] + pMcParticle.getMomentum()[1] * pMcParticle.getMomentum()[1] + pMcParticle.getMomentum()[2] * pMcParticle.getMomentum()[2] + pMcParticle.getMass() * pMcParticle.getMass());
                    mcParticleParameters.m_particleId = pMcParticle.getPDG();
                    mcParticleParameters.m_mcParticleType = pandora::MC_3D;
                    mcParticleParameters.m_pParentAddress = &pMcParticle;
                    unsigned int p_id = pMcParticle.id();
                    const edm4hep::MCParticle* p_mc = &pMcParticle;
                    (*m_id_pMC_map) [p_id]   = p_mc;
                    mcParticleParameters.m_momentum = pandora::CartesianVector(pMcParticle.getMomentum()[0], pMcParticle.getMomentum()[1],
                        pMcParticle.getMomentum()[2]);
                    mcParticleParameters.m_vertex = pandora::CartesianVector(pMcParticle.getVertex()[0], pMcParticle.getVertex()[1],
                        pMcParticle.getVertex()[2]);
                    mcParticleParameters.m_endpoint = pandora::CartesianVector(pMcParticle.getEndpoint()[0], pMcParticle.getEndpoint()[1],
                        pMcParticle.getEndpoint()[2]);

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*m_pPandora, mcParticleParameters));

                    // Create parent-daughter relationships
                    for(std::vector<edm4hep::ConstMCParticle>::const_iterator itDaughter = pMcParticle.daughters_begin(),
                        itDaughterEnd = pMcParticle.daughters_end(); itDaughter != itDaughterEnd; ++itDaughter)
                    {   
                        for (int ida = 0; ida < pMCParticleCollection.size(); ida++)
                        {
                           if(pMCParticleCollection.at(ida).id()==(*itDaughter).id())
                           {
                              
                                const edm4hep::MCParticle& dMcParticle = pMCParticleCollection.at(ida);
                                if(&pMcParticle == &dMcParticle){std::cout<< "error, mother and daughter are the same mc particle, don't save SetMCParentDaughterRelationship"<<std::endl;}
                                else PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetMCParentDaughterRelationship(*m_pPandora, &pMcParticle, &dMcParticle));
                                break;
                           }
                        }
                    }
                    
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    std::cout << "Failed to extract MCParticle: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    std::cout << "Failed to extract MCParticle: " <<  std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout << "Failed to extract MCParticles collection: " << *iter << ", " <<  std::endl;
        }
    }
    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
/*
pandora::StatusCode MCParticleCreator::CreateTrackToMCParticleRelationships(const EVENT::LCEvent *const pLCEvent, const TrackVector &trackVector) const
{
    for (StringVector::const_iterator iter = m_settings.m_lcTrackRelationCollections.begin(), iterEnd = m_settings.m_lcTrackRelationCollections.end();
         iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pMCRelationCollection = pLCEvent->getCollection(*iter);
            UTIL::LCRelationNavigator navigate(pMCRelationCollection);

            for (TrackVector::const_iterator trackIter = trackVector.begin(), trackIterEnd = trackVector.end();
                trackIter != trackIterEnd; ++trackIter)
            {
                try
                {
                    EVENT::Track *pTrack = *trackIter;
                    const EVENT::LCObjectVec &objectVec = navigate.getRelatedToObjects(*trackIter);

                    // Get reconstructed momentum at dca
                    const pandora::Helix helixFit(pTrack->getPhi(), pTrack->getD0(), pTrack->getZ0(), pTrack->getOmega(), pTrack->getTanLambda(), m_bField);
                    const float recoMomentum(helixFit.GetMomentum().GetMagnitude());

                    // Use momentum magnitude to identify best mc particle
                    MCParticle *pBestMCParticle = NULL;
                    float bestDeltaMomentum(std::numeric_limits<float>::max());

                    for (EVENT::LCObjectVec::const_iterator itRel = objectVec.begin(), itRelEnd = objectVec.end(); itRel != itRelEnd; ++itRel)
                    {
                        EVENT::MCParticle *pMCParticle = NULL;
                        pMCParticle = dynamic_cast<MCParticle *>(*itRel);

                        if (NULL == pMCParticle)
                            continue;

                        const float trueMomentum(pandora::CartesianVector(pMCParticle->getMomentum()[0], pMCParticle->getMomentum()[1],
                            pMCParticle->getMomentum()[2]).GetMagnitude());

                        const float deltaMomentum(std::fabs(recoMomentum - trueMomentum));

                        if (deltaMomentum < bestDeltaMomentum)
                        {
                            pBestMCParticle = pMCParticle;
                            bestDeltaMomentum = deltaMomentum;
                        }
                    }

                    if (NULL == pBestMCParticle)
                        continue;

                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackToMCParticleRelationship(*m_pPandora, pTrack,
                        pBestMCParticle));
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract track to mc particle relationship: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract track to mc particle relationship: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(DEBUG5) << "Failed to extract track to mc particle relationships collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}
*/
//------------------------------------------------------------------------------------------------------------------------------------------
/*
pandora::StatusCode MCParticleCreator::CreateCaloHitToMCParticleRelationships(const EVENT::LCEvent *const pLCEvent, const CalorimeterHitVector &calorimeterHitVector) const
{
    typedef std::map<MCParticle *, float> MCParticleToEnergyWeightMap;
    MCParticleToEnergyWeightMap mcParticleToEnergyWeightMap;

    for (StringVector::const_iterator iter = m_settings.m_lcCaloHitRelationCollections.begin(), iterEnd = m_settings.m_lcCaloHitRelationCollections.end();
         iter != iterEnd; ++iter)
    {
        try
        {
            const EVENT::LCCollection *pMCRelationCollection = pLCEvent->getCollection(*iter);
            UTIL::LCRelationNavigator navigate(pMCRelationCollection);

            for (CalorimeterHitVector::const_iterator caloHitIter = calorimeterHitVector.begin(),
                caloHitIterEnd = calorimeterHitVector.end(); caloHitIter != caloHitIterEnd; ++caloHitIter)
            {
                try
                {
                    mcParticleToEnergyWeightMap.clear();
                    const EVENT::LCObjectVec &objectVec = navigate.getRelatedToObjects(*caloHitIter);

                    for (EVENT::LCObjectVec::const_iterator itRel = objectVec.begin(), itRelEnd = objectVec.end(); itRel != itRelEnd; ++itRel)
                    {
                        EVENT::SimCalorimeterHit *pSimHit = dynamic_cast<SimCalorimeterHit *>(*itRel);

                        if (NULL == pSimHit)
                            continue;

                        for (int iCont = 0, iEnd = pSimHit->getNMCContributions(); iCont < iEnd; ++iCont)
                        {
                            mcParticleToEnergyWeightMap[pSimHit->getParticleCont(iCont)] += pSimHit->getEnergyCont(iCont);
                        }
                    }

                    for (MCParticleToEnergyWeightMap::const_iterator mcParticleIter = mcParticleToEnergyWeightMap.begin(),
                        mcParticleIterEnd = mcParticleToEnergyWeightMap.end(); mcParticleIter != mcParticleIterEnd; ++mcParticleIter)
                    {
                        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(*m_pPandora,
                            *caloHitIter, mcParticleIter->first, mcParticleIter->second));
                    }
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    streamlog_out(ERROR) << "Failed to extract calo hit to mc particle relationship: " << statusCodeException.ToString() << std::endl;
                }
                catch (EVENT::Exception &exception)
                {
                    streamlog_out(WARNING) << "Failed to extract calo hit to mc particle relationship: " << exception.what() << std::endl;
                }
            }
        }
        catch (EVENT::Exception &exception)
        {
            streamlog_out(DEBUG5) << "Failed to extract calo hit to mc particle relationships collection: " << *iter << ", " << exception.what() << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}
*/

//-------------- using sim calo hit and digi calo hit, no weight here---------------------------------------------------------------------------------------------------------------
pandora::StatusCode MCParticleCreator::CreateCaloHitToMCParticleRelationships(const CollectionMaps& collectionMaps, const CalorimeterHitVector &calorimeterHitVector) const
{
    std::cout<<"Do CreateCaloHitToMCParticleRelationships"<<std::endl;
    typedef std::map<const edm4hep::MCParticle *, float> MCParticleToEnergyWeightMap;
    MCParticleToEnergyWeightMap mcParticleToEnergyWeightMap;

    for (StringVector::const_iterator iter = m_settings.m_CaloHitRelationCollections.begin(), iterEnd = m_settings.m_CaloHitRelationCollections.end();
         iter != iterEnd; ++iter)
    {
        if(collectionMaps.collectionMap_CaloRel.find(*iter) == collectionMaps.collectionMap_CaloRel.end()) continue;
        try
        {
            const std::vector<edm4hep::MCRecoCaloAssociation>& pMCRecoCaloAssociationCollection = (collectionMaps.collectionMap_CaloRel.find(*iter))->second;

            for (unsigned i_calo=0; i_calo < calorimeterHitVector.size(); i_calo++)
            {
                try
                {
                    mcParticleToEnergyWeightMap.clear();
                    for(unsigned ic=0; ic < pMCRecoCaloAssociationCollection.size(); ic++)
                    {
                        if( pMCRecoCaloAssociationCollection.at(ic).getRec().id() != (*(calorimeterHitVector.at(i_calo))).id() ) continue;
                        
                        const edm4hep::ConstSimCalorimeterHit pSimHit = pMCRecoCaloAssociationCollection.at(ic).getSim();
                        for (int iCont = 0, iEnd = pSimHit.contributions_size(); iCont < iEnd; ++iCont)
                        {
                            edm4hep::ConstCaloHitContribution conb = pSimHit.getContributions(iCont);
                            const edm4hep::ConstMCParticle ipa = conb.getParticle();
                            float  ien = conb.getEnergy();
                            if( m_id_pMC_map->find(ipa.id()) == m_id_pMC_map->end() ) continue;
                            const edm4hep::MCParticle * p_tmp = (*m_id_pMC_map)[ipa.id()]; 
                            mcParticleToEnergyWeightMap[p_tmp] += ien;
                        }
                        
                    }

                    for (MCParticleToEnergyWeightMap::const_iterator mcParticleIter = mcParticleToEnergyWeightMap.begin(),
                        mcParticleIterEnd = mcParticleToEnergyWeightMap.end(); mcParticleIter != mcParticleIterEnd; ++mcParticleIter)
                    {
                        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetCaloHitToMCParticleRelationship(*m_pPandora,
                            calorimeterHitVector.at(i_calo), mcParticleIter->first, mcParticleIter->second));
                    }
                }
                catch (pandora::StatusCodeException &statusCodeException)
                {
                    std::cout<<"ERROR Failed to extract calo hit to mc particle relationship: " << statusCodeException.ToString() << std::endl;
                }
                catch (...)
                {
                    std::cout<<"WARNING Failed to extract calo hit to mc particle relationship " << std::endl;
                }
            }
        }
        catch (...)
        {
            std::cout<<"DEBUG5 Failed to extract calo hit to mc particle relationships collection: " << *iter << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}


//-------------- using sim tracker hit and tracker hit, no weight here---------------------------------------------------------------------------------------------------------------
pandora::StatusCode MCParticleCreator::CreateTrackToMCParticleRelationships(const CollectionMaps& collectionMaps, const TrackVector &trackVector) const
{
    std::cout<<"Do CreateTrackToMCParticleRelationships"<<std::endl;
    for (unsigned ik = 0; ik < trackVector.size(); ik++)
    {
        const edm4hep::Track *pTrack = trackVector.at(ik);
        // Get reconstructed momentum at dca
        const pandora::Helix helixFit(pTrack->getTrackStates(0).phi, pTrack->getTrackStates(0).D0, pTrack->getTrackStates(0).Z0, pTrack->getTrackStates(0).omega, pTrack->getTrackStates(0).tanLambda, m_bField);
        const float recoMomentum(helixFit.GetMomentum().GetMagnitude());
        // Use momentum magnitude to identify best mc particle
        edm4hep::MCParticle *pBestMCParticle = NULL;
        float bestDeltaMomentum(std::numeric_limits<float>::max());
        try
        {
            for (StringVector::const_iterator iter = m_settings.m_TrackRelationCollections.begin(), iterEnd = m_settings.m_TrackRelationCollections.end(); iter != iterEnd; ++iter)
            {
                if(collectionMaps.collectionMap_TrkRel.find(*iter) == collectionMaps.collectionMap_TrkRel.end()) continue;
                const std::vector<edm4hep::MCRecoTrackerAssociation>& pMCRecoTrackerAssociationCollection = (collectionMaps.collectionMap_TrkRel.find(*iter))->second;
                for(unsigned ith=0 ; ith<pTrack->trackerHits_size(); ith++)
                {
                    for(unsigned ic=0; ic < pMCRecoTrackerAssociationCollection.size(); ic++)
                    {
                        if( pMCRecoTrackerAssociationCollection.at(ic).getRec().id() != pTrack->getTrackerHits(ith).id() ) continue;
                        const edm4hep::ConstSimTrackerHit pSimHit = pMCRecoTrackerAssociationCollection.at(ic).getSim();
                        const edm4hep::ConstMCParticle ipa = pSimHit.getMCParticle();
                        if( m_id_pMC_map->find(ipa.id()) == m_id_pMC_map->end() ) continue;
                        const float trueMomentum(pandora::CartesianVector(ipa.getMomentum()[0], ipa.getMomentum()[1], ipa.getMomentum()[2]).GetMagnitude());
                        const float deltaMomentum(std::fabs(recoMomentum - trueMomentum));
                        if (deltaMomentum < bestDeltaMomentum)
                        {
                            pBestMCParticle =const_cast<edm4hep::MCParticle*>((*m_id_pMC_map)[ipa.id()]);
                            bestDeltaMomentum = deltaMomentum;
                        }
                    }
                }
            }
            
            if (NULL == pBestMCParticle) continue;
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::SetTrackToMCParticleRelationship(*m_pPandora, pTrack, pBestMCParticle));
        }
        catch (pandora::StatusCodeException &statusCodeException)
        {
            std::cout<<"ERROR Failed to extract track to mc particle relationship: " << statusCodeException.ToString() << std::endl;
        }
        catch (...)
        {
            std::cout<<"WARNING Failed to extract track to mc particle relationship: " << std::endl;
        }
    }

    return pandora::STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

MCParticleCreator::Settings::Settings()
{
}
