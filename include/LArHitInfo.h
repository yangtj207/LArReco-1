/**
 *  @file   LArReco/include/LArHitInfo.h
 *
 *  @brief  Header file for storing hit information (from various formats)
 *
 *  $Log: $
 */
#ifndef PANDORA_LAR_HIT_INFO_H
#define PANDORA_LAR_HIT_INFO_H 1

#include "Pandora/PandoraInputTypes.h"
#include "TG4HitSegment.h"
#include "TLorentzVector.h"

namespace lar_nd_reco
{

class LArHitInfo
{
public:
    /**
     *  @brief  Constructor using TG4HitSegments
     *
     *  @param  g4Hit The Geant4 hit segment (step)
     *  @param  lengthScale Scaling factor to use cm length dimensions
     *  @param  energyScale Scaling factor to use GeV energies
     */
    LArHitInfo(const TG4HitSegment &g4Hit, const float lengthScale, const float energyScale);

    /**
     *  @brief  Constructor using CartesianVectors
     *
     *  @param  g4Hit The Geant4 hit segment (step)
     *  @param  lengthScale Scaling factor to use cm length dimensions
     *  @param  energyScale Scaling factor to use GeV energies
     */
    LArHitInfo(const pandora::CartesianVector &start, const pandora::CartesianVector &stop, const float energy, const int trackID,
        const float lengthScale, const float energyScale);

    pandora::CartesianVector m_start; ///< Starting point of the hit step
    pandora::CartesianVector m_stop;  ///< End point of the hit step
    float m_energy;                   ///< The total energy of the step
    int m_trackID;                    ///< The ID of the (main) contributing particle
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArHitInfo::LArHitInfo(const TG4HitSegment &g4Hit, const float lengthScale, const float energyScale) :
    m_start({0.f, 0.f, 0.f}),
    m_stop({0.f, 0.f, 0.f}),
    m_energy(0.f),
    m_trackID(0)
{

    // Start and end positions
    const TLorentzVector &hitStart = g4Hit.GetStart() * lengthScale;
    const TLorentzVector &hitStop = g4Hit.GetStop() * lengthScale;

    m_start.SetValues(hitStart.X(), hitStart.Y(), hitStart.Z());
    m_stop.SetValues(hitStop.X(), hitStop.Y(), hitStop.Z());

    // Hit segment total energy
    m_energy = g4Hit.GetEnergyDeposit() * energyScale;

    // Get the trackID of the (main) contributing particle.
    // ATTN: this can very rarely be more than one track
    m_trackID = g4Hit.GetContributors()[0];
}

inline LArHitInfo::LArHitInfo(const pandora::CartesianVector &start, const pandora::CartesianVector &stop, const float energy,
    const int trackID, const float lengthScale, const float energyScale) :
    m_start(start * lengthScale),
    m_stop(stop * lengthScale),
    m_energy(energy * energyScale),
    m_trackID(trackID)
{
}
} // namespace lar_nd_reco

#endif
