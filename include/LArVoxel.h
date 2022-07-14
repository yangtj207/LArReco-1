/**
 *  @file   LArReco/include/LArRay.h
 *
 *  @brief  Header file for LArVoxel
 *
 *  $Log: $
 */
#ifndef PANDORA_LAR_VOXEL_H
#define PANDORA_LAR_VOXEL_H 1

#include "Pandora/PandoraInputTypes.h"
#include <vector>

namespace lar_nd_reco
{

class LArVoxel
{
public:
    /**
     *  @brief  Constructor
     *
     *  @param  voxelID Total bin number for the voxel (long integer, since it can be > 2^31)
     *  @param  energyInVoxel The total deposited energy in the voxel (GeV)
     *  @param  voxelPosVect Voxel position, set as the first corner of the voxel bin
     *  @param  trackID The Geant4 ID of the (first) contributing track to this voxel
     */
    LArVoxel(const long voxelID, const float energyInVoxel, const pandora::CartesianVector &voxelPosVect, const int trackID);

    /**
     *  @brief  Set voxel energy (GeV)
     *
     *  @param  E voxel energy
     */
    void SetEnergy(const float E);

    long m_voxelID;                          ///< The long integer ID of the voxel (can be larger than 2^31)
    float m_energyInVoxel;                   ///< The energy in the voxel (GeV)
    pandora::CartesianVector m_voxelPosVect; ///< Position vector (x,y,z) of the first voxel corner
    int m_trackID;                           ///< The Geant4 ID of the (first) contributing track to this voxel
};

typedef std::vector<LArVoxel> LArVoxelList;

//------------------------------------------------------------------------------------------------------------------------------------------

inline LArVoxel::LArVoxel(const long voxelID, const float energyInVoxel, const pandora::CartesianVector &voxelPosVect, const int trackID) :
    m_voxelID(voxelID),
    m_energyInVoxel(energyInVoxel),
    m_voxelPosVect(voxelPosVect),
    m_trackID(trackID)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void LArVoxel::SetEnergy(const float E)
{
    m_energyInVoxel = E;
}

} // namespace lar_nd_reco

#endif
