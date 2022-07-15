/**
 *  @file   LArRecoMP/test/PandoraInterface.cc
 *
 *  @brief  Implementation of the lar reco mp application
 *
 *  $Log: $
 */

#include "TFile.h"
#include "TTree.h"

#include "TGeoBBox.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoShape.h"
#include "TGeoVolume.h"

#include "TG4PrimaryVertex.h"

#include "Api/PandoraApi.h"
#include "Geometry/LArTPC.h"
#include "Helpers/XmlHelper.h"
#include "Managers/GeometryManager.h"
#include "Managers/PluginManager.h"
#include "Xml/tinyxml.h"

#include "larpandoracontent/LArContent.h"
#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"
#include "larpandoracontent/LArControlFlow/MultiPandoraApi.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArMCParticle.h"
#include "larpandoracontent/LArPlugins/LArPseudoLayerPlugin.h"
#include "larpandoracontent/LArPlugins/LArRotationalTransformationPlugin.h"

#include "LArRay.h"
#include "PandoraInterface.h"

#ifdef MONITORING
#include "TApplication.h"
#endif

#include <getopt.h>
#include <iostream>
#include <memory>
#include <random>
#include <string>

using namespace pandora;
using namespace lar_nd_reco;

int main(int argc, char *argv[])
{
    int errorNo(0);
    const Pandora *pPrimaryPandora(nullptr);

    try
    {
        Parameters parameters;

        if (!ParseCommandLine(argc, argv, parameters))
            return 1;

#ifdef MONITORING
        TApplication *pTApplication = new TApplication("LArReco", &argc, argv);
        pTApplication->SetReturnFromRun(kTRUE);
#endif

        pPrimaryPandora = new pandora::Pandora();

        if (!pPrimaryPandora)
            throw StatusCodeException(STATUS_CODE_FAILURE);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterAlgorithms(*pPrimaryPandora));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArContent::RegisterBasicPlugins(*pPrimaryPandora));

        MultiPandoraApi::AddPrimaryPandoraInstance(pPrimaryPandora);

        CreateGeometry(parameters, pPrimaryPandora);
        ProcessExternalParameters(parameters, pPrimaryPandora);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetPseudoLayerPlugin(*pPrimaryPandora, new lar_content::LArPseudoLayerPlugin));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraApi::SetLArTransformationPlugin(*pPrimaryPandora, new lar_content::LArRotationalTransformationPlugin));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ReadSettings(*pPrimaryPandora, parameters.m_settingsFile));

        ProcessEvents(parameters, pPrimaryPandora);
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cerr << "Pandora StatusCodeException: " << statusCodeException.ToString() << statusCodeException.GetBackTrace() << std::endl;
        errorNo = 1;
    }
    catch (...)
    {
        std::cerr << "Unknown exception: " << std::endl;
        errorNo = 1;
    }

    MultiPandoraApi::DeletePandoraInstances(pPrimaryPandora);
    return errorNo;
}

//------------------------------------------------------------------------------------------------------------------------------------------

namespace lar_nd_reco
{

void CreateGeometry(const Parameters &parameters, const Pandora *const pPrimaryPandora)
{
    // Get the geometry info from the appropriate ROOT file
    TFile *fileSource = TFile::Open(parameters.m_geomFileName.c_str(), "READ");
    if (!fileSource)
    {
        std::cout << "Error in CreateGeometry(): can't open file " << parameters.m_geomFileName << std::endl;
        return;
    }

    TGeoManager *pSimGeom = dynamic_cast<TGeoManager *>(fileSource->Get(parameters.m_geomManagerName.c_str()));
    if (!pSimGeom)
    {
        std::cout << "Could not find the geometry manager named " << parameters.m_geomManagerName << std::endl;
        return;
    }

    // Start by looking at the top level volume and move down to the one we need
    std::string name;
    const std::string neededNode(parameters.m_geometryVolName);
    TGeoNode *pCurrentNode = pSimGeom->GetCurrentNode();
    bool foundNode(false);

    // Initialise volume matrix using the master (top) volume.
    // This is updated as we go down the volume hierarchy
    std::unique_ptr<TGeoHMatrix> pVolMatrix = std::make_unique<TGeoHMatrix>(*pCurrentNode->GetMatrix());

    // Maximum number of nodes to search through, same as default TGeoManager counting node limit
    const int maxNodes(10000);
    int iNode(0);
    while (foundNode == false && iNode < maxNodes)
    {
        pCurrentNode = pSimGeom->GetCurrentNode();
        iNode++;
        name = pCurrentNode->GetName();
        std::unique_ptr<TGeoHMatrix> pCurrentMatrix = std::make_unique<TGeoHMatrix>(*pCurrentNode->GetMatrix());
        pVolMatrix->Multiply(pCurrentMatrix.get());

        int i1 = 0;
        for (int i = 0; i < pCurrentNode->GetNdaughters(); i++)
        {

            pSimGeom->CdDown(i1);
            TGeoNode *pNode = pSimGeom->GetCurrentNode();
            name = pNode->GetName();
            std::unique_ptr<TGeoHMatrix> pMatrix = std::make_unique<TGeoHMatrix>(*pNode->GetMatrix());

            if (name == neededNode)
            {
                foundNode = true;
                pVolMatrix->Multiply(pMatrix.get());
                break;
            }
            else if (i + 1 != pCurrentNode->GetNdaughters())
            {
                pSimGeom->CdUp();
                i1++;
            }
        }

        if (foundNode)
            break;
    }

    if (!foundNode)
    {
        std::cout << "Could not find the required placement geometry volume " << neededNode << std::endl;
        return;
    }

    // The current node should now be the placement volume we need to set the geometry parameters
    pCurrentNode = pSimGeom->GetCurrentNode();
    name = pCurrentNode->GetName();
    std::cout << "Current Node: " << name << std::endl;
    std::cout << "Current N daughters: " << pCurrentNode->GetVolume()->GetNdaughters() << std::endl;
    std::cout << "Current transformation matrix:" << std::endl;
    pVolMatrix->Print();

    // Get the BBox dimensions from the placement volume, which is assumed to be a cube
    TGeoVolume *pCurrentVol = pCurrentNode->GetVolume();
    TGeoShape *pCurrentShape = pCurrentVol->GetShape();
    pCurrentShape->InspectShape();
    TGeoBBox *pBox = dynamic_cast<TGeoBBox *>(pCurrentShape);

    // Now can get origin/width data from the BBox
    const float dx = pBox->GetDX() * parameters.m_lengthScale; // Note these are the half widths
    const float dy = pBox->GetDY() * parameters.m_lengthScale;
    const float dz = pBox->GetDZ() * parameters.m_lengthScale;
    const double *pOrigin = pBox->GetOrigin();

    std::cout << "Origin = (" << pOrigin[0] << ", " << pOrigin[1] << ", " << pOrigin[2] << ")" << std::endl;

    // Translate local origin to global coordinates
    double level1[3] = {0.0, 0.0, 0.0};
    pCurrentNode->LocalToMasterVect(pOrigin, level1);

    std::cout << "Level1 = (" << level1[0] << ", " << level1[1] << ", " << level1[2] << ")" << std::endl;

    // Can now create a geometry using the found parameters
    PandoraApi::Geometry::LArTPC::Parameters geoparameters;

    try
    {
        const double *pVolTrans = pVolMatrix->GetTranslation();
        geoparameters.m_centerX = (level1[0] + pVolTrans[0]) * parameters.m_lengthScale;
        geoparameters.m_centerY = (level1[1] + pVolTrans[1]) * parameters.m_lengthScale;
        geoparameters.m_centerZ = (level1[2] + pVolTrans[2]) * parameters.m_lengthScale;
        geoparameters.m_widthX = dx * 2.0;
        geoparameters.m_widthY = dy * 2.0;
        geoparameters.m_widthZ = dz * 2.0;

        // ATTN: parameters past here taken from uboone
        geoparameters.m_larTPCVolumeId = 0;
        geoparameters.m_wirePitchU = 0.300000011921;
        geoparameters.m_wirePitchV = 0.300000011921;
        geoparameters.m_wirePitchW = 0.300000011921;
        geoparameters.m_wireAngleU = 1.04719758034;
        geoparameters.m_wireAngleV = -1.04719758034;
        geoparameters.m_wireAngleW = 0.0;
        geoparameters.m_sigmaUVW = 1;
        geoparameters.m_isDriftInPositiveX = 0;
    }
    catch (const pandora::StatusCodeException &)
    {
        std::cout << "CreatePandoraLArTPCs - invalid tpc parameter provided" << std::endl;
    }

    try
    {
        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::Geometry::LArTPC::Create(*pPrimaryPandora, geoparameters));
    }
    catch (const pandora::StatusCodeException &)
    {
        std::cout << "CreatePandoraLArTPCs - unable to create tpc, insufficient or "
                     "invalid information supplied"
                  << std::endl;
    }

    fileSource->Close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessEvents(const Parameters &parameters, const Pandora *const pPrimaryPandora)
{
    if (parameters.m_dataFormat == Parameters::LArNDFormat::SED)
    {

        ProcessSEDEvents(parameters, pPrimaryPandora);
    }
    else
    {
        ProcessEDepSimEvents(parameters, pPrimaryPandora);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessEDepSimEvents(const Parameters &parameters, const Pandora *const pPrimaryPandora)
{

    TFile *fileSource = TFile::Open(parameters.m_inputFileName.c_str(), "READ");
    if (!fileSource)
    {
        std::cout << "Error in ProcessEDepSimEvents(): can't open file " << parameters.m_inputFileName << std::endl;
        return;
    }

    TTree *pEDepSimTree = dynamic_cast<TTree *>(fileSource->Get(parameters.m_inputTreeName.c_str()));
    if (!pEDepSimTree)
    {
        std::cout << "Could not find the event tree " << parameters.m_inputTreeName << std::endl;
        fileSource->Close();
        return;
    }

    TG4Event *pEDepSimEvent(nullptr);
    pEDepSimTree->SetBranchAddress("Event", &pEDepSimEvent);

    // Factory for creating LArCaloHits
    lar_content::LArCaloHitFactory m_larCaloHitFactory;

    // Detector volume for voxelising the hits
    const GeometryManager *geom = pPrimaryPandora->GetGeometry();
    const LArTPC &tpc = geom->GetLArTPC();

    const float botX = tpc.GetCenterX() - 0.5 * tpc.GetWidthX();
    const float botY = tpc.GetCenterY() - 0.5 * tpc.GetWidthY();
    const float botZ = tpc.GetCenterZ() - 0.5 * tpc.GetWidthZ();
    const float topX = botX + tpc.GetWidthX();
    const float topY = botY + tpc.GetWidthY();
    const float topZ = botZ + tpc.GetWidthZ();
    const float voxelWidth(parameters.m_voxelWidth);
    const LArGrid grid(pandora::CartesianVector(botX, botY, botZ), pandora::CartesianVector(topX, topY, topZ),
        pandora::CartesianVector(voxelWidth, voxelWidth, voxelWidth));

    std::cout << "Total grid volume: bot = " << grid.m_bottom << "\n top = " << grid.m_top << std::endl;
    std::cout << "Making voxels with size " << grid.m_binWidths << std::endl;

    // Total number of entries in the TTree
    const int nEntries(pEDepSimTree->GetEntries());

    // Starting event
    const int startEvt = parameters.m_nEventsToSkip > 0 ? parameters.m_nEventsToSkip : 0;
    // Number of events to process, up to nEntries
    const int nProcess = parameters.m_nEventsToProcess > 0 ? parameters.m_nEventsToProcess : nEntries;
    // End event, up to nEntries
    const int endEvt = (startEvt + nProcess) < nEntries ? startEvt + nProcess : nEntries;

    std::cout << "Start event is " << startEvt << " and end event is " << endEvt - 1 << std::endl;

    for (int iEvt = startEvt; iEvt < endEvt; iEvt++)
    {
        if (parameters.m_shouldDisplayEventNumber)
            std::cout << std::endl << "   PROCESSING EVENT: " << iEvt << std::endl << std::endl;

        pEDepSimTree->GetEntry(iEvt);

        if (!pEDepSimEvent)
            return;

        int hitCounter(0);

        // Create MCParticles from Geant4 trajectories
        const MCParticleEnergyMap MCEnergyMap = CreateEDepSimMCParticles(*pEDepSimEvent, pPrimaryPandora, parameters);

        // Loop over (EDep) hits, which are stored in the hit segment detectors.
        // Only process hits from the detector we are interested in
        for (TG4HitSegmentDetectors::iterator detector = pEDepSimEvent->SegmentDetectors.begin();
             detector != pEDepSimEvent->SegmentDetectors.end(); ++detector)
        {
            if (detector->first != parameters.m_sensitiveDetName)
            {
                std::cout << "Skipping sensitive detector " << detector->first << "; expecting " << parameters.m_sensitiveDetName << std::endl;
                continue;
            }

            std::cout << "Show hits for " << detector->first << " (" << detector->second.size() << " hits)" << std::endl;
            std::cout << "                                 " << std::endl;

            LArVoxelList voxelList;

            // Loop over hit segments and create voxels from them
            for (TG4HitSegment &g4Hit : detector->second)
            {
                const LArHitInfo hitInfo(g4Hit, parameters.m_lengthScale, parameters.m_energyScale);
                const LArVoxelList currentVoxelList = MakeVoxels(hitInfo, grid, parameters);

                for (const LArVoxel &voxel : currentVoxelList)
                    voxelList.emplace_back(voxel);
            }

            std::cout << "Produced " << voxelList.size() << " voxels from " << detector->second.size() << " hit segments." << std::endl;

            // Merge voxels with the same IDs
            const LArVoxelList mergedVoxels = MergeSameVoxels(voxelList);

            std::cout << "Produced " << mergedVoxels.size() << " merged voxels from " << voxelList.size() << " voxels." << std::endl;

            // Stop processing the event if we have too many voxels: reco takes too long
            if (parameters.m_maxMergedVoxels > 0 && mergedVoxels.size() > parameters.m_maxMergedVoxels)
            {
                std::cout << "SKIPPING EVENT: number of merged voxels " << mergedVoxels.size() << " > " << parameters.m_maxMergedVoxels << std::endl;
                break;
            }

            // Loop over the voxels and make them into caloHits
            for (const LArVoxel &voxel : mergedVoxels)
            {
                const pandora::CartesianVector voxelPos(voxel.m_voxelPosVect);
                const float voxelE = voxel.m_energyInVoxel;
                const float MipE = 0.00075;
                const float voxelMipEquivalentE = voxelE / MipE;

                if (voxelMipEquivalentE > parameters.m_minVoxelMipEquivE)
                {
                    lar_content::LArCaloHitParameters caloHitParameters;
                    caloHitParameters.m_positionVector = voxelPos;
                    caloHitParameters.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f);
                    caloHitParameters.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
                    caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
                    caloHitParameters.m_cellSize0 = voxelWidth;
                    caloHitParameters.m_cellSize1 = voxelWidth;
                    caloHitParameters.m_cellThickness = voxelWidth;
                    caloHitParameters.m_nCellRadiationLengths = 1.f;
                    caloHitParameters.m_nCellInteractionLengths = 1.f;
                    caloHitParameters.m_time = 0.f;
                    caloHitParameters.m_inputEnergy = voxelE;
                    caloHitParameters.m_mipEquivalentEnergy = voxelMipEquivalentE;
                    caloHitParameters.m_electromagneticEnergy = voxelE;
                    caloHitParameters.m_hadronicEnergy = voxelE;
                    caloHitParameters.m_isDigital = false;
                    caloHitParameters.m_hitType = pandora::TPC_3D;
                    caloHitParameters.m_hitRegion = pandora::SINGLE_REGION;
                    caloHitParameters.m_layer = 0;
                    caloHitParameters.m_isInOuterSamplingLayer = false;
                    caloHitParameters.m_pParentAddress = (void *)(static_cast<uintptr_t>(++hitCounter));
                    caloHitParameters.m_larTPCVolumeId = 0;
                    caloHitParameters.m_daughterVolumeId = 0;

                    if (parameters.m_use3D)
                        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitParameters, m_larCaloHitFactory));

                    if (parameters.m_useLArTPC)
                    {
                        // Create U, V and W views assuming x is the common drift coordinate
                        const float x0_cm(voxelPos.GetX());
                        const float y0_cm(voxelPos.GetY());
                        const float z0_cm(voxelPos.GetZ());

                        lar_content::LArCaloHitParameters caloHitPars_UView(caloHitParameters);
                        caloHitPars_UView.m_hitType = pandora::TPC_VIEW_U;
                        const float upos_cm(pPrimaryPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoU(y0_cm, z0_cm));
                        caloHitPars_UView.m_positionVector = pandora::CartesianVector(x0_cm, 0.f, upos_cm);

                        lar_content::LArCaloHitParameters caloHitPars_VView(caloHitParameters);
                        caloHitPars_VView.m_hitType = pandora::TPC_VIEW_V;
                        const float vpos_cm(pPrimaryPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoV(y0_cm, z0_cm));
                        caloHitPars_VView.m_positionVector = pandora::CartesianVector(x0_cm, 0.f, vpos_cm);

                        lar_content::LArCaloHitParameters caloHitPars_WView(caloHitParameters);
                        caloHitPars_WView.m_hitType = pandora::TPC_VIEW_W;
                        const float wpos_cm(pPrimaryPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoW(y0_cm, z0_cm));
                        caloHitPars_WView.m_positionVector = pandora::CartesianVector(x0_cm, 0.f, wpos_cm);

                        // Create LArCaloHits for U, V and W views
                        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitPars_UView, m_larCaloHitFactory));
                        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitPars_VView, m_larCaloHitFactory));
                        PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                            PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitPars_WView, m_larCaloHitFactory));
                    }

                    // Set calo hit voxel to MCParticle relation using trackID
                    const int trackID = voxel.m_trackID;

                    // Find the energy fraction: voxelHitE/MCParticleE
                    float energyFrac(0.f), MCEnergy(0.f);
                    MCParticleEnergyMap::const_iterator mapIter = MCEnergyMap.find(trackID);
                    if (mapIter != MCEnergyMap.end())
                        MCEnergy = mapIter->second;

                    if (MCEnergy > 0.0)
                        energyFrac = voxelE / MCEnergy;

                    PandoraApi::SetCaloHitToMCParticleRelationship(
                        *pPrimaryPandora, (void *)((intptr_t)hitCounter), (void *)((intptr_t)trackID), energyFrac);
                }
            } // end voxel loop

        } // end segment detector loop

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pPrimaryPandora));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pPrimaryPandora));
    }

    // Close input file
    fileSource->Close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessSEDEvents(const Parameters &parameters, const Pandora *const pPrimaryPandora)
{

    TFile *fileSource = TFile::Open(parameters.m_inputFileName.c_str(), "READ");
    if (!fileSource)
    {
        std::cout << "Error in ProcessSEDEvents(): can't open file " << parameters.m_inputFileName << std::endl;
        return;
    }

    TTree *ndsim = dynamic_cast<TTree *>(fileSource->Get(parameters.m_inputTreeName.c_str()));
    if (!ndsim)
    {
        std::cout << "Could not find the event tree " << parameters.m_inputTreeName << std::endl;
        fileSource->Close();
        return;
    }

    const LArSED larsed(ndsim);

    // Factory for creating LArCaloHits
    lar_content::LArCaloHitFactory m_larCaloHitFactory;

    // Detector volume for voxelising the hits
    const GeometryManager *geom = pPrimaryPandora->GetGeometry();
    const LArTPC &tpc = geom->GetLArTPC();

    const float botX = tpc.GetCenterX() - 0.5 * tpc.GetWidthX();
    const float botY = tpc.GetCenterY() - 0.5 * tpc.GetWidthY();
    const float botZ = tpc.GetCenterZ() - 0.5 * tpc.GetWidthZ();
    const float topX = botX + tpc.GetWidthX();
    const float topY = botY + tpc.GetWidthY();
    const float topZ = botZ + tpc.GetWidthZ();
    const float voxelWidth(parameters.m_voxelWidth);
    const LArGrid grid(pandora::CartesianVector(botX, botY, botZ), pandora::CartesianVector(topX, topY, topZ),
        pandora::CartesianVector(voxelWidth, voxelWidth, voxelWidth));

    std::cout << "Total grid volume: bot = " << grid.m_bottom << "\n top = " << grid.m_top << std::endl;
    std::cout << "Making voxels with size " << grid.m_binWidths << std::endl;

    // Total number of entries in the TTree
    const int nEntries(ndsim->GetEntries());

    // Starting event
    const int startEvt = parameters.m_nEventsToSkip > 0 ? parameters.m_nEventsToSkip : 0;
    // Number of events to process, up to nEntries
    const int nProcess = parameters.m_nEventsToProcess > 0 ? parameters.m_nEventsToProcess : nEntries;
    // End event, up to nEntries
    const int endEvt = (startEvt + nProcess) < nEntries ? startEvt + nProcess : nEntries;

    std::cout << "Start event is " << startEvt << " and end event is " << endEvt - 1 << std::endl;

    for (int iEvt = startEvt; iEvt < endEvt; iEvt++)
    {
        if (parameters.m_shouldDisplayEventNumber)
            std::cout << std::endl << "   PROCESSING EVENT: " << iEvt << std::endl << std::endl;

        ndsim->GetEntry(iEvt);

        int hitCounter(0);

        // Create MCParticles from Geant4 trajectories
        MCParticleEnergyMap MCEnergyMap;
        for (size_t imcp = 0; imcp < larsed.mcp_id->size(); ++imcp)
        {
            MCEnergyMap[(*larsed.mcp_id)[imcp]] = (*larsed.mcp_energy)[imcp];
        }
        CreateSEDMCParticles(larsed, pPrimaryPandora, parameters);

        LArVoxelList voxelList;

        // Loop over the energy deposits and create voxels
        for (size_t ised = 0; ised < larsed.sed_det->size(); ++ised)
        {
            if ((*larsed.sed_det)[ised] == parameters.m_sensitiveDetName) // usually volTPCActive
            {
                const float startx = (*larsed.sed_startx)[ised];
                const float starty = (*larsed.sed_starty)[ised];
                const float startz = (*larsed.sed_startz)[ised];
                const float endx = (*larsed.sed_endx)[ised];
                const float endy = (*larsed.sed_endy)[ised];
                const float endz = (*larsed.sed_endz)[ised];
                const float energy = (*larsed.sed_energy)[ised];
                const int g4id = std::abs((*larsed.sed_id)[ised]);

                const pandora::CartesianVector start(startx, starty, startz);
                const pandora::CartesianVector end(endx, endy, endz);

                const LArHitInfo hitInfo(start, end, energy, g4id, parameters.m_lengthScale, parameters.m_energyScale);
                const LArVoxelList currentVoxelList = MakeVoxels(hitInfo, grid, parameters);

                for (const LArVoxel &voxel : currentVoxelList)
                    voxelList.emplace_back(voxel);
            }
        }

        std::cout << "Produced " << voxelList.size() << " voxels from " << larsed.sed_det->size() << " hit segments." << std::endl;

        // Merge voxels with the same IDs
        const LArVoxelList mergedVoxels = MergeSameVoxels(voxelList);

        std::cout << "Produced " << mergedVoxels.size() << " merged voxels from " << voxelList.size() << " voxels." << std::endl;

        // Stop processing the event if we have too many voxels: reco takes too long
        if (parameters.m_maxMergedVoxels > 0 && mergedVoxels.size() > parameters.m_maxMergedVoxels)
        {
            std::cout << "SKIPPING EVENT: number of merged voxels " << mergedVoxels.size() << " > " << parameters.m_maxMergedVoxels << std::endl;
            break;
        }

        // Loop over the voxels and make them into caloHits
        for (const LArVoxel &voxel : mergedVoxels)
        {
            const pandora::CartesianVector voxelPos(voxel.m_voxelPosVect);
            const float voxelE = voxel.m_energyInVoxel;
            const float MipE = 0.00075;
            const float voxelMipEquivalentE = voxelE / MipE;

            if (voxelMipEquivalentE > parameters.m_minVoxelMipEquivE)
            {
                lar_content::LArCaloHitParameters caloHitParameters;
                caloHitParameters.m_positionVector = voxelPos;
                caloHitParameters.m_expectedDirection = pandora::CartesianVector(0.f, 0.f, 1.f);
                caloHitParameters.m_cellNormalVector = pandora::CartesianVector(0.f, 0.f, 1.f);
                caloHitParameters.m_cellGeometry = pandora::RECTANGULAR;
                caloHitParameters.m_cellSize0 = voxelWidth;
                caloHitParameters.m_cellSize1 = voxelWidth;
                caloHitParameters.m_cellThickness = voxelWidth;
                caloHitParameters.m_nCellRadiationLengths = 1.f;
                caloHitParameters.m_nCellInteractionLengths = 1.f;
                caloHitParameters.m_time = 0.f;
                caloHitParameters.m_inputEnergy = voxelE;
                caloHitParameters.m_mipEquivalentEnergy = voxelMipEquivalentE;
                caloHitParameters.m_electromagneticEnergy = voxelE;
                caloHitParameters.m_hadronicEnergy = voxelE;
                caloHitParameters.m_isDigital = false;
                caloHitParameters.m_hitType = pandora::TPC_3D;
                caloHitParameters.m_hitRegion = pandora::SINGLE_REGION;
                caloHitParameters.m_layer = 0;
                caloHitParameters.m_isInOuterSamplingLayer = false;
                caloHitParameters.m_pParentAddress = (void *)(static_cast<uintptr_t>(++hitCounter));
                caloHitParameters.m_larTPCVolumeId = 0;
                caloHitParameters.m_daughterVolumeId = 0;

                if (parameters.m_use3D)
                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                        PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitParameters, m_larCaloHitFactory));

                if (parameters.m_useLArTPC)
                {
                    // Create U, V and W views assuming x is the common drift coordinate
                    const float x0_cm(voxelPos.GetX());
                    const float y0_cm(voxelPos.GetY());
                    const float z0_cm(voxelPos.GetZ());

                    lar_content::LArCaloHitParameters caloHitPars_UView(caloHitParameters);
                    caloHitPars_UView.m_hitType = pandora::TPC_VIEW_U;
                    const float upos_cm(pPrimaryPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoU(y0_cm, z0_cm));
                    caloHitPars_UView.m_positionVector = pandora::CartesianVector(x0_cm, 0.f, upos_cm);

                    lar_content::LArCaloHitParameters caloHitPars_VView(caloHitParameters);
                    caloHitPars_VView.m_hitType = pandora::TPC_VIEW_V;
                    const float vpos_cm(pPrimaryPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoV(y0_cm, z0_cm));
                    caloHitPars_VView.m_positionVector = pandora::CartesianVector(x0_cm, 0.f, vpos_cm);

                    lar_content::LArCaloHitParameters caloHitPars_WView(caloHitParameters);
                    caloHitPars_WView.m_hitType = pandora::TPC_VIEW_W;
                    const float wpos_cm(pPrimaryPandora->GetPlugins()->GetLArTransformationPlugin()->YZtoW(y0_cm, z0_cm));
                    caloHitPars_WView.m_positionVector = pandora::CartesianVector(x0_cm, 0.f, wpos_cm);

                    // Create LArCaloHits for U, V and W views
                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                        PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitPars_UView, m_larCaloHitFactory));
                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                        PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitPars_VView, m_larCaloHitFactory));
                    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                        PandoraApi::CaloHit::Create(*pPrimaryPandora, caloHitPars_WView, m_larCaloHitFactory));
                }

                // Set calo hit voxel to MCParticle relation using trackID
                const int trackID = voxel.m_trackID;

                // Find the energy fraction: voxelHitE/MCParticleE
                float energyFrac(0.f), MCEnergy(0.f);
                MCParticleEnergyMap::const_iterator mapIter = MCEnergyMap.find(trackID);
                if (mapIter != MCEnergyMap.end())
                    MCEnergy = mapIter->second;

                if (MCEnergy > 0.0)
                    energyFrac = voxelE / MCEnergy;

                PandoraApi::SetCaloHitToMCParticleRelationship(*pPrimaryPandora, (void *)((intptr_t)hitCounter), (void *)((intptr_t)trackID), energyFrac);
            }
        } // end voxel loop

    } // end segment detector loop

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::ProcessEvent(*pPrimaryPandora));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::Reset(*pPrimaryPandora));

    fileSource->Close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

MCParticleEnergyMap CreateEDepSimMCParticles(const TG4Event &event, const pandora::Pandora *const pPrimaryPandora, const Parameters &parameters)
{
    // Create map of trackID and energy
    MCParticleEnergyMap energyMap;

    if (!pPrimaryPandora)
    {
        std::cout << "Could not create MC particles, since pPrimaryPandora is null" << std::endl;
        return energyMap;
    }

    lar_content::LArMCParticleFactory mcParticleFactory;

    // Create the primary MC neutrino, linked to the trajectories below
    int neutrinoID(999999), neutrinoPDG(14), nuanceCode(1000);
    TLorentzVector neutrinoVtx, neutrinoP4;

    // Get the initial primary vertex
    if (event.Primaries.size() > 0)
    {
        const TG4PrimaryVertex &g4PrimaryVtx = event.Primaries[0];
        neutrinoVtx = g4PrimaryVtx.GetPosition() * parameters.m_lengthScale;
        std::cout << "Neutrino vertex = " << neutrinoVtx.X() << ", " << neutrinoVtx.Y() << ", " << neutrinoVtx.Z() << std::endl;

        const std::string reaction(g4PrimaryVtx.GetReaction());
        nuanceCode = GetNuanceCode(reaction);

        // Get the primary vertex particle information
        if (g4PrimaryVtx.Informational.size() > 0)
        {
            const TG4PrimaryVertex &g4Info = g4PrimaryVtx.Informational[0];

            // Get the first primary particle, which should be the neutrino.
            // Other primaries would be nuclei etc.
            if (g4Info.Particles.size() > 0)
            {
                const TG4PrimaryParticle &g4Primary = g4Info.Particles[0];

                neutrinoID = g4Primary.GetTrackId();
                neutrinoPDG = g4Primary.GetPDGCode();
                neutrinoP4 = g4Primary.GetMomentum() * parameters.m_energyScale;

                std::cout << "Neutrino ID = " << neutrinoID << ", PDG = " << neutrinoPDG << ", E = " << neutrinoP4.E()
                          << ", px = " << neutrinoP4.Px() << ", py = " << neutrinoP4.Py() << ", pz = " << neutrinoP4.Pz() << std::endl;
            }
        }
    }

    lar_content::LArMCParticleParameters mcNeutrinoParameters;
    mcNeutrinoParameters.m_nuanceCode = nuanceCode;
    mcNeutrinoParameters.m_process = lar_content::MC_PROC_INCIDENT_NU;
    mcNeutrinoParameters.m_energy = neutrinoP4.E();
    mcNeutrinoParameters.m_momentum = pandora::CartesianVector(neutrinoP4.Px(), neutrinoP4.Py(), neutrinoP4.Pz());
    mcNeutrinoParameters.m_vertex = pandora::CartesianVector(neutrinoVtx.X(), neutrinoVtx.Y(), neutrinoVtx.Z());
    mcNeutrinoParameters.m_endpoint = pandora::CartesianVector(neutrinoVtx.X(), neutrinoVtx.Y(), neutrinoVtx.Z());
    mcNeutrinoParameters.m_particleId = neutrinoPDG;
    mcNeutrinoParameters.m_mcParticleType = pandora::MC_3D;
    mcNeutrinoParameters.m_pParentAddress = (void *)((intptr_t)neutrinoID);
    PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPrimaryPandora, mcNeutrinoParameters, mcParticleFactory));

    std::cout << "Creating MC Particles from the Geant4 event trajectories" << std::endl;

    // Loop over trajectories
    for (const TG4Trajectory &g4Traj : event.Trajectories)
    {
        // LArMCParticle parameters
        lar_content::LArMCParticleParameters mcParticleParameters;

        // Initial momentum and energy in GeV (Geant4 uses MeV)
        const TLorentzVector initMtm(g4Traj.GetInitialMomentum() * parameters.m_energyScale);
        const float energy(initMtm.E());
        mcParticleParameters.m_energy = energy;
        mcParticleParameters.m_momentum = pandora::CartesianVector(initMtm.X(), initMtm.Y(), initMtm.Z());

        // Particle codes
        mcParticleParameters.m_particleId = g4Traj.GetPDGCode();
        mcParticleParameters.m_mcParticleType = pandora::MC_3D;
        mcParticleParameters.m_nuanceCode = mcNeutrinoParameters.m_nuanceCode.Get();

        // Set unique parent integer address using trackID
        const int trackID = g4Traj.GetTrackId();
        mcParticleParameters.m_pParentAddress = (void *)((intptr_t)trackID);

        // Start and end points in cm (Geant4 uses mm)
        const std::vector<TG4TrajectoryPoint> trajPoints = g4Traj.Points;
        const int nPoints(trajPoints.size());

        if (nPoints > 1)
        {
            const TG4TrajectoryPoint start = trajPoints[0];
            const TLorentzVector vertex = start.GetPosition() * parameters.m_lengthScale;
            mcParticleParameters.m_vertex = pandora::CartesianVector(vertex.X(), vertex.Y(), vertex.Z());

            const TG4TrajectoryPoint end = trajPoints[nPoints - 1];
            const TLorentzVector endPos = end.GetPosition() * parameters.m_lengthScale;
            mcParticleParameters.m_endpoint = pandora::CartesianVector(endPos.X(), endPos.Y(), endPos.Z());
            // Process ID
            mcParticleParameters.m_process = start.GetProcess();
        }
        else
        {
            // Should not reach here, but set sensible values just in case
            mcParticleParameters.m_vertex = pandora::CartesianVector(0.f, 0.f, 0.f);
            mcParticleParameters.m_endpoint = pandora::CartesianVector(0.f, 0.f, 0.f);
            mcParticleParameters.m_process = lar_content::MC_PROC_UNKNOWN;
        }

        // Create MCParticle
        PANDORA_THROW_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPrimaryPandora, mcParticleParameters, mcParticleFactory));

        // Set parent relationships
        const int parentID = g4Traj.GetParentId();

        if (parentID < 0) // link to mc neutrino
        {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                PandoraApi::SetMCParentDaughterRelationship(*pPrimaryPandora, (void *)((intptr_t)neutrinoID), (void *)((intptr_t)trackID)));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                PandoraApi::SetMCParentDaughterRelationship(*pPrimaryPandora, (void *)((intptr_t)parentID), (void *)((intptr_t)trackID)));
        }

        // Store particle energy for given trackID
        energyMap[trackID] = energy;
    }

    return energyMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CreateSEDMCParticles(const LArSED &larsed, const pandora::Pandora *const pPrimaryPandora, const Parameters &parameters)
{

    lar_content::LArMCParticleFactory mcParticleFactory;

    const int nuidoffset(100000000);

    // Create MC neutrinos
    for (size_t i = 0; i < larsed.nuPDG->size(); ++i)
    {
        const int neutrinoID = nuidoffset + i;
        const int neutrinoPDG = (*larsed.nuPDG)[i];
        const std::string reaction = GetNuanceReaction((*larsed.ccnc)[i], (*larsed.mode)[i]);
        const int nuanceCode = GetNuanceCode(reaction);

        const float nuVtxX = (*larsed.nuvtxx)[i] * parameters.m_lengthScale;
        const float nuVtxY = (*larsed.nuvtxy)[i] * parameters.m_lengthScale;
        const float nuVtxZ = (*larsed.nuvtxz)[i] * parameters.m_lengthScale;

        const float nuE = (*larsed.enu)[i] * parameters.m_energyScale;
        const float nuPx = nuE * (*larsed.nu_dcosx)[i];
        const float nuPy = nuE * (*larsed.nu_dcosy)[i];
        const float nuPz = nuE * (*larsed.nu_dcosz)[i];

        lar_content::LArMCParticleParameters mcNeutrinoParameters;
        mcNeutrinoParameters.m_nuanceCode = nuanceCode;
        mcNeutrinoParameters.m_process = lar_content::MC_PROC_INCIDENT_NU;

        mcNeutrinoParameters.m_energy = nuE;
        mcNeutrinoParameters.m_momentum = pandora::CartesianVector(nuPx, nuPy, nuPz);
        mcNeutrinoParameters.m_vertex = pandora::CartesianVector(nuVtxX, nuVtxY, nuVtxZ);
        mcNeutrinoParameters.m_endpoint = pandora::CartesianVector(nuVtxX, nuVtxY, nuVtxZ);

        mcNeutrinoParameters.m_particleId = neutrinoPDG;
        mcNeutrinoParameters.m_mcParticleType = pandora::MC_3D;
        mcNeutrinoParameters.m_pParentAddress = (void *)((intptr_t)neutrinoID);

        PANDORA_THROW_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPrimaryPandora, mcNeutrinoParameters, mcParticleFactory));
    }

    // Create MC particles
    for (size_t i = 0; i < larsed.mcp_id->size(); ++i)
    {
        // LArMCParticle parameters
        lar_content::LArMCParticleParameters mcParticleParameters;

        // Initial momentum and energy in GeV
        const float px = (*larsed.mcp_px)[i] * parameters.m_energyScale;
        const float py = (*larsed.mcp_py)[i] * parameters.m_energyScale;
        const float pz = (*larsed.mcp_pz)[i] * parameters.m_energyScale;
        const float energy = (*larsed.mcp_energy)[i] * parameters.m_energyScale;
        mcParticleParameters.m_energy = energy;
        mcParticleParameters.m_momentum = pandora::CartesianVector(px, py, pz);

        // Particle codes
        mcParticleParameters.m_particleId = (*larsed.mcp_pdg)[i];
        mcParticleParameters.m_mcParticleType = pandora::MC_3D;

        // Neutrino info
        const int nuid = (*larsed.mcp_nuid)[i];
        const int neutrinoID = nuid + nuidoffset;
        const std::string reaction = GetNuanceReaction((*larsed.ccnc)[nuid], (*larsed.mode)[nuid]);
        mcParticleParameters.m_nuanceCode = GetNuanceCode(reaction);

        // Set unique parent integer address using trackID
        const int trackID = (*larsed.mcp_id)[i];
        mcParticleParameters.m_pParentAddress = (void *)((intptr_t)trackID);

        // Start and end points in cm
        const float startx = (*larsed.mcp_startx)[i] * parameters.m_lengthScale;
        const float starty = (*larsed.mcp_starty)[i] * parameters.m_lengthScale;
        const float startz = (*larsed.mcp_startz)[i] * parameters.m_lengthScale;
        mcParticleParameters.m_vertex = pandora::CartesianVector(startx, starty, startz);

        const float endx = (*larsed.mcp_endx)[i] * parameters.m_lengthScale;
        const float endy = (*larsed.mcp_endy)[i] * parameters.m_lengthScale;
        const float endz = (*larsed.mcp_endz)[i] * parameters.m_lengthScale;
        mcParticleParameters.m_endpoint = pandora::CartesianVector(endx, endy, endz);

        // Process ID
        mcParticleParameters.m_process = lar_content::MC_PROC_UNKNOWN;

        // Create MCParticle
        PANDORA_THROW_RESULT_IF(
            pandora::STATUS_CODE_SUCCESS, !=, PandoraApi::MCParticle::Create(*pPrimaryPandora, mcParticleParameters, mcParticleFactory));

        // Set parent relationships
        const int parentID = (*larsed.mcp_mother)[i];

        if (parentID == 0) // link to mc neutrino
        {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                PandoraApi::SetMCParentDaughterRelationship(*pPrimaryPandora, (void *)((intptr_t)neutrinoID), (void *)((intptr_t)trackID)));
        }
        else
        {
            PANDORA_THROW_RESULT_IF(pandora::STATUS_CODE_SUCCESS, !=,
                PandoraApi::SetMCParentDaughterRelationship(*pPrimaryPandora, (void *)((intptr_t)parentID), (void *)((intptr_t)trackID)));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

int GetNuanceCode(const std::string &reaction)
{
    // The GENIE reaction string (also stored by edep-sim) is created using
    // https://github.com/GENIE-MC/Generator/blob/master/src/Framework/Interaction/Interaction.cxx#L249
    // String format is "nu:PDGId;tgt:PDGId;N:PDGId;proc:interactionType,scattering;", e.g.
    //                  "nu:14;tgt:1000180400;N:2112;proc:Weak[CC],QES;"

    // GENIE scattering codes:
    // https://github.com/GENIE-MC/Generator/blob/master/src/Framework/Interaction/ScatteringType.h

    // Nuance codes: https://internal.dunescience.org/doxygen/MCNeutrino_8h_source.html

    // GENIE conversion code for RooTracker output files:
    // https://github.com/GENIE-MC/Generator/blob/master/src/contrib/t2k/neut_code_from_rootracker.C
    // Similar code is available here (Neut reaction code):
    // https://internal.dunescience.org/doxygen/namespacegenie_1_1utils_1_1ghep.html

    // For now, just set the basic reaction types, excluding any specific final states:
    // https://github.com/GENIE-MC/Generator/blob/master/src/contrib/t2k/neut_code_from_rootracker.C#L276
    int code(1000);

    const bool is_cc = (reaction.find("Weak[CC]") != std::string::npos); // weak charged-current
    const bool is_nc = (reaction.find("Weak[NC]") != std::string::npos); // weak neutral-current
    // const bool is_charm = (reaction.find("charm")    != std::string::npos); // charm production
    const bool is_qel = (reaction.find("QES") != std::string::npos);   // quasi-elastic scattering
    const bool is_dis = (reaction.find("DIS") != std::string::npos);   // deep inelastic scattering
    const bool is_res = (reaction.find("RES") != std::string::npos);   // resonance
    const bool is_cohpi = (reaction.find("COH") != std::string::npos); // coherent pi
    const bool is_ve = (reaction.find("NuEEL") != std::string::npos);  // nu e elastic
    const bool is_imd = (reaction.find("IMD") != std::string::npos);   // inverse mu decay
    const bool is_mec = (reaction.find("MEC") != std::string::npos);   // meson exchange current

    if (is_qel)
    {
        code = 0;
        if (is_cc)
            code = 1001;
        else if (is_nc)
            code = 1002;
    }
    else if (is_dis)
    {
        code = 2;
        if (is_cc)
            code = 1091;
        else if (is_nc)
            code = 1092;
    }
    else if (is_res)
        code = 1;
    else if (is_cohpi)
    {
        code = 3;
        if (is_qel)
            code = 4;
    }
    else if (is_ve)
        code = 1098;
    else if (is_imd)
        code = 1099;
    else if (is_mec)
        code = 10;

    std::cout << "Reaction " << reaction << " has code = " << code << std::endl;

    return code;
}

std::string GetNuanceReaction(const int ccnc, const int mode)
{
    std::string reaction("");

    if (mode == 0)
    {
        reaction = "QES";
    }
    else if (mode == 2)
    {
        reaction = "DIS";
    }
    else if (mode == 1)
    {
        reaction = "RES";
    }
    else if (mode == 3)
    {
        reaction = "COH";
    }
    else if (mode == 5)
    {
        reaction = "NuEEL";
    }
    else if (mode == 6)
    {
        reaction = "IMD";
    }
    else if (mode == 10)
    {
        reaction = "MEC";
    }

    if (ccnc == 0)
    {
        reaction += "Weak[CC]";
    }
    else if (ccnc == 1)
    {
        reaction += "Weak[NC]";
    }

    return reaction;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArVoxelList MakeVoxels(const LArHitInfo &hitInfo, const LArGrid &grid, const Parameters &parameters)
{
    // Code based on https://github.com/chenel/larcv2/tree/edepsim-formattruth/larcv/app/Supera/Voxelize.cxx
    // which is made available under the MIT license (which is fully compatible with Pandora's GPLv3 license).

    LArVoxelList currentVoxelList;

    // Start and end positions
    const pandora::CartesianVector start(hitInfo.m_start);
    const pandora::CartesianVector stop(hitInfo.m_stop);

    // Direction vector and hit segment length
    const pandora::CartesianVector dir = stop - start;
    const float hitLength(dir.GetMagnitude());

    // Check hit length is greater than epsilon limit
    if (hitLength < std::numeric_limits<float>::epsilon())
        return currentVoxelList;

    // Hit segment total energy in GeV (Geant4 uses MeV)
    const float g4HitEnergy(hitInfo.m_energy);

    // Check hit energy is greater than epsilon limit
    if (g4HitEnergy < std::numeric_limits<float>::epsilon())
        return currentVoxelList;

    // Get the trackID of the (main) contributing particle.
    // ATTN: this can very rarely be more than one track
    const int trackID = hitInfo.m_trackID;

    // Define ray trajectory, which checks dirMag (hitLength) >= epsilon limit
    const pandora::CartesianVector dirNorm = dir.GetUnitVector();
    LArRay ray(start, dirNorm);

    // We need to shuffle along the hit segment path and create voxels as we go.
    // There are 4 cases for the start and end points inside the voxelisation region.
    // Case 1: start & stop are both inside the voxelisation boundary
    // Case 2: start & stop are both outside, but path direction intersects boundary
    // Case 3: start is inside boundary, stop = intersection at region boundary
    // Case 4: end is inside boundary, start = intersection at region boundary

    double t0(0.0), t1(0.0);
    pandora::CartesianVector point1(0.f, 0.f, 0.f), point2(0.f, 0.f, 0.f);

    // Check if the start and end points are inside the voxelisation region
    const bool inStart = grid.Inside(start);
    const bool inStop = grid.Inside(stop);

    if (inStart && inStop)
    {
        // Case 1: Start and end points are inside boundary
        point1 = start;
        point2 = stop;
    }
    else if (!inStart && !inStop)
    {
        // Case 2: Start and end points are outside boundary
        if (grid.Intersect(ray, t0, t1))
        {
            point1 = ray.GetPoint(t0);
            point2 = ray.GetPoint(t1);
        }
        else
            return currentVoxelList;
    }
    else if (inStart && !inStop)
    {
        // Case 3: Start inside boundary
        point1 = start;
        if (grid.Intersect(ray, t0, t1))
            point2 = ray.GetPoint(t1);
        else
            return currentVoxelList;
    }
    else if (!inStart && inStop)
    {
        // Case 4: End inside boundary
        point2 = stop;
        if (grid.Intersect(ray, t0, t1))
            point1 = ray.GetPoint(t0);
        else
            return currentVoxelList;
    }

    // Now create voxels between point1 and point2.
    // Ray direction will be the same, but update starting point
    ray.UpdateOrigin(point1);

    bool shuffle(true);

    // Keep track of total voxel path length so far
    float totalPath(0.0);
    int loop(0);

    while (shuffle)
    {
        // Get point along path to define voxel bin (bottom corner)
        const pandora::CartesianVector voxelPoint = ray.GetPoint(parameters.m_voxelPathShift);

        // Grid 3d bin containing this point; 4th element is the total bin number
        const LongBin4Array gridBins = grid.GetBinIndices(voxelPoint);
        const long voxelID = gridBins[3];
        const long xBin = gridBins[0];
        const long yBin = gridBins[1];
        const long zBin = gridBins[2];

        // Voxel bottom and top corners
        const pandora::CartesianVector voxBot = grid.GetPoint(xBin, yBin, zBin);
        const pandora::CartesianVector voxTop = grid.GetPoint(xBin + 1, yBin + 1, zBin + 1);

        // Voxel box
        const LArBox vBox(voxBot, voxTop);

        // Get ray intersections with this box: t0 and t1 are set as the start
        // and end intersection pathlengths relative to the current ray point.
        // If we can't find t0 and t1, then stop shuffling along the path
        if (!vBox.Intersect(ray, t0, t1))
            shuffle = false;

        // Voxel extent = intersection path difference
        double dL(t1 - t0);
        // For the first path length, use the distance from the
        // starting ray point to the 2nd intersection t1
        if (loop == 0)
            dL = t1;

        // Stop processing if we are not moving along the path
        if (dL < parameters.m_voxelPathShift)
            shuffle = false;

        totalPath += dL;

        // Stop adding voxels if we have enough
        if (totalPath > hitLength)
        {
            shuffle = false;
            // Adjust final path according to hit segment total length
            dL = hitLength - totalPath + dL;
        }

        // Voxel energy (GeV) using path length fraction w.r.t hit length.
        // Here, hitLength is guaranteed to be greater than zero
        const float voxelEnergy(g4HitEnergy * dL / hitLength);

        // Store voxel object in vector
        const LArVoxel voxel(voxelID, voxelEnergy, voxBot, trackID);
        currentVoxelList.emplace_back(voxel);

        // Update ray starting position using intersection path difference
        const pandora::CartesianVector newStart = ray.GetPoint(dL);
        ray.UpdateOrigin(newStart);
        loop++;
    }

    return currentVoxelList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArVoxelList MergeSameVoxels(const std::vector<LArVoxel> &voxelList)
{
    std::cout << "Merging voxels with the same IDs" << std::endl;
    LArVoxelList mergedVoxels;

    const int nVoxels = voxelList.size();
    std::vector<bool> processed(nVoxels, false);

    for (int i = 0; i < nVoxels; i++)
    {
        // Skip voxel if it was already used in a merge
        if (processed[i])
            continue;

        LArVoxel voxel1 = voxelList[i];
        float voxE1 = voxel1.m_energyInVoxel;
        const long id1 = voxel1.m_voxelID;
        const int trackid1 = voxel1.m_trackID;

        // Loop over other voxels (from i+1) and check if we have an ID match.
        // If so, add their energies and only store the combined voxel at the end
        for (int j = i + 1; j < nVoxels; j++)
        {
            // Skip voxel if it was already used in a merge
            if (processed[j])
                continue;

            const LArVoxel voxel2 = voxelList[j];
            const long id2 = voxel2.m_voxelID;
            const int trackid2 = voxel2.m_trackID;

            if (id2 == id1 && trackid1 == trackid2)
            {
                // IDs match. Add energy and set processed integer
                voxE1 += voxel2.m_energyInVoxel;
                processed[j] = true;
            }
        }

        // Add combined (or untouched) voxel to the merged list
        voxel1.SetEnergy(voxE1);
        mergedVoxels.emplace_back(voxel1);

        // We have processed the ith voxel
        processed[i] = true;
    }

    return mergedVoxels;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ParseCommandLine(int argc, char *argv[], Parameters &parameters)
{
    if (1 == argc)
        return PrintOptions();

    int cOpt(0);

    std::string recoOption("");
    std::string viewOption("3d");
    std::string formatOption("EDepSim");
    std::string geomFileName("");
    std::string inputTreeName("");

    while ((cOpt = getopt(argc, argv, "r:i:e:k:f:g:t:v:d:n:s:j:w:m:c:pNh")) != -1)
    {
        switch (cOpt)
        {
            case 'r':
                recoOption = optarg;
                break;
            case 'i':
                parameters.m_settingsFile = optarg;
                break;
            case 'e':
                parameters.m_inputFileName = optarg;
                break;
            case 'k':
                inputTreeName = optarg;
                break;
            case 'f':
                formatOption = optarg;
                break;
            case 'g':
                geomFileName = optarg;
                break;
            case 't':
                parameters.m_geomManagerName = optarg;
                break;
            case 'v':
                parameters.m_geometryVolName = optarg;
                break;
            case 'd':
                parameters.m_sensitiveDetName = optarg;
                break;
            case 'n':
                parameters.m_nEventsToProcess = atoi(optarg);
                break;
            case 's':
                parameters.m_nEventsToSkip = atoi(optarg);
                break;
            case 'p':
                parameters.m_printOverallRecoStatus = true;
                break;
            case 'j':
                viewOption = optarg;
                break;
            case 'w':
                parameters.m_voxelWidth = atof(optarg);
                break;
            case 'm':
                parameters.m_maxMergedVoxels = atoi(optarg);
                break;
            case 'c':
                parameters.m_minVoxelMipEquivE = atof(optarg);
                break;
            case 'N':
                parameters.m_shouldDisplayEventNumber = true;
                break;
            case 'h':
            default:
                return PrintOptions();
        }
    }

    ProcessViewOption(viewOption, parameters);
    ProcessFormatOption(formatOption, inputTreeName, geomFileName, parameters);
    return ProcessRecoOption(recoOption, parameters);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool PrintOptions()
{
    std::cout << std::endl
              << "./bin/PandoraInterface " << std::endl
              << "    -r RecoOption          (required) [Full, AllHitsCR, AllHitsNu, CRRemHitsSliceCR, CRRemHitsSliceNu, AllHitsSliceCR, AllHitsSliceNu]"
              << std::endl
              << "    -i Settings            (required) [algorithm description: xml]" << std::endl
              << "    -e EventsFile          (required) [input data ROOT file containing events (& geometry)]" << std::endl
              << "    -k EventsTree          (optional) [name of the input ROOT TTree (default = EDepSimEvents)]" << std::endl
              << "    -f DataFormat          (optional) [EDepSim (default, rooTracker format), SED (LArSoft-like)]" << std::endl
              << "    -g GeometryFile        (optional) [TGeoManager ROOT file (default = input EventsFile for EDepSim)]" << std::endl
              << "    -t TGeoManagerName     (optional) [TGeoManager name (default = EDepSimGeometry)]" << std::endl
              << "    -n NEventsToProcess    (optional) [no. of events to process]" << std::endl
              << "    -s NEventsToSkip       (optional) [no. of events to skip in "
                 "first file]"
              << std::endl
              << "    -p                     (optional) [print status]" << std::endl
              << "    -N                     (optional) [print event numbers]" << std::endl
              << "    -j Projection          (optional) [3D (default), LArTPC, Both]" << std::endl
              << "    -w width               (optional) [voxel bin width (cm), default = 0.4 cm]" << std::endl
              << "    -m maxMergedVoxels     (optional) [skip events that have N(merged voxels) > maxMergedVoxels (default = no events skipped)]"
              << std::endl
              << "    -c minMipEquivE        (optional) [minimum MIP equivalent energy, default = 0.3]" << std::endl
              << "    -v geometryVolName     (optional) [Geant4 geometry placement detector volume name, default = volArgonCubeDetector_PV_0]"
              << std::endl
              << "    -d sensitiveDetName    (optional) [Geant4 sensitive hits detector name, default = ArgonCube]" << std::endl
              << std::endl;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessViewOption(const std::string &viewOption, Parameters &parameters)
{
    std::string chosenViewOption(viewOption);
    std::transform(chosenViewOption.begin(), chosenViewOption.end(), chosenViewOption.begin(), ::tolower);

    if (chosenViewOption == "3d")
    {
        // 3D hits only
        std::cout << "Using 3D hits" << std::endl;
        parameters.m_useLArTPC = false;
        parameters.m_use3D = true;
    }
    else if (chosenViewOption == "both")
    {
        // LArTPC and 3D hits
        std::cout << "Using LArTPC projections _and_ 3D hits" << std::endl;
        parameters.m_useLArTPC = true;
        parameters.m_use3D = true;
    }
    else
    {
        // LArTPC hits only
        std::cout << "Using LArTPC projections" << std::endl;
        parameters.m_useLArTPC = true;
        parameters.m_use3D = false;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProcessRecoOption(const std::string &recoOption, Parameters &parameters)
{
    std::string chosenRecoOption(recoOption);
    std::transform(chosenRecoOption.begin(), chosenRecoOption.end(), chosenRecoOption.begin(), ::tolower);

    if ("full" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = true;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = true;
    }
    else if ("allhitscr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = false;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("nostitchingcr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = false;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("allhitsnu" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = false;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("crremhitsslicecr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = true;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("crremhitsslicenu" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = true;
        parameters.m_shouldRunStitching = true;
        parameters.m_shouldRunCosmicHitRemoval = true;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("allhitsslicecr" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = false;
        parameters.m_shouldRunCosmicRecoOption = true;
        parameters.m_shouldPerformSliceId = false;
    }
    else if ("allhitsslicenu" == chosenRecoOption)
    {
        parameters.m_shouldRunAllHitsCosmicReco = false;
        parameters.m_shouldRunStitching = false;
        parameters.m_shouldRunCosmicHitRemoval = false;
        parameters.m_shouldRunSlicing = true;
        parameters.m_shouldRunNeutrinoRecoOption = true;
        parameters.m_shouldRunCosmicRecoOption = false;
        parameters.m_shouldPerformSliceId = false;
    }
    else
    {
        std::cout << "LArReco, Unrecognized reconstruction option: " << recoOption << std::endl << std::endl;
        return PrintOptions();
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessFormatOption(const std::string &formatOption, const std::string &inputTreeName, const std::string &geomFileName, Parameters &parameters)
{
    std::string chosenFormatOption(formatOption);
    std::transform(chosenFormatOption.begin(), chosenFormatOption.end(), chosenFormatOption.begin(), ::tolower);

    if (chosenFormatOption == "sed")
    {
        // LArSoft-type SimEnergyDeposit (SED) ROOT format
        parameters.m_dataFormat = Parameters::LArNDFormat::SED;
        // Set the geometry file name
        parameters.m_geomFileName = geomFileName;
        // All lengths are already in cm, so don't rescale
        parameters.m_lengthScale = 1.0f;
        // All energies are already in GeV, so don't rescale
        parameters.m_energyScale = 1.0f;
        // Set expected input TTree name for SED data
        parameters.m_inputTreeName = "simdump/ndsim";
    }
    else
    {
        // Assume EDepSim rooTracker format
        parameters.m_dataFormat = Parameters::LArNDFormat::EDepSim;
        // TGeoManager is stored in the input rooTracker file containing the hits
        parameters.m_geomFileName = parameters.m_inputFileName;
        // All lengths are in mm, so we need to convert them to cm
        parameters.m_lengthScale = parameters.m_mm2cm;
        // All energies are in MeV, so we need to convert them to GeV
        parameters.m_energyScale = parameters.m_MeV2GeV;
    }

    // Default input tree name is "EDepSimEvents"
    if (inputTreeName.size() > 0)
        parameters.m_inputTreeName = inputTreeName;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ProcessExternalParameters(const Parameters &parameters, const Pandora *const pPandora)
{
    auto *const pEventSteeringParameters = new lar_content::MasterAlgorithm::ExternalSteeringParameters;
    pEventSteeringParameters->m_shouldRunAllHitsCosmicReco = parameters.m_shouldRunAllHitsCosmicReco;
    pEventSteeringParameters->m_shouldRunStitching = parameters.m_shouldRunStitching;
    pEventSteeringParameters->m_shouldRunCosmicHitRemoval = parameters.m_shouldRunCosmicHitRemoval;
    pEventSteeringParameters->m_shouldRunSlicing = parameters.m_shouldRunSlicing;
    pEventSteeringParameters->m_shouldRunNeutrinoRecoOption = parameters.m_shouldRunNeutrinoRecoOption;
    pEventSteeringParameters->m_shouldRunCosmicRecoOption = parameters.m_shouldRunCosmicRecoOption;
    pEventSteeringParameters->m_shouldPerformSliceId = parameters.m_shouldPerformSliceId;
    pEventSteeringParameters->m_printOverallRecoStatus = parameters.m_printOverallRecoStatus;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraApi::SetExternalParameters(*pPandora, "LArMaster", pEventSteeringParameters));
}

} // namespace lar_nd_reco
