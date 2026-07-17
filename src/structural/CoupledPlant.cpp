#include "CoupledPlant.h"

#include <stdexcept>
#include <string>
#include <vector>

namespace CPlantBox {

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

CoupledPlant::CoupledPlant(const std::string& participantName,
                           const std::string& configFile,
                           int rank, int size)
    : Plant()
    , participant_(std::make_unique<precice::Participant>(participantName, configFile, rank, size))
{
}

// ---------------------------------------------------------------------------
// Lifecycle
// ---------------------------------------------------------------------------

void CoupledPlant::initialize(bool verbose, std::string mode) {
    // 1. Initialise the plant structure (roots, stems, leaves).
    Plant::initialize(verbose, mode);

    // 2. Register the initial segment midpoints as preCICE coupling vertices.
    resetCouplingMesh();

    // 3. Start the preCICE coupling session.
    participant_->initialize();

    // Allocate exchange buffers to match the registered vertex count.
    const std::size_t nSeg = vertexIDs_.size();
    soilPressures_.assign(nSeg, 0.0);
    waterFluxes_.assign(nSeg, 0.0);
}

void CoupledPlant::simulate(double dt, bool verbose) {
    if (!participant_) {
        throw std::runtime_error(
            "CoupledPlant::simulate: participant not initialised — "
            "call initialize() first.");
    }

    // 1. Read soil matric pressure heads from the coupled solver.
    participant_->readData(
        meshName_, pressureDataName_,
        vertexIDs_,
        dt,                 // relativeReadTime: end-of-window values
        soilPressures_);

    // 2. Advance the plant growth simulation.
    Plant::simulate(dt, verbose);

    // 3. Compute water-uptake fluxes using the freshly read soil pressures.
    waterFluxes_ = computeWaterFluxes(soilPressures_);

    if (waterFluxes_.size() != vertexIDs_.size()) {
        throw std::runtime_error(
            "CoupledPlant::simulate: computeWaterFluxes() returned " +
            std::to_string(waterFluxes_.size()) + " values but " +
            std::to_string(vertexIDs_.size()) + " were expected.");
    }

    // 4. Write water-uptake fluxes back to the coupled solver.
    participant_->writeData(
        meshName_, fluxDataName_,
        vertexIDs_,
        waterFluxes_);

    // 5. Advance the preCICE coupling step.
    participant_->advance(dt);
}

void CoupledPlant::finalize() {
    if (participant_) {
        participant_->finalize();
    }
}

// ---------------------------------------------------------------------------
// Coupling state queries
// ---------------------------------------------------------------------------

bool CoupledPlant::isCouplingOngoing() const {
    return participant_ && participant_->isCouplingOngoing();
}

double CoupledPlant::getMaxTimeStepSize() const {
    if (!participant_) {
        throw std::runtime_error(
            "CoupledPlant::getMaxTimeStepSize: participant not initialised.");
    }
    return participant_->getMaxTimeStepSize();
}

// ---------------------------------------------------------------------------
// Protected overridables
// ---------------------------------------------------------------------------

std::vector<double> CoupledPlant::computeWaterFluxes(
    const std::vector<double>& soilPressures)
{
    // Default: zero flux everywhere.  Override to attach a hydraulic model.
    return std::vector<double>(soilPressures.size(), 0.0);
}

void CoupledPlant::resetCouplingMesh() {
    const int meshDim = participant_->getMeshDimensions(meshName_);
    if (meshDim != 3) {
        throw std::runtime_error(
            "CoupledPlant: preCICE mesh \"" + meshName_ + "\" has dimension " +
            std::to_string(meshDim) + " but 3 is required.");
    }

    const std::vector<double> coords = segmentMidpointCoords_();
    const std::size_t nSeg = coords.size() / static_cast<std::size_t>(meshDim);
    vertexIDs_.resize(nSeg);

    participant_->setMeshVertices(meshName_, coords, vertexIDs_);
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

std::vector<double> CoupledPlant::segmentMidpointCoords_() const {
    const auto nodes    = getNodes();
    const auto segments = getSegments();

    std::vector<double> coords;
    coords.reserve(segments.size() * 3);

    for (const auto& seg : segments) {
        const Vector3d& a = nodes.at(static_cast<std::size_t>(seg.x));
        const Vector3d& b = nodes.at(static_cast<std::size_t>(seg.y));
        coords.push_back(0.5 * (a.x + b.x));
        coords.push_back(0.5 * (a.y + b.y));
        coords.push_back(0.5 * (a.z + b.z));
    }
    return coords;
}

} // namespace CPlantBox
