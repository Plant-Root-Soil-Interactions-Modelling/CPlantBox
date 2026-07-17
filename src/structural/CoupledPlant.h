#ifndef COUPLED_PLANT_H_
#define COUPLED_PLANT_H_

#include "Plant.h"

#include <precice/precice.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace CPlantBox {

/**
 * CoupledPlant
 *
 * Adaptor that integrates a CPlantBox Plant simulation into a preCICE
 * multi-physics coupling (https://precice.org).
 *
 * The coupling mesh is built from the midpoints of all plant segments.
 * Per coupling step the adaptor:
 *   - reads soil matric pressure head [cm] from the coupled solver, and
 *   - writes water-uptake flux [cm³/day] per segment back to it.
 *
 * Water flux computation is delegated to the virtual method
 * computeWaterFluxes().  The default returns zero for all segments; derive
 * from CoupledPlant and override that method to attach a hydraulic model.
 *
 * Typical usage
 * -------------
 *   CoupledPlant plant("PlantSolver", "precice-config.xml");
 *   plant.readParameters("myplant.xml");
 *   plant.initialize();                    // sets up Plant AND preCICE mesh
 *   while (plant.isCouplingOngoing()) {
 *       double dt = plant.getMaxTimeStepSize();
 *       plant.simulate(dt);                // read → grow → write → advance
 *   }
 *   plant.finalize();
 *
 * CMake
 * -----
 * Link against the preCICE library in CMakeLists.txt:
 *   find_package(precice 3 REQUIRED CONFIG)
 *   target_link_libraries(<target> PRIVATE precice::precice)
 *
 * Notes
 * -----
 * - The coupling mesh is registered once during initialize().  Segments that
 *   grow after that call are not represented in preCICE.  Override simulate()
 *   and call resetCouplingMesh() each step if a fully dynamic mesh is needed
 *   (requires preCICE to be configured for mesh reset between coupling steps).
 * - Requires preCICE >= 3.0.
 */
class CoupledPlant : public Plant {
public:
    /**
     * @param participantName  Name of this participant as defined in the
     *                         preCICE configuration file.
     * @param configFile       Path to the preCICE XML configuration file.
     * @param rank             MPI rank of this process (default 0).
     * @param size             Total number of MPI ranks (default 1).
     */
    CoupledPlant(const std::string& participantName,
                 const std::string& configFile,
                 int rank = 0, int size = 1);

    virtual ~CoupledPlant() = default;

    // -----------------------------------------------------------------------
    // Configuration — call before initialize()
    // -----------------------------------------------------------------------

    /** Name of the preCICE mesh that carries the plant coupling vertices.
     *  Default: "PlantMesh". Must match the precice-config.xml. */
    void setMeshName(const std::string& name) { meshName_ = name; }

    /** preCICE data name for the soil matric pressure head [cm] that is
     *  read from the coupled solver.  Default: "SoilPressure". */
    void setSoilPressureDataName(const std::string& name) { pressureDataName_ = name; }

    /** preCICE data name for the water-uptake flux [cm³/day] that is written
     *  to the coupled solver.  Default: "WaterFlux". */
    void setWaterFluxDataName(const std::string& name) { fluxDataName_ = name; }

    // -----------------------------------------------------------------------
    // Lifecycle
    // -----------------------------------------------------------------------

    /** Initialises the Plant (calls Plant::initialize), builds the preCICE
     *  coupling mesh from the initial segment midpoints, and calls
     *  precice::Participant::initialize(). */
    void initialize(bool verbose = true, std::string mode = "") override;

    /** Executes one coupled time step:
     *  1. Reads soil pressure heads from preCICE.
     *  2. Advances the Plant simulation by dt.
     *  3. Computes water fluxes via computeWaterFluxes().
     *  4. Writes water fluxes to preCICE.
     *  5. Calls precice::Participant::advance(dt).
     *
     *  @param dt  Time step size [days].  Should be at most
     *             getMaxTimeStepSize() to respect the coupling scheme. */
    void simulate(double dt, bool verbose = false) override;

    /** Calls precice::Participant::finalize() and releases resources.
     *  Must be called once after the coupling loop exits. */
    void finalize();

    // -----------------------------------------------------------------------
    // Coupling state queries
    // -----------------------------------------------------------------------

    /** Returns true while the preCICE coupling has not yet completed. */
    bool isCouplingOngoing() const;

    /** Returns the maximum time step size allowed by the current preCICE
     *  coupling scheme. */
    double getMaxTimeStepSize() const;

    // -----------------------------------------------------------------------
    // Access to exchanged data (read-only)
    // -----------------------------------------------------------------------

    /** Soil matric pressure heads [cm] per segment, last read from preCICE.
     *  The i-th entry corresponds to the i-th segment returned by
     *  Organism::getSegments(). */
    const std::vector<double>& getSoilPressures() const { return soilPressures_; }

    /** Water-uptake fluxes [cm³/day] per segment, last written to preCICE. */
    const std::vector<double>& getWaterFluxes() const { return waterFluxes_; }

protected:
    /**
     * Compute water-uptake fluxes [cm³/day] per segment given soil matric
     * pressure heads [cm].
     *
     * Override this in a derived class to connect a hydraulic model (e.g.
     * based on MappedPlant / XylemFlux).  The default implementation returns
     * zero flux for every segment.
     *
     * @param soilPressures  Matric pressure head [cm] per segment.
     * @return               Water-uptake flux [cm³/day] per segment.
     *                       Must have the same length as @p soilPressures.
     */
    virtual std::vector<double> computeWaterFluxes(const std::vector<double>& soilPressures);

    /** (Re-)registers segment midpoints with preCICE.  May be called before
     *  precice::Participant::initialize() to set up or refresh the mesh. */
    void resetCouplingMesh();

private:
    std::unique_ptr<precice::Participant> participant_;

    std::string meshName_        = "PlantMesh";
    std::string pressureDataName_ = "SoilPressure";
    std::string fluxDataName_    = "WaterFlux";

    std::vector<precice::VertexID> vertexIDs_;
    std::vector<double>            soilPressures_;
    std::vector<double>            waterFluxes_;

    /// Returns a flat array [x0,y0,z0, x1,y1,z1, ...] of segment midpoints.
    std::vector<double> segmentMidpointCoords_() const;
};

} // namespace CPlantBox

#endif /* COUPLED_PLANT_H_ */
