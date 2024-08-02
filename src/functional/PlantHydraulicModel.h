// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef PLANT_HYDRAULIC_MODEL_
#define PLANT_HYDRAULIC_MODEL_

#include "MappedOrganism.h"
#include "PlantHydraulicParameters.h"

namespace CPlantBox {

/**
 * Hydraulic model using parameters stored in PlantHydraulicParameters
 *
 * Hybrid solver (Meunier et al. 2017)
 * Doussan solver (Doussan et al. 2006)
 * *
 */
class PlantHydraulicModel
{
public:

    PlantHydraulicModel(std::shared_ptr<CPlantBox::MappedSegments> ms, std::shared_ptr<CPlantBox::PlantHydraulicParameters> params);

    virtual ~PlantHydraulicModel() { }

    void linearSystemMeunier(double simTime, const std::vector<double> sx, bool cells = true); ///< builds linear system (simTime is needed for age dependent conductivities)

    std::vector<double> getRadialFluxes(double simTime, const std::vector<double> rx, const std::vector<double> sx, bool approx = false, bool cells = false) const; // for each segment in [cm3/day]
    std::map<int,double> sumSegFluxes(const std::vector<double> segFluxes); ///< sums segment fluxes over soil cells,  soilFluxes = sumSegFluxes(segFluxes), [cm3/day]

    std::shared_ptr<CPlantBox::MappedSegments> ms;
    std::shared_ptr<CPlantBox::PlantHydraulicParameters> params;

    std::vector<int> aI; // to assemble the sparse matrix on the Python side
    std::vector<int> aJ;
    std::vector<double> aV;
    std::vector<double> aB;

protected:

    virtual size_t fillVectors(size_t k, int i, int j, double bi, double cii, double cij, double psi_s) ; ///< fills row k of Meunier matrix
	virtual double getPsiOut(bool cells, int si, const std::vector<double>& sx_) const; ///< get the outer water potential [cm]

};

} // namespace

#endif
