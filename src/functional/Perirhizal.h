// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#ifndef PERIRHIZAL_H_
#define PERIRHIZAL_H_

#include "MappedOrganism.h"

#include <vector>

namespace CPlantBox {

/**
 * Wraps a MappedSegments (or specialisations MappedPlant, MappedRootsystem)
 * and adds functions to help modelling the perirhizal zone.
 *
 * Currently rather useless, main part is in Python (in Perirhizal.py),
 * but perfomance critical methods could be implemented in C++
 */
class Perirhizal
{
public:

    Perirhizal() { }
    Perirhizal(std::shared_ptr<MappedSegments> ms) :ms(ms) { }

	std::map<int,double> sumSegFluxes(const std::vector<double>& segFluxes); ///< sums segment fluxes over soil cells,  soilFluxes = sumSegFluxes(segFluxes), [cm3/day]
    std::vector<double> splitSoilFluxes(const std::vector<double>& soilFluxes, int type = 0) const; ///< splits soil fluxes (per cell) into segment fluxes
    std::vector<double> segOuterRadii(int type, const std::vector<double>& vols = std::vector<double>(0)) const; ///< outer cylinder radii to match cell volume

    std::shared_ptr<MappedSegments> ms;
    std::vector<double> adapt_values(std::vector<double> val_new_, 
                     double minVal_, double maxVal_, 
                     const std::vector<double>& volumes_, 
                     bool divideEqually_, bool verbose_);///<  Update val_new_ to remain between bounds while maintaining as much as possible the gradient
                     
    std::vector<double> distributeValSolute_(
        std::vector<double> seg_values_content, 
        const std::vector<double>& volumes, 
        double source, 
        double dt) ;
        
    std::vector<double> distributeValWater_(
        std::vector<double> seg_values_perVol, 
        const std::vector<double>& volumes, 
        double source, 
        double dt, 
        double theta_S, 
        double theta_wilting_point) ;
        
    std::vector<double> splitSoilVals(
        const std::vector<double>& soilVals, 
        const std::vector<int>& cellIds, 
        bool isWater, 
        const std::vector<double>& seg_values, 
        const std::vector<double>& seg_volume, 
        double dt, 
        double theta_S,
        double theta_wilting_point) ;
        

protected:
    void redistribute_excess();
    void redistribute_deficit();
    void handle_excess();
    void handle_deficit() ;
    
    std::vector<double> val_new; 
    double minVal; 
    double maxVal;
    std::vector<double> volumes;
    bool divideEqually;    
    bool verbose;
};

}

#endif
