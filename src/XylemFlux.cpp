// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "XylemFlux.h"

#include <algorithm>
#include <set>

namespace CPlantBox {



XylemFlux::XylemFlux(std::shared_ptr<CPlantBox::MappedSegments> rs): rs(rs)
    {
    size_t length_leaf = std::count(rs->organTypes.begin(), rs->organTypes.end(), 4);
    gs.resize(length_leaf);
    pg.resize(length_leaf);
    }


/**
 * Assembles the linear system as sparse matrix, given by public member variables,
 * indices aI, aJ, and corresponding values aV; and load aB
 *
 * @param simTime[day]  	current simulation time, needed for age dependent conductivities,
 *                  		to calculate the age from the creation times (age = sim_time - segment creation time).
 * @param sx [cm]			soil matric potential in the cells or around the segments, given per cell or per segment
 * @param cells 			sx per cell (true), or segments (false)
 * @param soil_k [day-1]    optionally, soil conductivities can be prescribed per segment,
 *                          conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)
 */
void XylemFlux::linearSystem(double simTime, const std::vector<double>& sx, bool cells, const std::vector<double> soil_k)
{
    int Ns = rs->segments.size(); // number of segments
    aI.resize(4*Ns);
    aJ.resize(4*Ns);
    aV.resize(4*Ns);
    int N = rs->nodes.size(); // number of nodes
    aB.resize(N);
    std::fill(aB.begin(), aB.end(), 0.);
    std::fill(aV.begin(), aV.end(), 0.);
    std::fill(aI.begin(), aI.end(), 0);
    std::fill(aJ.begin(), aJ.end(), 0);
    size_t k=0;
    size_t numleaf = 0;
    for (int si = 0; si<Ns; si++) {

        int i = rs->segments[si].x;
        int j = rs->segments[si].y;

        double psi_s;
        int organType = rs->organTypes[si];
        if (cells) { // soil matric potential given per cell
            int cellIndex = rs->seg2cell[j-1];
            if (cellIndex>=0) {
                if(sx.size()>1) {
                    psi_s = sx.at(cellIndex);
                } else {
                    psi_s = sx.at(0);
                }
            } else {
                psi_s = airPressure;
            }
        } else {
            psi_s = sx.at(j-1); // j-1 = segIdx = s.y-1
        }
        if (organType == 4 && pg.at(0)!= 0){
            psi_s = pg.at(numleaf);
        }
        double a = rs->radii[si]; // si is correct, with ordered and unordered segmetns
        double age = simTime - rs->nodeCTs[j];
        int subType = rs->subTypes[si];
        double kx = 0.;
        double  kr = 0.;

        try {
            kx = kx_f(si, age, subType, organType);
            kr = kr_f(si, age, subType, organType, numleaf);
        } catch(...) {
            std::cout << "\n XylemFlux::linearSystem: conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }

        if(organType == 4) {
            numleaf +=1;
        }
        //        if (age<=0) {
        //            std::cout << si << ", " << j <<" age leq 0 " << age << ", " << kx <<  ", " << kr << ", time "<< simTime << ", " << rs->nodeCTs[j] << "\n";
        //        }
        if (soil_k.size()>0) {
            kr = std::min(kr, soil_k[si]);
        }

        auto n1 = rs->nodes[i];
        auto n2 = rs->nodes[j];
        auto v = n2.minus(n1);
        double l = v.length();
        if (l<1.e-5) {
            std::cout << "XylemFlux::linearSystem: warning segment length smaller 1.e-5 \n"; // quick fix?
            l = 1.e-5;
        }
        double vz = v.z / l; // normed direction

        double tau = std::sqrt(2.*a * M_PI * kr / kx); // Eqn (2)
        double delta = std::exp(-tau * l) - std::exp(tau * l); // Eqn (5)
        double idelta = 1. / delta;

        double cii = -kx * idelta * tau * (std::exp(-tau * l) + std::exp(tau * l)); // Eqn (16)
        double cij = 2 * kx * idelta * tau;  // Eqn 17
        double bi = kx * vz; //  # Eqn 18

        aB[i] += ( bi + cii * psi_s +cij * psi_s) ;
        aI[k] = i; aJ[k]= i; aV[k] = cii;
        k += 1;
        aI[k] = i; aJ[k] = j;  aV[k] = cij;
        k += 1;

        int ii = i;
        i = j;  j = ii; // edge ji
        aB[i] += ( -bi + cii * psi_s +cij * psi_s) ; // (-bi) Eqn (14) with changed sign
        aI[k] = i; aJ[k]= i; aV[k] = cii;
        k += 1;
        aI[k] = i; aJ[k] = j;  aV[k] = cij;
        k += 1;
    }
}

/**
 * Assembles the linear system as sparse matrix, given by public member variables,
 * indices aI, aJ, and corresponding values aV; and load aB
 *
 * @param simTime[day]      current simulation time, needed for age dependent conductivities,
 *                          to calculate the age from the creation times (age = sim_time - segment creation time).
 * @param sx [cm]           soil matric potential in the cells or around the segments, given per cell or per segment
 * @param cells             sx per cell (true), or segments (false)
 * @param soil_k [day-1]    optionally, soil conductivities can be prescribed per segment,
 *                          conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)
 */
void XylemFlux::linearSystem_detached(double simTime, const std::vector<double>& sx, bool cells, const std::vector<double> soil_k)
{
    int Ns = rs->segments.size(); // number of segments
    aI.resize(4*Ns);
    aJ.resize(4*Ns);
    aV.resize(4*Ns);
    int N = rs->nodes.size(); // number of nodes
    aB.resize(N);
    std::fill(aB.begin(), aB.end(), 0.);
    std::fill(aV.begin(), aV.end(), 0.);
    std::fill(aI.begin(), aI.end(), 0);
    std::fill(aJ.begin(), aJ.end(), 0);
    size_t k=0;
    size_t numleaf = 0;
    for (int si = 0; si<Ns; si++) {

        int i = rs->segments[si].x;
        int j = rs->segments[si].y;

        double psi_s;
        int organType = rs->organTypes[si];
        if (cells) { // soil matric potential given per cell
            int cellIndex = rs->seg2cell[j-1];
            if (cellIndex>=0) {
                if(sx.size()>1) {
                    psi_s = sx.at(cellIndex);
                } else {
                    psi_s = sx.at(0);
                }
            } else {
                psi_s = airPressure;
            }
        } else {
            psi_s = sx.at(j-1); // j-1 = segIdx = s.y-1
        }
        if (organType == 4 && pg.at(0)!= 0){
            psi_s = pg.at(numleaf);
        }
        double a = rs->radii[si]; // si is correct, with ordered and unordered segmetns
        double age = simTime - rs->nodeCTs[j];
        int subType = rs->subTypes[si];
        double kx = 0.;
        double  kr = 0.;

        try {
            kx = kx_f(si, age, subType, organType);
            kr = kr_f(si, age, subType, organType, numleaf);
        } catch(...) {
            std::cout << "\n XylemFlux::linearSystem: conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }

        if(organType == 4) {
            numleaf +=1;
        }
        //        if (age<=0) {
        //            std::cout << si << ", " << j <<" age leq 0 " << age << ", " << kx <<  ", " << kr << ", time "<< simTime << ", " << rs->nodeCTs[j] << "\n";
        //        }
        if (soil_k.size()>0) {
            kr = std::min(kr, soil_k[si]);
        }

        auto n1 = rs->nodes[i];
        auto n2 = rs->nodes[j];
        auto v = n2.minus(n1);
        double l = v.length();
        if (l<1.e-5) {
            std::cout << "XylemFlux::linearSystem: warning segment length smaller 1.e-5 \n"; // quick fix?
            l = 1.e-5;
        }
        double vz = v.z / l; // normed direction

        double tau = std::sqrt(2.*a * M_PI * kr / kx); // Eqn (2)
        double delta = std::exp(-tau * l) - std::exp(tau * l); // Eqn (5)
        double idelta = 1. / delta;

        double cii = -kx * idelta * tau * (std::exp(-tau * l) + std::exp(tau * l)); // Eqn (16)
        double cij = 2 * kx * idelta * tau;  // Eqn 17
        double bi = kx * vz; //  # Eqn 18

        i = rs->segments[si].y; // <---------------------------------- that' it (everything gets connected to the first segment)

        aB[i] += ( bi + cii * psi_s +cij * psi_s) ;
        aI[k] = i; aJ[k]= i; aV[k] = cii;
        k += 1;
        aI[k] = i; aJ[k] = j;  aV[k] = cij;
        k += 1;

        int ii = i;
        i = j;  j = ii; // edge ji
        aB[i] += ( -bi + cii * psi_s +cij * psi_s) ; // (-bi) Eqn (14) with changed sign
        aI[k] = i; aJ[k]= i; aV[k] = cii;
        k += 1;
        aI[k] = i; aJ[k] = j;  aV[k] = cij;
        k += 1;
    }
}


/**
 * Fluxes from root segments into soil cells
 *
 * @param simTime   [days] current simulation time is needed for age dependent conductivities,
 *                  to calculate the age from the creation times (age = sim_time - segment creation time).
 * @param rx        [cm] root xylem matric potential
 * @param sx        [cm] soil matric potential for each cell
 * @param approx    approximate or exact (default = false, i.e. exact)
 *
 * @return hash map with cell indices as keys and fluxes as values [cm3/day]
 */
std::map<int,double> XylemFlux::soilFluxes(double simTime, const std::vector<double>& rx, const std::vector<double>& sx,
    bool approx, const std::vector<double> soil_k)
{
    return sumSegFluxes(segFluxes(simTime,  rx, sx, approx, true, soil_k));
}

/**
 * Volumetric fluxes for each segment according to a given solution @param rx and @param sx
 *
 * @param simTime   [days] current simulation time is needed for age dependent conductivities,
 *                  to calculate the age from the creation times (age = sim_time - segment creation time).
 * @param rx        [cm] root xylem matric potential
 * @param sx        [cm] soil matric potential for each segment
 * @param approx    approximate or exact (default = false, i.e. exact)
 * @param cells     sx per cell (true), or segments (false)
 * @param soil_k    [day-1] optionally, soil conductivities can be prescribed per segment,
 *                          conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)
 *
 * @return Volumetric fluxes for each segment [cm3/day]
 */
std::vector<double> XylemFlux::segFluxes(double simTime, const std::vector<double>& rx, const std::vector<double>& sx,
    bool approx, bool cells, const std::vector<double> soil_k)
{
    std::vector<double> fluxes = std::vector<double>(rs->segments.size());
    size_t numleaf = 0;
    for (int si = 0; si<rs->segments.size(); si++) {

        int i = rs->segments[si].x;
        int j = rs->segments[si].y;
        int organType = rs->organTypes[si];

        double psi_s;
        if (cells) { // soil matric potential given per cell
            int cellIndex = rs->seg2cell[j-1];
            if (cellIndex>=0) {
                if(sx.size()>1) {
                    psi_s = sx.at(cellIndex);
                } else {
                    psi_s = sx.at(0);
                }
            } else {
                psi_s = airPressure;
            }
        } else {
            psi_s = sx.at(j-1); // j-1 = segIdx = s.y-1
        }

        if (organType == 4 && pg.at(0)!= 0){
            psi_s = pg.at(numleaf);
        }

        double a = rs->radii[si]; // si is correct, with ordered and unordered segments
        double age = simTime - rs->nodeCTs[j];
        int subType = rs->subTypes[si];

        double kx = 0.;
        double kr = 0.;
        try {
            kx = kx_f(si, age, subType, organType);
            kr = kr_f(si, age, subType, organType, numleaf);
        } catch(...) {
            std::cout << "\n XylemFlux::segFluxes: conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }
        if (soil_k.size()>0) {
            kr = std::min(kr, soil_k[si]);
        }

        if (organType == 4) {
            numleaf +=1;
        }

        Vector3d n1 = rs->nodes[i];
        Vector3d n2 = rs->nodes[j];
        double l = (n2.minus(n1)).length();

        double f = -2*a*M_PI*kr; // flux is proportional to f // *rho*g
        double fApprox = f*l*(psi_s - rx[j]); // cm3 / day

        double tau = std::sqrt(2*a*M_PI*kr/kx); // sqrt(c) [cm-1]
        double d = std::exp(-tau*l)-std::exp(tau*l); // det
        double fExact = -f*(1./(tau*d))*(rx[i]-psi_s+rx[j]-psi_s)*(2.-std::exp(-tau*l)-std::exp(tau*l));

        double flux = fExact*(!approx)+approx*fApprox;
        fluxes[si] = flux;

    }
    return fluxes;
}

/**
 * Sums segment fluxes over each cell
 *
 * @param segFluxes 	segment fluxes given per segment index [cm3/day]
 * @return hash map with cell indices as keys and fluxes as values [cm3/day]
 */
std::map<int,double> XylemFlux::sumSegFluxes(const std::vector<double>& segFluxes)
{
    std::map<int,double> fluxes;
    for (int si = 0; si<rs->segments.size(); si++) {
        int j = rs->segments[si].y;
        int segIdx = j-1;
        if (rs->seg2cell.count(segIdx)>0) {
            int cellIdx = rs->seg2cell[segIdx];
            if (cellIdx>=0) {
                if (fluxes.count(cellIdx)==0) {
                    fluxes[cellIdx] = segFluxes[segIdx];
                } else {
                    fluxes[cellIdx] = fluxes[cellIdx] + segFluxes[segIdx]; // sum up fluxes per cell
                }
            }
        }
    }
    return fluxes;
}

/**
 * Splits soil fluxes per cell to the segments within the cell, so that the summed fluxes agree, @see sumSoilFluxes()
 *
 * @param soilFluxes 	cell fluxes per global index [cm3/day]
 * @param type 			split flux proportional to 0: segment volume, 1: segment surface, 2: segment length
 * @return fluxes for each segment [cm3/day]
 */
std::vector<double> XylemFlux::splitSoilFluxes(const std::vector<double>& soilFluxes, int type) const
{
    auto lengths =  this->segLength();
    std::vector<double> fluxes = std::vector<double>(rs->segments.size());
    std::fill(fluxes.begin(), fluxes.end(), 0.);
    auto map = rs->cell2seg;
    for(auto iter = map.begin(); iter != map.end(); ++iter) {
        int cellId =  iter->first;
        auto segs = map.at(cellId);
        double v = 0.;  // calculate sum over cell
        for (int i : segs) {
            if (type==0) { // volume
                v += M_PI*(rs->radii[i]*rs->radii[i])*lengths[i];
            } else if (type==1) { // surface
                v += 2*M_PI*rs->radii[i]*lengths[i];
            } else if (type==2) { // length
                v += lengths[i];
            }
        }
        for (int i : segs) { // calculate outer radius
            double t =0.; // proportionality factor (must sum up to == 1 over cell)
            if (type==0) { // volume
                t = M_PI*(rs->radii[i]*rs->radii[i])*lengths[i]/v;
            } else if (type==1) { // surface
                t = 2*M_PI*rs->radii[i]*lengths[i]/v;
            } else if (type==2) { // length
                t = lengths[i]/v;
            }
            fluxes[i] = t*soilFluxes.at(cellId);
        }
    }
    return fluxes;
}

/**
 *  Sets the radial conductivity in [1 day-1]
 * TODO: make deprecated: in the examples, replace setKr[Kr] by setKr[[Kr]]
 */
//either age or type/subtype dependent
void XylemFlux::setKr(std::vector<double> values, std::vector<double> age) {
    kr =  { values }; //because kr is std::vector<std::vector<double>>
    kr_t = { age };
    if (age.size()==0) {
        if (values.size()==1) {
            kr_f = std::bind(&XylemFlux::kr_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
            std::cout << "Kr is constant " << values[0] << " cm2 day g-1 \n";
        } else {
            kr_f  = std::bind(&XylemFlux::kr_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
            std::cout << "Kr is constant per type, type 0 = " << values[0] << " cm2 day g-1 \n";
        }
    } else {
        kr_f  = std::bind(&XylemFlux::kr_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
        std::cout << "Kr is age dependent\n";
    }
}

/**
 *  Sets the radial conductivity in [1 day-1]
 * in case of organ_type specific kr
 */
//either age or type/subtype dependent
void XylemFlux::setKr(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age) {
    kr = values;
    kr_t = age;
    if (age.size()==0) {
        if (values.size()==1) {
            if (values[0].size()==1) {
                kr_f = std::bind(&XylemFlux::kr_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
                std::cout << "Kr is constant " << values[0][0] << " cm2 day g-1 \n";
            } else {
                kr_f  = std::bind(&XylemFlux::kr_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
                std::cout << "Kr is constant per subtype, subtype 0 = " << values[0][0] << " cm2 day g-1 \n";
            }
        } else {
            if (values[0].size()==1) {
                kr_f = std::bind(&XylemFlux::kr_perOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
                std::cout << "Kr is constant per organ type, organ type 2 (root) = " << values[0][0] << " cm2 day g-1 \n";
            } else {
                kr_f  = std::bind(&XylemFlux::kr_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
                std::cout << "Kr is constant per subtype of organ type, for root, subtype 0 = " << values[0][0] << " cm2 day g-1 \n";
            }
        }
    } else {
        kr_f  = std::bind(&XylemFlux::kr_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
        std::cout << "Kr is equal for all organs and age dependent\n";
    }
}

/**
 *  Sets the axial conductivity in [cm3 day-1]
 * TODO: make deprecated: in the examples, replace setKx[Kx] by setKx[[Kx]]
 */
//either age or type/subtype dependent
void XylemFlux::setKx(std::vector<double> values, std::vector<double> age) {
    kx = { values }; // because kx is std::vector<std::vector<double>>
    kx_t = { age };
    if (age.size()==0) {
        if (values.size()==1) {
            kx_f = std::bind(&XylemFlux::kx_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
            std::cout << "Kx is constant " << values[0] << " cm2 day g-1 \n";
        } else {
            kx_f  = std::bind(&XylemFlux::kx_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
            std::cout << "Kx is constant per subtype, subtype 0 = " << values[0] << " cm2 day g-1 \n";
        }
    } else {
        kx_f  = std::bind(&XylemFlux::kx_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        std::cout << "Kx is age dependent\n";
    }
}

//either age or type/subtype dependent
void XylemFlux::setKx(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age) {
    kx = values;
    kx_t = age;
    if (age.size()==0) {
        if (values.size()==1) {
            if (values[0].size()==1) {
                kx_f = std::bind(&XylemFlux::kx_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                std::cout << "Kx is constant " << values[0][0] << " cm2 day g-1 \n";
            } else {
                kx_f  = std::bind(&XylemFlux::kx_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                std::cout << "Kx is constant per subtype, subtype 0 = " << values[0][0] << " cm2 day g-1 \n";
            }
        } else {
            if (values[0].size()==1) {
                kx_f = std::bind(&XylemFlux::kx_perOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                std::cout << "Kx is constant per organ type, organ type 2 (root) = " << values[0][0] << " cm2 day g-1 \n";
            } else {
                kx_f  = std::bind(&XylemFlux::kx_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                std::cout << "Kx is constant per subtype of organ type, for root, subtype 0 = " << values[0][0] << " cm2 day g-1 \n";
            }
        }
    } else {
        kx_f  = std::bind(&XylemFlux::kx_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        std::cout << "Kx is equal for all organs and age dependent\n";
    }
}

/**
 * Sets the radial conductivity in [1 day-1]
 *	@param values 			kr values for age (its linearly interpolated between these values) for each root type
 *	@param age 				ages for the given values for each root type
 *
 * TODO: make deprecated (i would leave it in for now, used in pyhton_modules/root_conductivities.py)
 */
//both age and type/subtype dependent
void XylemFlux::setKrTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age) {
    krs= { values };
    krs_t = { age };
    kr_f = std::bind(&XylemFlux::kr_tablePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
    std::cout << "Kr is age dependent per root type\n";
}

/**
 *  Sets the axial conductivity in [cm3 day-1]
 *	@param values 			kx values for age (its linearly interpolated between these values) for each root type
 *	@param age 				ages for the given values for each root type
 *
 * TODO: make deprecated (i would leave it in for now, used in pyhton_modules/root_conductivities.py)
 */
//both age and type/subtype dependent
void XylemFlux::setKxTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age) {
    kxs = {values};
    kxs_t = {age};
    kx_f = std::bind(&XylemFlux::kx_tablePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    std::cout << "Kx is age dependent per root type\n";
}

/**
 *  Sets the radial conductivity in [1 day-1] age, organ type and sub-type dependent
 *
 *	@param values 			kr values for age (its linearly interpolated between these values) for each organ type and each sub-type
 *	@param age 				ages for the given values for each organ type and for each sub type
 */
void XylemFlux::setKrTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age) {
    krs = values;
    krs_t = age;
    if (age[0].size()==1) {
        kr_f = std::bind(&XylemFlux::kr_tablePerOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
        std::cout << "Kr is age dependent per organ type\n";
    }
    else{
        kr_f = std::bind(&XylemFlux::kr_tablePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);
        std::cout << "Kr is age dependent per organ type and sub type\n";
    }
}

/**
 *  Sets the axial conductivity in [cm3 day-1] age, organ type and sub-type dependent
 *
 *	@param values 			kx values for age (its linearly interpolated between these values) for each organ type and each sub type
 *	@param age 				ages for the given values for each organ type and for each sub-type
 */
void XylemFlux::setKxTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age) {
    kxs= values;
    kxs_t = age;
    if (age[0].size()==1) {
        kx_f = std::bind(&XylemFlux::kx_tablePerOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        std::cout << "Kx is age dependent per organ type\n";
    }
    else {
        kx_f = std::bind(&XylemFlux::kx_tablePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        std::cout << "Kx is age dependent per organ type and sub type\n";
    }
}

/**
 *
 */
void XylemFlux::setKrValues(std::vector<double> values) {
    assert(values.size() == rs->segments.size() && "XylemFlux::setKrValues: values size must equal number of segments");
    kr.clear();
    kr_t.clear();
    kr.push_back(values);
}

/**
 *
 */
void XylemFlux::setKxValues(std::vector<double> values) {
    assert(values.size() == rs->segments.size() && "XylemFlux::setKrValues: values size must equal number of segments");
    kx.clear();
    kx_t.clear();
    kx.push_back(values);
}

} // namespace
