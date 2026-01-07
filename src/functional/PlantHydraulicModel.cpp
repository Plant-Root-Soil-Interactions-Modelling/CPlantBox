// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "PlantHydraulicModel.h"

#include <algorithm>
#include <set>

namespace CPlantBox {

PlantHydraulicModel::PlantHydraulicModel(std::shared_ptr<CPlantBox::MappedSegments> ms, std::shared_ptr<CPlantBox::PlantHydraulicParameters> params):
    ms(ms), params(params) { params->ms = ms;}


/**
 * Assembles the linear system as sparse matrix, given by public member variables,
 * indices aI, aJ, and corresponding values aV; and load aB
 *
 * @param simTime[day]  	current simulation time, needed for age dependent conductivities,
 *                  		to calculate the age from the creation times (age = sim_time - segment creation time).
 * @param sx [cm]			soil matric potential in the cells or around the segments, given per cell or per segment
 * @param cells 			sx per cell (true), or segments (false)
 * @param soil_k    [day-1] optionally, soil conductivities can be prescribed per segment,
 *                          conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)
 */
void PlantHydraulicModel::linearSystemMeunier(double simTime, const std::vector<double> sx, bool cells,
			const std::vector<double> soil_k)
{
    int Ns = ms->segments.size(); // number of segments
    aI.resize(4*Ns);
    aJ.resize(4*Ns);
    aV.resize(4*Ns);
    int N = ms->nodes.size(); // number of nodes
    aB.resize(N);
    std::fill(aB.begin(), aB.end(), 0.);
    std::fill(aV.begin(), aV.end(), 0.);
    std::fill(aI.begin(), aI.end(), 0);
    std::fill(aJ.begin(), aJ.end(), 0);
    size_t k=0;

    for (int si = 0; si<Ns; si++) {

        int i = ms->segments[si].x;
        int j = ms->segments[si].y;

        double psi_s = getPsiOut(cells, si, sx);
        double age = simTime - ms->nodeCTs[j];
        int organType = ms->organTypes[si];
        int subType = ms->subTypes[si];
        double kx = 0.;
        double kr = 0.;

        try {
            kx = params->kx_f(si, age, subType, organType);
            kr = params->kr_f(si, age, subType, organType);
        } catch(...) {
            std::cout << "\n XylemFlux::linearSystem: conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }
        if (soil_k.size()>0) {
            kr = std::min(kr, soil_k[si]);
        }

        auto n1 = ms->nodes[i];
        auto n2 = ms->nodes[j];
        auto v = n2.minus(n1);
        double l = v.length();
        if (l<1.e-5) {
            // std::cout << "XylemFlux::linearSystem: warning segment length smaller 1.e-5 \n";
            l = 1.e-5; // valid quick fix? (also in segFluxes)
        }
		double perimeter = ms->getPerimeter(si, l);//perimeter of exchange surface
        double vz = v.z / l; // normed direction

        double cii, cij, bi;

        if (perimeter * kr>1.e-16) {
            double tau = std::sqrt(perimeter * kr / kx); // Eqn (6)
            double delta = std::exp(-tau * l) - std::exp(tau * l); // Eqn (12)
            double idelta = 1. / delta;
            cii = -kx * idelta * tau * (std::exp(-tau * l) + std::exp(tau * l)); // Eqn (23)
            cij = 2 * kx * idelta * tau;  // Eqn 24
            bi = kx * vz; //  # Eqn 25
        } else { // solution for a=0, or kr = 0
            cii = kx/l;
            cij = -kx/l;
            bi = kx * vz;
            psi_s = 0;//
        }

		k = fillVectors(k, i, j, bi, cii, cij, psi_s);
    }
}

/**
 * Volumetric radial fluxes for each segment according to a given solution @param rx and @param sx
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
std::vector<double> PlantHydraulicModel::getRadialFluxes(double simTime, const std::vector<double> rx, const std::vector<double> sx,
    bool approx, bool cells, const std::vector<double> soil_k) const
{
    std::vector<double> fluxes = std::vector<double>(ms->segments.size());
    for (int si = 0; si<ms->segments.size(); si++) {

        int i = ms->segments.at(si).x;
        int j = ms->segments.at(si).y;
        int organType = ms->organTypes.at(si);

        double psi_s = getPsiOut(cells, si, sx);


        // si is correct, with ordered and unordered segments
        double age = simTime - ms->nodeCTs[j];
        int subType = ms->subTypes[si];

        double kx = 0.;
        double kr = 0.;
        try {
            kx = params->kx_f(si, age, subType, organType);
            kr = params->kr_f(si, age, subType, organType);
        } catch(...) {
            std::cout << "\n XylemFlux::segFluxes: conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }
        if (soil_k.size()>0) {
            kr = std::min(kr, soil_k[si]);
        }

        auto n1 = ms->nodes.at(i);
        auto n2 = ms->nodes.at(j);
        auto v = n2.minus(n1);
        double l = v.length();
        if (l<1.e-5) {
            // std::cout << "XylemFlux::linearSystem: warning segment length smaller 1.e-5 \n";
            l = 1.e-5; // valid quick fix? (also in segFluxes)
        }

		double perimeter = ms->getPerimeter(si, l);//perimeter of exchange surface

        if (perimeter * kr>1.e-16) { // only relevant for exact solution
            double f = -perimeter*kr; // flux is proportional to f // *rho*g
            if (approx) {
                fluxes[si] = f*l*(psi_s - rx[j]); // cm3 / day
            } else {
                double tau = std::sqrt(perimeter*kr/kx); // sqrt(c) [cm-1]
                double d = std::exp(-tau*l)-std::exp(tau*l); // det
                double fExact = -f*(1./(tau*d))*(rx[i]-psi_s+rx[j]-psi_s)*(2.-std::exp(-tau*l)-std::exp(tau*l));
                if (!std::isfinite(fExact)) {
                    std::cout << "XylemFlux::segFluxes: nan or Inf fExact. segIdx "<<si<<" organType "<<organType<<" subType "<<subType;
                    std::cout <<" tau " << tau << ", l " << l << ", d "<<" perimeter "<<perimeter<<" kr "<<kr;
                    std::cout<< d << ", rx "<< rx[i] << ", psi_s " << psi_s << ", f " << f << "\n";
                    throw std::runtime_error("XylemFlux::segFluxes: nan or Inf fExact");
                }
                fluxes[si] = fExact;
            }
        } else {
            fluxes[si] = 0.;
        }

    }
    return fluxes;
}

/**
 * Sums segment fluxes over each cell
 *
 * @param segFluxes 	segment fluxes given per segment index [cm3/day]
 * @return hash map with cell indices as keys and fluxes as values [cm3/day]
 */
std::map<int,double> PlantHydraulicModel::sumSegFluxes(const std::vector<double> segFluxes)
{
    std::map<int,double> fluxes;
    for (int si = 0; si<ms->segments.size(); si++) {
        int j = ms->segments[si].y;
        int segIdx = j-1;
        if (ms->seg2cell.count(segIdx)>0) {
            int cellIdx = ms->seg2cell[segIdx];
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
 * fill the matrices to be solved. Overloaded by @see Photosynthesis::fillVectors
 * @param k				index for the row- and column-index vectors
 * @param i, j			indexes of the non-zero elements of the sparse matrix
 * @param psi_s 		outer water potential [cm]
 * @param bi			value of variable b at row i [cm3/d]
 * @param cii			value of variable c at row i col i [cm2/d]
 * @param cij			value of variable c at row i col j [cm2/d]
 * @return k			next index for the row- and column-index vectors
 */

size_t PlantHydraulicModel::fillVectors(size_t k, int i, int j, double bi, double cii, double cij, double psi_s)
{
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
	return k;
}

/**
 *  give outer water potential [cm] overloaded by @see Photosynthesis::getPsiOut
 * @param cells 		sx per cell (true), or segments (false)
 * @param si 			segment index
 * @param sx        [cm] soil matric potential for each cell
 */

double PlantHydraulicModel::getPsiOut(bool cells, int si, const std::vector<double>& sx_) const
{
	int organType = ms->organTypes.at(si);
    double psi_s;
	if (cells) { // soil matric potential given per cell
		int cellIndex = ms->seg2cell.at(si);
		if (cellIndex>=0) {
			if(organType ==Organism::ot_leaf){ //add a runtime error?
				std::cout<<"XylemFlux::linearSystem: Leaf segment n#"<<si<<" below ground. OrganType: ";
				std::cout<<organType<<" cell Index: "<<cellIndex<<std::endl;
			}
			if(sx_.size()>1) {
				psi_s = sx_.at(cellIndex);
			} else {
				psi_s = sx_.at(0);
			}
		} else {
			if(organType == Organism::ot_root) //add a runtime error?
			{
				std::cout<<"XylemFlux::linearSystem: Root segment n#"<<si<<" aboveground. OrganType: ";
				std::cout<<organType<<" cell Index: "<<cellIndex<<std::endl;
			}
			psi_s = params->psi_air;
		}
	} else {
		psi_s = sx_.at(si); // j-1 = segIdx = s.y-1
	}
	return psi_s;
}


} // namespace
