// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "XylemFlux.h"

#include <algorithm>
#include <set>

namespace CPlantBox {



XylemFlux::XylemFlux(std::shared_ptr<CPlantBox::MappedSegments> rs): rs(rs){}


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
void XylemFlux::linearSystem(double simTime, const std::vector<double>& sx, bool cells, const std::vector<double> soil_k,
								bool verbose)
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

    for (int si = 0; si<Ns; si++) {

        int i = rs->segments[si].x;
        int j = rs->segments[si].y;

        double psi_s = getPsiOut(cells, si, sx, verbose);
        double age = simTime - rs->nodeCTs[j];
        int organType = rs->organTypes[si];
        int subType = rs->subTypes[si];
        double kx = 0.;
        double  kr = 0.;

        try {
            kx = kx_f(si, age, subType, organType);
            kr = kr_f_wrapped(si, age, subType, organType, cells);
        } catch(...) {
            std::cout << "\n XylemFlux::linearSystem: conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }
        if (soil_k.size()>0) {
            kr = std::min(kr, soil_k[si]);
        }

        auto n1 = rs->nodes[i];
        auto n2 = rs->nodes[j];
        auto v = n2.minus(n1);
        double l = v.length();
        if (l<1.e-5) {
            // std::cout << "XylemFlux::linearSystem: warning segment length smaller 1.e-5 \n";
            l = 1.e-5; // valid quick fix? (also in segFluxes)
        }
		double perimeter = rs->getPerimeter(si, l);//perimeter of exchange surface
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
            psi_s = 0;
        }

		k = fillVectors(k, i, j, bi, cii, cij, psi_s);
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
    return rs->sumSegFluxes(segFluxes(simTime,  rx, sx, approx, true, soil_k));
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
    bool approx, bool cells, const std::vector<double> soil_k, bool verbose) const
{
    std::vector<double> fluxes = std::vector<double>(rs->segments.size());
    for (int si = 0; si<rs->segments.size(); si++) {

        int i = rs->segments.at(si).x;
        int j = rs->segments.at(si).y;
        int organType = rs->organTypes.at(si);

        double psi_s = getPsiOut(cells, si, sx, verbose);


        // si is correct, with ordered and unordered segments
        double age = simTime - rs->nodeCTs[j];
        int subType = rs->subTypes[si];

        double kx = 0.;
        double kr = 0.;
        try {
            kx = kx_f(si, age, subType, organType);
            kr = kr_f_wrapped(si, age, subType, organType, cells);
        } catch(...) {
            std::cout << "\n XylemFlux::segFluxes: conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }
        if (soil_k.size()>0) {
            kr = std::min(kr, soil_k[si]);
        }
        auto n1 = rs->nodes.at(i);
        auto n2 = rs->nodes.at(j);
        auto v = n2.minus(n1);
        double l = v.length();
        if (l<1.e-5) {
            // std::cout << "XylemFlux::linearSystem: warning segment length smaller 1.e-5 \n";
            l = 1.e-5; // valid quick fix? (also in segFluxes)
        }

		double perimeter = rs->getPerimeter(si, l);//perimeter of exchange surface


        if (perimeter * kr>1.e-16) { // only relevant for exact solution
            double f = -perimeter*kr; // flux is proportional to f // *rho*g
            double fApprox = f*l*(psi_s - rx[j]); // cm3 / day

            double tau = std::sqrt(perimeter*kr/kx); // sqrt(c) [cm-1]
            double d = std::exp(-tau*l)-std::exp(tau*l); // det
            double fExact = -f*(1./(tau*d))*(rx[i]-psi_s+rx[j]-psi_s)*(2.-std::exp(-tau*l)-std::exp(tau*l));
            if(!std::isfinite(fExact)) {
            	std::cout << "XylemFlux::segFluxes: nan or Inf fExact. segIdx "<<si<<" organType "<<organType<<" subType "<<subType;
				std::cout <<" tau " << tau << ", l " << l << ", d "<<" perimeter "<<perimeter<<" kr "<<kr;
				std::cout<< d << ", rx "<< rx[i] << ", psi_s " << psi_s << ", f " << f << "\n";
				throw std::runtime_error("XylemFlux::segFluxes: nan or Inf fExact");
			}
            double flux = fExact*(!approx)+approx*fApprox;
            fluxes[si] = flux;
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
std::map<int,double> XylemFlux::sumSegFluxes(const std::vector<double>& segFluxes)
{
    std::cout << "XylemFlux::sumSegFluxes: DEPRICATED use class MappedSegments instead (ms.sumSegFluxes) \n" << std::flush;
    std::map<int,double> fluxes;
    for (int si = 0; si<rs->segments.size(); si++) {
        int j = rs->segments[si].y;
        int segIdx = j-1;

        if (rs->seg2cell.count(segIdx)>0)
		{
			int cellIdx = rs->seg2cell[segIdx];
			if (cellIdx>=0)
			{
				if(rs->organTypes[segIdx] == Organism::ot_root)//only divid the fluxes between the root segments
				{
					if (fluxes.count(cellIdx)==0) {
						fluxes[cellIdx] = segFluxes[segIdx];
					} else {
						fluxes[cellIdx] = fluxes[cellIdx] + segFluxes[segIdx]; // sum up fluxes per cell
					}
				}else{
					if(segFluxes[segIdx] != 0.)
					{
						std::stringstream errMsg;
						errMsg<<"XylemFlux::sumSegFluxes. ot:"<<rs->organTypes[segIdx]<<" segIdx:"<<segIdx
						<<" cellIdx:"<<cellIdx<<" segFluxes[segIdx] :"
						<<segFluxes[segIdx]<<"=> shoot segment bellow ground ans exchanges water" <<std::endl;

						throw std::runtime_error(errMsg.str().c_str());
					}
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
    std::cout << "XylemFlux::splitSoilFluxes: DEPRICATED use class MappedSegments instead (ms.splitSoilFluxes) \n" << std::flush;
    auto lengths =  this->rs->segLength();
    std::vector<double> fluxes = std::vector<double>(rs->segments.size());
    std::fill(fluxes.begin(), fluxes.end(), 0.);
    auto map = rs->cell2seg;
	double fluxesTotTot =0;
    for(auto iter = map.begin(); iter != map.end(); ++iter) {
        int cellId =  iter->first;
        auto segs = map.at(cellId);
		if (cellId>=0) {
			double v = 0.;  // calculate sum over cell
			for (int i : segs) {

				if(rs->organTypes[i] == Organism::ot_root)//only divid the fluxes between the root segments
				{
					if (type==0) { // volume
						v += M_PI*(rs->radii[i]*rs->radii[i])*lengths[i];
					} else if (type==1) { // surface
						v += 2*M_PI*rs->radii[i]*lengths[i];
					} else if (type==2) { // length
						v += lengths[i];
					}
				}
			}
			double fluxesTot = 0;
			for (int i : segs) { // calculate outer radius

				if(rs->organTypes[i] == Organism::ot_root)
				{
					double t =0.; // proportionality factor (must sum up to == 1 over cell)
					if (type==0) { // volume
						t = M_PI*(rs->radii[i]*rs->radii[i])*lengths[i]/v;
					} else if (type==1) { // surface
						t = 2*M_PI*rs->radii[i]*lengths[i]/v;
					} else if (type==2) { // length
						t = lengths[i]/v;
					}
					if(fluxes[i] !=0){std::cout<<"fluxes "<<i<<" already set "<<std::endl;assert(false);}
					fluxes[i] = t*soilFluxes.at(cellId);
					fluxesTot +=  t*soilFluxes.at(cellId);
					fluxesTotTot += t*soilFluxes.at(cellId);
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

size_t XylemFlux::fillVectors(size_t k, int i, int j, double bi, double cii, double cij, double psi_s)
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

double XylemFlux::getPsiOut(bool cells, int si, const std::vector<double>& sx_, bool verbose) const
{
	int organType = rs->organTypes.at(si);
    double psi_s;
	if (cells) { // soil matric potential given per cell
		int cellIndex = rs->seg2cell.at(si);
		if (cellIndex>=0) {
			if((organType ==Organism::ot_leaf) && verbose){ //add a runtime error?
				std::cout<<"XylemFlux::linearSystem: Leaf segment n#"<<si<<" below ground. OrganType: ";
				std::cout<<organType<<" cell Index: "<<cellIndex<<std::endl;
			}
			if(sx_.size()>1) {
				psi_s = sx_.at(cellIndex);
			} else {
				psi_s = sx_.at(0);
			}
		} else {
			if((organType == Organism::ot_root) && verbose) //add a runtime error?
			{
				std::cout<<"XylemFlux::linearSystem: Root segment n#"<<si<<" aboveground. OrganType: ";
				std::cout<<organType<<" cell Index: "<<cellIndex<<std::endl;
			}
			psi_s = psi_air;
		}
	} else {
		psi_s = sx_.at(si); // j-1 = segIdx = s.y-1
	}
	return psi_s;
}



/**
 *  Sets the radial conductivity in [1 day-1]
 * TODO: make deprecated: in the examples, replace setKr[Kr] by setKr[[Kr]]
 */
//either age or type/subtype dependent
void XylemFlux::setKr(std::vector<double> values, std::vector<double> age, bool verbose) {
    kr =  { values }; //because kr is std::vector<std::vector<double>>
    kr_t = { age };
    if (age.size()==0) {
        if (values.size()==1) {
            kr_f = std::bind(&XylemFlux::kr_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
            if(verbose)
			{
				std::cout << "Kr is constant " << values[0] << " 1 day-1 \n";
            }
        } else {
            kr_f  = std::bind(&XylemFlux::kr_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
            if(verbose)
			{
				std::cout << "Kr is constant per type, type 0 = " << values[0] << " 1 day-1 \n";
            }
        }
    } else {
        kr_f  = std::bind(&XylemFlux::kr_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if(verbose)
        {
            std::cout << "Kr is age dependent\n";
        }
    }
}

/**
 *  Sets the radial conductivity in [1 day-1]
 * in case of organ_type specific kr
 * @param values 		kr per organ pr/and organ type) or/and per age [cm-1]
 * @param age 			ages if kr per age
 * @param kr_length_ 	exchange zone in root, where kr > 0 [cm from root tip], default = -1.0, i.e., no kr_length
 */
//either age or type/subtype dependent
void XylemFlux::setKr(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, double kr_length_, bool verbose) {
    kr = values;
    kr_t = age;
    if (age.size()==0) {
        if (values.size()==1) {
            if (values[0].size()==1) {
                kr_f = std::bind(&XylemFlux::kr_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose)
				{
					std::cout << "Kr is constant " << values[0][0] << " 1 day-1 \n";
                }
            } else {
                kr_f  = std::bind(&XylemFlux::kr_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose)
				{
					std::cout << "Kr is constant per subtype, subtype 0 = " << values[0][0] << " 1 day-1 \n";
                }
            }
        } else {
            if (values[0].size()==1) {
                kr_f = std::bind(&XylemFlux::kr_perOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose)
				{
					std::cout << "Kr is constant per organ type, organ type 2 (root) = " << values[0][0] << " 1 day-1 \n";
                }
            } else {
				if(kr_length_ > 0.){
					if(verbose)
					{
						std::cout << "Exchange zone in roots: kr > 0 until "<< kr_length_<<"cm from root tip"<<std::endl;
                    }
					rs->kr_length = kr_length_; //in MappedPlant. define distance to root tipe where kr > 0 as cannot compute distance from age in case of carbon-limited growth
					rs->calcExchangeZoneCoefs();	//computes coefficient used by XylemFlux::kr_RootExchangeZonePerType
					kr_f  = std::bind(&XylemFlux::kr_RootExchangeZonePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
				}else{
					kr_f  = std::bind(&XylemFlux::kr_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
				}
                if(verbose)
				{
					std::cout << "Kr is constant per subtype of organ type, for root, subtype 0 = " << values[0][0] << " 1 day-1 \n";
                }
            }
        }
    } else {
        kr_f  = std::bind(&XylemFlux::kr_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if(verbose)
        {
            std::cout << "Kr is equal for all organs and age dependent\n";
        }
    }
}

/**
 *  Sets the axial conductivity in [cm3 day-1]
 * TODO: make deprecated: in the examples, replace setKx[Kx] by setKx[[Kx]]
 */
//either age or type/subtype dependent
void XylemFlux::setKx(std::vector<double> values, std::vector<double> age, bool verbose) {
    kx = { values }; // because kx is std::vector<std::vector<double>>
    kx_t = { age };
    if (age.size()==0) {
        if (values.size()==1) {
            kx_f = std::bind(&XylemFlux::kx_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
            if(verbose)
			{
				std::cout << "Kx is constant " << values[0] << " cm3 day-1 \n";
            }
        } else {
            kx_f  = std::bind(&XylemFlux::kx_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
            if(verbose)
			{
				std::cout << "Kx is constant per subtype, subtype 0 = " << values[0] << " cm3 day-1 \n";
            }
        }
    } else {
        kx_f  = std::bind(&XylemFlux::kx_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if(verbose)
        {
            std::cout << "Kx is age dependent\n";
        }
    }
}

//either age or type/subtype dependent
void XylemFlux::setKx(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, bool verbose)
{
    kx = values;
    kx_t = age;
    if (age.size()==0)
    {
        if (values.size()==1)
        {
            if (values[0].size()==1)
            {
                kx_f = std::bind(&XylemFlux::kx_const, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose)
                {
                    std::cout << "Kx is constant " << values[0][0] << " cm3 day-1 \n";
                }
            } else
            {
                kx_f  = std::bind(&XylemFlux::kx_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose)
                {
                    std::cout << "Kx is constant per subtype, subtype 0 = " << values[0][0] << " cm3 day-1 \n";
                }
            }
        } else
        {
            if (values[0].size()==1)
            {
                kx_f = std::bind(&XylemFlux::kx_perOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose)
                {
                    std::cout << "Kx is constant per organ type, organ type 2 (root) = " << values[0][0] << " cm3 day-1 \n";
                }
            } else
            {
                kx_f  = std::bind(&XylemFlux::kx_perType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
                if(verbose)
                {
                    std::cout << "Kx is constant per subtype of organ type, for root, subtype 0 = " << values[0][0] << " cm3 day-1 \n";
                }
            }
        }
    } else
    {
        kx_f  = std::bind(&XylemFlux::kx_table, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if(verbose)
        {
            std::cout << "Kx is equal for all organs and age dependent\n";
        }
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
void XylemFlux::setKrTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, bool verbose, bool ageBased) {
    krs= { values };
    krs_t = { age };
    kr_f = std::bind(&XylemFlux::kr_tablePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    if(verbose)
	{
		std::cout << "Kr is age dependent per root type\n";
    }
}

/**
 *  Sets the axial conductivity in [cm3 day-1]
 *	@param values 			kx values for age (its linearly interpolated between these values) for each root type
 *	@param age 				ages for the given values for each root type
 *
 * TODO: make deprecated (i would leave it in for now, used in pyhton_modules/root_conductivities.py)
 */
//both age and type/subtype dependent
void XylemFlux::setKxTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age, bool verbose) {
    kxs = {values};
    kxs_t = {age};
    kx_f = std::bind(&XylemFlux::kx_tablePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    if(verbose)
	{
		std::cout << "Kx is age dependent per root type\n";
    }
}

/**
 *  Sets the radial conductivity in [1 day-1] age, organ type and sub-type dependent
 *
 *	@param values 			kr values for age (its linearly interpolated between these values) for each organ type and each sub-type
 *	@param age 				ages for the given values for each organ type and for each sub type
 */
void XylemFlux::setKrTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age, bool verbose,  bool ageBased) {
    krs = values;
    krs_t = age;
    if (age[0].size()==1) {
        kr_f = std::bind(&XylemFlux::kr_tablePerOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if(verbose)
        {
            std::cout << "Kr is age dependent per organ type\n";
        }
    }
    else{
		if (ageBased)
		{
			kr_f = std::bind(&XylemFlux::kr_tablePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		} else {
			rs->kr_length = 100000.; //in MappedPlant. define distance to root tipe where kr > 0 as cannot compute distance from age in case of carbon-limited growth
			rs->calcExchangeZoneCoefs();	//computes coefficient used by XylemFlux::kr_RootExchangeZonePerType
			kr_f  = std::bind(&XylemFlux::kr_tablePerType_distance, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
		}
        if(verbose)
        {
            std::cout << "Kr is age dependent per organ type and sub type\n";
        }
    }
}

/**
 *  Sets the axial conductivity in [cm3 day-1] age, organ type and sub-type dependent
 *
 *	@param values 			kx values for age (its linearly interpolated between these values) for each organ type and each sub type
 *	@param age 				ages for the given values for each organ type and for each sub-type
 */
void XylemFlux::setKxTables(std::vector<std::vector<std::vector<double>>> values, std::vector<std::vector<std::vector<double>>> age, bool verbose) {
    kxs= values;
    kxs_t = age;
    if (age[0].size()==1) {
        kx_f = std::bind(&XylemFlux::kx_tablePerOrgType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if(verbose)
        {
            std::cout << "Kx is age dependent per organ type\n";
        }
    }
    else {
        kx_f = std::bind(&XylemFlux::kx_tablePerType, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        if(verbose)
        {
            std::cout << "Kx is age dependent per organ type and sub type\n";
        }
    }
}

/**
 * Sets the radial conductivity conductivity [1 day-1] per segment (e.g. constant value per segment)
 */
void XylemFlux::setKrValues(std::vector<double> values, bool verbose) {
    assert(values.size() == rs->segments.size() && "XylemFlux::setKrValues: values size must equal number of segments");
    kr.clear();
    kr_t.clear();
    kr.push_back(values);
    kr_f = std::bind(&XylemFlux::kr_valuePerSegment, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    if(verbose)
    {
		std::cout << "Kr is given per segment\n";
    }
}

/**
 * Sets the axial conductivity [cm3 day-1] per segment (e.g. constant value per segment)
 */
void XylemFlux::setKxValues(std::vector<double> values, bool verbose) {
    assert(values.size() == rs->segments.size() && "XylemFlux::setKxValues: values size must equal number of segments");
    kx.clear();
    kx_t.clear();
    kx.push_back(values);
    kx_f = std::bind(&XylemFlux::kx_valuePerSegment, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
    if(verbose)
	{
		std::cout << "Kx is given per segment\n";
    }
}

/**
 * Returns radial conductivities per segment multiplied by segment surface for a specific simulation time (TODO numleaf is ingored)
 */
std::vector<double> XylemFlux::getEffKr(double simtime) {
    std::vector<double> kr = std::vector<double>(rs->segments.size());
    for (int si = 0; si<rs->segments.size(); si++) {
        int i = rs->segments[si].x;
        int j = rs->segments[si].y;
        Vector3d n1 = rs->nodes[i];
        Vector3d n2 = rs->nodes[j];
        double l = (n2.minus(n1)).length();
        double a = rs->radii[si];
        int organType = rs->organTypes[si];
        double age = simtime - rs->nodeCTs[j];
        int subType = rs->subTypes[si];
        try {
            kr[si] = 2.*M_PI *a*l*kr_f(si, age, subType, organType);
        } catch(...) {
            std::cout << "\n XylemFlux::segFluxes: radial conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }
    }
    return kr;
}

/**
 * Returns radial conductivities per segment multiplied by segment surface for a specific simulation time (TODO numleaf is ingored)
 */
std::vector<double> XylemFlux::getKr(double simtime) {
    std::vector<double> kr = std::vector<double>(rs->segments.size());
    for (int si = 0; si<rs->segments.size(); si++) {
        int j = rs->segments[si].y;
        int organType = rs->organTypes[si];
        double age = simtime - rs->nodeCTs[j];
        int subType = rs->subTypes[si];
        try {
            kr[si] = kr_f(si, age, subType, organType);
        } catch(...) {
            std::cout << "\n XylemFlux::segFluxes: radial conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }
    }
    return kr;
}


/**
 * Returns radial conductivities per segment for a specific simulation time
 */
std::vector<double> XylemFlux::getKx(double simtime) {
    std::vector<double> kx = std::vector<double>(rs->segments.size());
    for (int si = 0; si<rs->segments.size(); si++) {
        int j = rs->segments[si].y;
        int organType = rs->organTypes[si];
        double age = simtime - rs->nodeCTs[j];
        int subType = rs->subTypes[si];
        try {
            kx[si] = kx_f(si, age, subType, organType);
        } catch(...) {
            std::cout << "\n XylemFlux::segFluxes: axial conductivities failed" << std::flush;
            std::cout  << "\n organ type "<<organType<< " subtype " << subType <<std::flush;
        }
    }
    return kx;
}

/**
 * Returns soil matric potential per segment, for a given soil sx connected gy the mapper rs->seg2cell
 */
std::vector<double> XylemFlux::getHs(const std::vector<double>& sx) {
    std::vector<double> hs = std::vector<double>(rs->segments.size());
    for (int si = 0; si<rs->segments.size(); si++) {
        int cellIndex = rs->seg2cell[si];
        if (cellIndex>=0) {
            if(sx.size()>1) {
                hs[si] = sx.at(cellIndex);
            } else {
                hs[si] = sx.at(0);
            }
        } else {
            hs[si] = psi_air;
        }
    }
    return hs;
}


/**
 * Returns kr of roots belowground or leaves aboveground, overwise returns 0
 */
double XylemFlux::kr_f_wrapped(int si, double age, int subType, int organType, bool cells) const
{
	
	int cellIndex = 0;
	if (cells){
		cellIndex= rs->seg2cell.at(si);
	}
	if (cells&&(((cellIndex>=0)&&(organType !=Organism::ot_root))||((cellIndex < 0)&&(organType ==Organism::ot_root))))
	{
			return 0.;
	}else
	{
		return kr_f(si, age, subType, organType);
	}
}

} // namespace
