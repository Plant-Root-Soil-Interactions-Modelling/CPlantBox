// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "XylemFlux.h"

#include <algorithm>

// #include "../external/gauss_legendre/gauss_legendre.h"

namespace CPlantBox {

/**
 * Assembles the linear system as sparse matrix, given by public member variables,
 * indices aI, aJ, and corresponding values aV; and load aB
 *
 * @param simTime         	[days] current simulation time is needed for age dependent conductivities,
 *                  		to calculate the age from the creation times (age = sim_time - segment creation time).
 * @param sx 				soil matric potential in the cells or around the segments
 * @param cells 			sx per cell (true), or segments (false)
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
	for (int si = 0; si<Ns; si++) {

		int i = rs->segments[si].x;
		int j = rs->segments[si].y;

		double psi_s = 0.;
		if (cells) { // soil matric potential given per cell
			try {
		    psi_s = sx.at(rs->seg2cell[j-1]); // segIdx = s.y-1
			} catch(...) {
			  std::cout << "mapping failed\n" << std::flush;
			}
		} else {
			psi_s = sx.at(j-1); // segIdx = s.y-1
		}

		double a = rs->radii[si]; // si is correct, with ordered and unordered segmetns
		double age = simTime - rs->nodeCTs[j];
		int type = rs->types[si];
        double kx = 0.;
        double  kr = 0.;
		try {
		    kx = kx_f(age, type);
		    kr = kr_f(age, type);
        } catch(...) {
            std::cout << "conductivities failed\n" << std::flush;
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
 * Fluxes from root segments into soil cells
 *
 * @param simTime   [days] current simulation time is needed for age dependent conductivities,
 *                  to calculate the age from the creation times (age = sim_time - segment creation time).
 * @param rx        [cm] root xylem matric potential
 * @param sx        [cm] soil matric potential for eadh cell
 * @param approx    approximate or exact (default = false, i.e. exact)
 *
 * @return Hash map with cell indices as keys and fluxes as values [cm3/day]
 */
std::map<int,double> XylemFlux::soilFluxes(double simTime, const std::vector<double>& rx, const std::vector<double>& sx, bool approx)
{
	std::map<int,double> fluxes;

	for (int si = 0; si<rs->segments.size(); si++) {

		int i = rs->segments[si].x;
		int j = rs->segments[si].y;
		int segIdx = j-1;

		if (rs->seg2cell.count(segIdx)>0) {

			int cellIdx = rs->seg2cell[segIdx];
			double psi_s = sx.at(cellIdx);

			double a = rs->radii[si]; // si is correct, with ordered and unordered segments
			double age = simTime - rs->nodeCTs[j];
			int type = rs->types[si];
			double  kr = kr_f(age, type);
			double  kz = kx_f(age, type);

			auto n1 = rs->nodes[i];
			auto n2 = rs->nodes[j];
			double l = (n2.minus(n1)).length();

			double f =  -2*a*M_PI*kr; // flux is proportional to f // *rho*g
			double fApprox = f*l*(psi_s - rx[j]); // cm3 / day

			double tau = std::sqrt(2*a*M_PI*kr/kz); // sqrt(c) [cm-1]
			double d = std::exp(-tau*l)-std::exp(tau*l); // det
			double fExact = -f*(1./(tau*d))*(rx[i]-psi_s+rx[j]-psi_s)*(2.-std::exp(-tau*l)-std::exp(tau*l));

			double flux = fExact*(!approx)+approx*fApprox;
			// std::cout << cellIdx << ", " << fExact << ", " << fApprox << ", psi_s " << psi_s << "\n" << std::flush;
			if (fluxes.count(cellIdx)==0) {
				fluxes[cellIdx] = flux;
			} else {
				fluxes[cellIdx] += flux; // sum up fluxes per cell
			}

		} else {
			std::cout << "XylemFlux::soilFluxes: Warning! unmapped segments with index " << segIdx << "\n";
		}

	}
	return fluxes;
}

/**
 * Fluxes from root segments into soil cells
 *
 * @param simTime   [days] current simulation time is needed for age dependent conductivities,
 *                  to calculate the age from the creation times (age = sim_time - segment creation time).
 * @param rx        [cm] root xylem matric potential
 * @param sx        [cm] soil matric potential for each segment
 * @param approx    approximate or exact (default = false, i.e. exact)
 *
 * @return Hash map with cell indices as keys and fluxes as values [cm3/day]
 */
std::vector<double> XylemFlux::segFluxes(double simTime, const std::vector<double>& rx, const std::vector<double>& sx, bool approx)
{
	std::vector<double> fluxes = std::vector<double>(rs->segments.size());
	for (int si = 0; si<rs->segments.size(); si++) {

		int i = rs->segments[si].x;
		int j = rs->segments[si].y;

		if (j>=rs->nodeCTs.size()) {
		    std::cout << j << " index j is broken \n"<< std::flush;
		}

		double psi_s = sx.at(si);

		double a = rs->radii[si]; // si is correct, with ordered and unordered segments
		double age = simTime - rs->nodeCTs[j];
		int type = rs->types[si];
		double  kr = kr_f(age, type);
		double  kz = kx_f(age, type);

        Vector3d n1 = rs->nodes[i];
        Vector3d n2 = rs->nodes[j];
		double l = (n2.minus(n1)).length();

		double f =  -2*a*M_PI*kr; // flux is proportional to f // *rho*g
		double fApprox = f*l*(psi_s - rx[j]); // cm3 / day

		double tau = std::sqrt(2*a*M_PI*kr/kz); // sqrt(c) [cm-1]
		double d = std::exp(-tau*l)-std::exp(tau*l); // det
		double fExact = -f*(1./(tau*d))*(rx[i]-psi_s+rx[j]-psi_s)*(2.-std::exp(-tau*l)-std::exp(tau*l));

		double flux = fExact*(!approx)+approx*fApprox;
		fluxes[si] = flux;
	}
	return fluxes;
}

/**
 * Applies source term (operator splitting)
 * Only apllies soilFluxes if not stressed according to approximation of Schröder et al.
 * If stressed applies Flux based assuming critP in xylem
 */
std::vector<double> XylemFlux::segFluxesSchroeder(double simTime, std::vector<double> rx, const std::vector<double>& sx, double critP, std::function<double(double)> mfp) {
	auto outerRadii = this->segOuterRadii();
	auto fluxes = this->segFluxes(simTime, rx, sx, false);
	for (int i = 0; i<rs->segments.size(); i++) { // modify rx in case of stress
		double q_root = fluxes.at(i); // UNITS?
		double q_out = 0.;
		double r_root = rs->radii.at(i);
		double r_out = outerRadii.at(i);
		double rho = r_out/r_root;
		double r = r_root;
		double p = 0;
		if (rs->seg2cell.count(i)>0) {
			int cellIdx = rs->seg2cell.at(i);
			p = sx[cellIdx];
		} else {
			std::cout << "XylemFlux::segFluxesSchroeder: Unmapped segment \n";
		}

		double mfp_ = mfp(p);
		double noStress = mfp_ + (q_root*r_root-q_out*r_out)*((r*r)/(r_root*r_root)/(2*(1-rho*rho)))+
				(rho*rho)/(1-(rho*rho)*(log(r_out/r)-0.5)) + q_out*r_out*log(r/r_out);

		if (noStress<0) { // in case of stress, modify xylem pressure
			rx[i] = critP; // SHOULD BE SX
		}
	}
	return this->segFluxes(simTime, rx, sx, false); // exact fluxes, due to new rx
}

/**
 * Calculates outer segment radii, so that the summed segment volumes per cell equal the cell volume
 */
std::vector<double> XylemFlux::segOuterRadii(int type) const {
	auto lengths =  this->segLength();
    auto width = rs->maxBound.minus(rs->minBound);
	double cellVolume = width.x*width.y*width.z/rs->resolution.x/rs->resolution.y/rs->resolution.z; // TODO only true for equidistant rectangular grid
	std::vector<double> radii = std::vector<double>(rs->segments.size());
	std::fill(radii.begin(), radii.end(), 0.);
	auto& map = rs->cell2seg;
	for(auto iter = map.begin(); iter != map.end(); ++iter) {
		int cellId =  iter->first;
		auto segs = map.at(cellId);
		double v = 0.;  // calculate sum of root volumes or surfaces over cell
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
			double l = lengths[i];
			double t =0.; // proportionality factor (must sum up to == 1 over cell)
			if (type==0) { // volume
				t = M_PI*(rs->radii[i]*rs->radii[i])*l/v;
			} else if (type==1) { // surface
				t = 2*M_PI*rs->radii[i]*l/v;
			} else if (type==2) { // length
                t = l/v;
            }
			double targetV = t * cellVolume;  // target volume
			radii[i] = sqrt(targetV/(M_PI*l)+rs->radii[i]*rs->radii[i]);
		}
	}
	return radii;
}


/**
 * Sums segment fluxes over each cell, @see splitSoilFluxes()
 *
 * @param segFluxes 	segment fluxes [cm3/day]
 * @return fluxes for each cell idx [cm3/day]
 */
std::map<int,double> XylemFlux::sumSoilFluxes(const std::vector<double>& segFluxes)
{
	std::map<int,double> fluxes;
	for (int si = 0; si<rs->segments.size(); si++) {
		int j = rs->segments[si].y;
		int segIdx = j-1;
		if (rs->seg2cell.count(segIdx)>0) {
			int cellIdx = rs->seg2cell[segIdx];
			if (fluxes.count(cellIdx)==0) {
				fluxes[cellIdx] = segFluxes[segIdx];
			} else {
				fluxes[cellIdx] += segFluxes[segIdx]; // sum up fluxes per cell
			}
		}
	}
	return fluxes;
}

/**
 * Splits soil fluxes per cell to the segments within the cell, so that the summed fluxes agree, @see sumSoilFluxes()
 *
 * @param soilFluxes 	cell fluxes per global index [cm3/day]
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
 * Calculates segment lengths
 */
std::vector<double> XylemFlux::segLength() const {
	std::vector<double> lengths = std::vector<double>(rs->segments.size());
	for(int i=0; i<lengths.size(); i++) {
		auto n1 = rs->nodes[rs->segments[i].x];
		auto n2 = rs->nodes[rs->segments[i].y];
		lengths[i] = (n2.minus(n1)).length();
	}
	return lengths;
}

/**
 *  Sets the radial conductivity in [1 day-1]
 */
void XylemFlux::setKr(std::vector<double> values, std::vector<double> age) {
	kr = values;
	kr_t = age;
	if (age.size()==0) {
		if (values.size()==1) {
			kr_f = std::bind(&XylemFlux::kr_const, this, std::placeholders::_1, std::placeholders::_2);
			std::cout << "Kr is constant " << values[0] << " cm2 day g-1 \n";
		} else {
			kr_f  = std::bind(&XylemFlux::kr_perType, this, std::placeholders::_1, std::placeholders::_2);
			std::cout << "Kr is constant per type, type 0 = " << values[0] << " cm2 day g-1 \n";
		}
	} else {
		kr_f  = std::bind(&XylemFlux::kr_table, this, std::placeholders::_1, std::placeholders::_2);
		std::cout << "Kr is age dependent\n";
	}
}


/**
 *  Sets the axial conductivity in [cm3 day-1]
 */
void XylemFlux::setKx(std::vector<double> values, std::vector<double> age) {
	kx = values;
	kx_t = age;
	if (age.size()==0) {
		if (values.size()==1) {
			kx_f = std::bind(&XylemFlux::kx_const, this, std::placeholders::_1, std::placeholders::_2);
			std::cout << "Kx is constant " << values[0] << " cm2 day g-1 \n";
		} else {
			kx_f  = std::bind(&XylemFlux::kx_perType, this, std::placeholders::_1, std::placeholders::_2);
			std::cout << "Kx is constant per type, type 0 = " << values[0] << " cm2 day g-1 \n";
		}
	} else {
		kx_f  = std::bind(&XylemFlux::kx_table, this, std::placeholders::_1, std::placeholders::_2);
		std::cout << "Kx is age dependent\n";
	}
}

/**
 *  Sets the radial conductivity in [1 day-1]
 */
void XylemFlux::setKrTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age) {
	krs.resize(0);
	for (auto v :values) {
		krs.push_back(v);
	}
	krs_t = age;
	kr_f = std::bind(&XylemFlux::kr_tablePerType, this, std::placeholders::_1, std::placeholders::_2);
	std::cout << "Kr is age dependent per root type\n";
}

/**
 *  Sets the axial conductivity in [cm3 day-1]
 */
void XylemFlux::setKxTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age) {
	kxs.resize(0);
	for (auto v :values) {
		kxs.push_back(v);
	}
	kxs_t = age;
	kx_f = std::bind(&XylemFlux::kx_tablePerType, this, std::placeholders::_1, std::placeholders::_2);
	std::cout << "Kx is age dependent per root type\n";
}

/**
 * Linear interpolation
 */
double XylemFlux::interp1(double ip, std::vector<double> x, std::vector<double> y) {
	if (ip > x.back()) return y.back(); // check bounds
	if (ip < x[0]) return y[0];

	// if we are within bounds find the index of the lower bound
	const auto lookUpIndex = std::distance(x.begin(), std::lower_bound(x.begin(), x.end(), ip));
	if (lookUpIndex == 0) {
		return y[0];
	}
	double ip_ = (ip - x[lookUpIndex-1])/(x[lookUpIndex] - x[lookUpIndex-1]);
	return y[lookUpIndex-1]*(1.0 - ip_)  + y[lookUpIndex]*ip_;
}

} // namespace
