// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "XylemFlux.h"

#include <algorithm>

namespace CPlantBox {

/**
 * Assembles the linear system as sparse matrix, given by public member variables,
 * indices aI, aJ, and corresponding values aV; and load aB
 *
 * @param simTime         	[days] current simulation time is needed for age dependent conductivities,
 *                  		to calculate the age from the creation times.
 * @param sx 				soil matric potential in the cells or around the segments
 * @param cells 			sx per cell (true), or segments (false)
 */
void XylemFlux::linearSystem(double simTime, const std::vector<double>& sx, bool cells)
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

		double psi_s;
		if (cells) { // soil matric potential given per cell
			psi_s = sx.at(rs->seg2cell[j-1]); // segIdx = s.y-1
		} else {
			psi_s = sx.at(j-1); // segIdx = s.y-1
		}

		double a = rs->radii[si]; // si is correct, with ordered and unordered segmetns
		double age = simTime - rs->nodeCTs[j];
		int type = rs->types[si];
		double kx = kx_f(age, type);
		double  kr = kr_f(age, type);

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
 * @param simTime   [days] current simulation time (to calculate age dependent conductivities)
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

			double f =  -2*a*M_PI*(kr*rho*g); // flux is proportional to f
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
 * @param simTime   [days] current simulation time (to calculate age dependent conductivities)
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
		int segIdx = j-1;
		double psi_s = sx.at(segIdx);

		double a = rs->radii[si]; // si is correct, with ordered and unordered segments
		double age = simTime - rs->nodeCTs[j];
		int type = rs->types[si];
		double  kr = kr_f(age, type);
		double  kz = kx_f(age, type);

		auto n1 = rs->nodes[i];
		auto n2 = rs->nodes[j];
		double l = (n2.minus(n1)).length();

		double f =  -2*a*M_PI*(kr*rho*g); // flux is proportional to f
		double fApprox = f*l*(psi_s - rx[j]); // cm3 / day

		double tau = std::sqrt(2*a*M_PI*kr/kz); // sqrt(c) [cm-1]
		double d = std::exp(-tau*l)-std::exp(tau*l); // det
		double fExact = -f*(1./(tau*d))*(rx[i]-psi_s+rx[j]-psi_s)*(2.-std::exp(-tau*l)-std::exp(tau*l));

		double flux = fExact*(!approx)+approx*fApprox;
		fluxes[i] = flux;
	}
	return fluxes;
}

/**
 * TODO
 */
std::map<int,double> XylemFlux::sumSoilFluxes(std::vector<double> segFluxes)
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
 *  Sets the radial conductivity in [1 day-1], converts to [cm2 day g-1] by dividing by rho*g
 */
void XylemFlux::setKr(std::vector<double> values, std::vector<double> age) {
	std::transform(values.begin(), values.end(), values.begin(), std::bind1st(std::multiplies<double>(),1./(rho * g)));
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
 *  Sets the axial conductivity in [cm3 day-1], converts to [cm5 day g-1] by dividing by rho*g
 */
void XylemFlux::setKx(std::vector<double> values, std::vector<double> age) {
	std::transform(values.begin(), values.end(), values.begin(), std::bind1st(std::multiplies<double>(),1./(rho * g)));
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
 *  Sets the radial conductivity in [1 day-1], converts to [cm2 day g-1] by dividing by rho*g
 */
void XylemFlux::setKrTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age) {
	krs.resize(0);
	for (auto v :values) {
		std::transform(v.begin(), v.end(), v.begin(), std::bind1st(std::multiplies<double>(),1./(rho * g)));
		krs.push_back(v);
	}
	krs_t = age;
	kr_f = std::bind(&XylemFlux::kr_tablePerType, this, std::placeholders::_1, std::placeholders::_2);
	std::cout << "Kr is age dependent per root type\n";
}

/**
 *  Sets the axial conductivity in [cm3 day-1], converts to [cm5 day g-1] by dividing by rho*g
 */
void XylemFlux::setKxTables(std::vector<std::vector<double>> values, std::vector<std::vector<double>> age) {
	kxs.resize(0);
	for (auto v :values) {
		std::transform(v.begin(), v.end(), v.begin(), std::bind1st(std::multiplies<double>(),1./(rho * g)));
		kxs.push_back(v);
	}
	kxs_t = age;
	kx_f = std::bind(&XylemFlux::kx_tablePerType, this, std::placeholders::_1, std::placeholders::_2);
	std::cout << "Kx is age dependent per root type\n";
}

/**
 *
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
