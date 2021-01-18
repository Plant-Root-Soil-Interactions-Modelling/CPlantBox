#ifndef EXUDATIONMODEL2_H
#define EXUDATIONMODEL2_H

#include "external/gauss_legendre/gauss_legendre.h"
#include "soil.h"
#include "sdf_rs.h"
#include "RootSystem.h"

#include <functional>


namespace CPlantBox {

/**
 * docme
 */
class ExudationModel2 {
public:

	enum IntegrationType { mps_straight = 0, mps = 1, mls = 2 };

	/*
	 * Model parameters (same for all roots)
	 */
	double Q = 1e-5;
	double Dl = 1e-5; // cm2 / day
	double theta = 0.3;
	double R = 1;
	double k = 1e-6;
	double l = 0.1; // cm

	/*
	 *  Numerical parameters
	 */
	EquidistantGrid3D grid;
	int type = mps;
	int n0 = 5; // integration points per [cm]
	double thresh13 = 1.e-15; // threshold for Eqn 13
	bool calc13 = true; // turns Eqn 13 on and off
	double observationRadius = 5; //  limits computational domain around roots [cm]
	double eps = 1.e-6; // minimal distance from grid point

	/**
	 * Constructors
	 *
	 */
	ExudationModel2(double width, double depth, int n, std::shared_ptr<RootSystem> rs) :ExudationModel2(width, width, depth, n, n, n, rs) { }

	ExudationModel2(double length, double width, double depth, int nx, int ny, int nz, std::shared_ptr<RootSystem> rs) :grid(EquidistantGrid3D(length, width, depth, nx, ny, nz, false)) {

		dx3 = (length/nx)*(width/ny)*(depth/nz); // for integration of eqn 13
		roots = rs->getRoots();

		for (const auto& r : roots) {
			if (r->getNumberOfNodes()>1) { // started growing
				// time when the root stopped growing
				double sTime = r->getNodeCT(r->getNumberOfNodes()-1);
				if (r->getLength()<(r->param()->getK()-0.1)) { // is active?
					stopTime.push_back(0);
				} else {
					stopTime.push_back(sTime);
				}
				// root tip
				Vector3d t = r->getNode(r->getNumberOfNodes()-1);
				tip.push_back(t);
				// direction towards root base
				Vector3d base = r->getNode(0);
				double a = r->getNodeCT(r->getNumberOfNodes()-1) - r->getNodeCT(0);
				v.push_back(base.minus(t).times(1./a));
			}
			sdfs.push_back(SDF_RootSystem(*r, observationRadius));
		}

		n_size = nx*ny*nz;
	}

	/**
	 * makes a list of voxels for each root, that lie within the observation radius of the root
	 */
	void makeVoxelLists(int i0 = 0, int iend = -1) {

		voxelList.clear();

		if (iend==-1) {
			iend = roots.size();
		}
		i0_ = i0;
		iend_ = iend;

		for (size_t ri = i0; ri< iend; ri++) {
			int id = roots[ri]->getId();
			//            // find bounding box of root id
			//            int n = roots[ri]->getNumberOfNodes();
			//            std::vector<double> x(n),y(n),z(n);
			//            for (int i=0; i<n; i++) {
			//                auto n = roots[ri]->getNode(i);
			//                x[i] = n.x;
			//                y[i] = n.y;
			//                z[i] = n.z;
			//            }
			//            auto minmaxX = std::minmax_element(x.begin(), x.end());
			//            auto minmaxY = std::minmax_element(y.begin(), y.end());
			//            auto minmaxZ = std::minmax_element(z.begin(), z.end());
			//            int minX = std::max(int(*minmaxX.first - observationRadius + 0.5),0);
			//            int minY = std::max(int(*minmaxY.first - observationRadius + 0.5),0);
			//            int minZ = std::max(int(*minmaxZ.first - observationRadius + 0.5),0);
			//            int maxX = std::min(int(*minmaxX.second + observationRadius + 0.5),(int)grid.nx);
			//            int maxY = std::min(int(*minmaxY.second + observationRadius + 0.5),(int)grid.ny);
			//            int maxZ = std::min(int(*minmaxZ.second + observationRadius + 0.5),(int)grid.nz);

			//            std::cout << "creating voxel list for root " << id << "\n";
			//            for (size_t i = minX; i<maxX; i++) {
			//                for(size_t j = minY; j<maxY; j++) {
			//                    for (size_t k = minZ; k<maxZ; k++) {
			for (size_t i = 0; i<grid.nx; i++) {
				for(size_t j = 0; j<grid.ny; j++) {
					for (size_t k = 0; k<grid.nz; k++) {
						size_t lind = k*(grid.nx*grid.ny)+j*grid.nx+i; // k*(nx*ny)+j*nx+i, same ordering as RectilinearGrid3D
						x_ = grid.getGridPoint(i,j,k);
						if (-sdfs[ri].getDist(x_)<observationRadius) {
							if (voxelList.count(id)==0) {
								voxelList[id] = std::vector<int>();
							}
							voxelList.at(id).push_back(lind);
						}
					}
				}
			}
			std::cout << "made list " << id <<"/" << iend << " length " << roots[ri]->getLength() <<
					" cm , vol "<< roots[ri]->getLength()*M_PI*M_PI*observationRadius << " cm3, having " << voxelList.at(id).size() << " voxels\n";
		}
	}

	/**
	 * adds voxel results to solution c
	 */
	std::vector<double> addResults(std::vector<double>& c) {
		for (int lind=0; lind<n_size; lind++) {
			if (c_.count(lind)>0) {
				c[lind] += c_[lind];
			}
		}
		return c;
	}

	/**
	 * For each root for each grid point
	 * @param tend      final simulation time
	 * @param i0        optionally, initial root index (default = 0)
	 * @param iend      optionally, final root index (default = roots.size())
	 */
	std::vector<double> calculate(double tend, int i0 = 0, int iend = -1) { //

		if (iend==-1) {
			iend = roots.size();
		}

		c_.clear();
		g_.clear();

		for (size_t ri = i0; ri< iend; ri++) {

			r_ = roots[ri];
			age_ = std::min(r_->getNodeCT(r_->getNumberOfNodes()-1),tend) - r_->getNodeCT(0);

			if (age_>0) {

				// per root (passed to integrands)
				n_ = int(n0*r_->getLength()); // number of integration points eq 11
				v_ = v[ri]; // for mps_straight, eq 11
				tip_ = tip[ri]; // for mps_straight, eq 11
				st_ = stopTime[ri]; // eq 13
				st_ *= calc13;

				double age;
				if (st_>0) {
					double stopAge = st_ - r_->getNodeCT(0);
					age = std::min(age_, stopAge);
				} else {
					age = age_;
				}

				std::cout << "Root #" << ri << "/" << roots.size() << ", age_ "<< age_ <<  ", age "<< age <<
						", st_ "<< st_ << ", nodeCT(0) "<< r_->getNodeCT(0) << ", tend "<< tend <<
						", res "<< n_ << " \n"; // for debugging

				// EQN 11
				int id = r_->getId();
				for (int lind : voxelList[id]) {
					// different flavors of Eqn (11)
					x_ = grid.getGridPoint(lind);
					double c = eqn11(0, age, 0, l); // needs x_!
					if (c_.count(lind)==0) {
						c_[lind] = 0.;
					}
					c_[lind] += c;
					g_[lind] = c;
				}

				// EQN 13
				if ((st_>0+1.e-10) && (st_+1.e-10<tend)) { // has stopped growing
					std::cout << "13!";
					for (int lind : voxelList[id]) {
						if (g_[lind] > thresh13) {
							x_ = grid.getGridPoint(lind);
							c_[lind] += integrate13(tend, id); // needs x_!
						}
					}
					std::cout << "13\n";
				}
			} // if ages.at(i)>0
		}

		std::vector<double> c = std::vector<double>(n_size);
		std::fill(c.begin(), c.end(), 0.);
		return addResults(c);
	}

	double eqn11(double x0, double xend, double y0, double yend) {
		switch (type) {
		case mps_straight: {
			return gauss_legendre(n_, integrandMPS_straight, this, x0, xend);
		}
		case mps: {
			return gauss_legendre(n_, integrandMPS, this, x0, xend);
		}
		case mls: {
			return gauss_legendre_2D_cube(n_, integrandMLS, this, x0, xend, y0, yend);
		}
		}
		std::cout << "Unknown integration type \n";
		return 0.;
	}

	// simplistic integration in 3d
	double integrate13(double tend, int id) { // called for fixed x_
		double c = 0;
		for (int lind : voxelList[id]) {
			Vector3d y = grid.getGridPoint(lind);
			c += integrand13(y, lind, tend)*dx3; // depends on (fixed) x_
		}
		//		for (int lind= 0; lind<n_size; lind++) { // full integral
		//			Vector3d y = grid.getGridPoint(lind);
		//			c += integrand13(y, lind, t)*dx3; // depends on (fixed) x_
		//		}
		return c;
	}

	// integrand Eqn 13
	double integrand13(Vector3d& y, size_t lind, double tend) {
		double dt = tend-st_;
		double c = to32(R)*g_[lind] / to32(4*Dl*M_PI*dt);
		Vector3d z = x_.minus(y);
		double l = z.times(z);
		l = std::max(eps,l);
		return c*exp(-R/(4*Dl*dt) * l - k*dt/R); // note that dt -> 0 => integrand13 -> 0 if t.times(z)>0
	}

	// Returns the linearly interpolated position along the root r at age a
	static Vector3d pointAtAge(std::shared_ptr<Root> r, double a) {
		a = std::max(0.,a);
		double et = r->getNodeCT(0)+a; // age -> emergence time
		size_t i=0;
		while (i<r->getNumberOfNodes()) {
			if (r->getNodeCT(i)>et) { // first index bigger than emergence time, interpolate i-1, i
				break;
			}
			i++;
		}
		if (i == r->getNumberOfNodes()) { // this happens if a root has stopped growing
			std::cout << "pointAtAge(): warning age is older than the root \n";
			return r->getNode(i-1);
		}
		Vector3d n1 = r->getNode(i-1);
		Vector3d n2 = r->getNode(i);
		double t = (et - r->getNodeCT(i - 1)) / (r->getNodeCT(i) - r->getNodeCT(i - 1)); // t in (0,1]
		return (n1.times(1. - t)).plus(n2.times(t));
	}

	static double to32(double x) { return sqrt(x*x*x); }

	static double to3(double x) { return x*x*x; }

	// point source, root is represented by a single straight line (substituted)
	static double integrandMPS_straight(double t, void* param) {
		ExudationModel2* p = (ExudationModel2*) param;
		double c = -p->R / ( 4*p->Dl*t );
		double d = 8*(p->theta)*ExudationModel2::to32(M_PI*p->Dl*t);

		Vector3d xtip = p->tip_.plus(p->v_.times(t)); // for t=0 at tip, at t=age at base, as above
		Vector3d z = p->x_.minus(xtip);

		double l_ = z.times(z);
		l_ = std::max(p->eps,l_);

		return ((p->Q)*sqrt(p->R))/d *exp(c*l_- p->k/p->R * t); // Eqn (11)
	}

	// moving line source, root is represented by a straight segments
	static double integrandMLS(double t, double l, void* param) {
		ExudationModel2* p = (ExudationModel2*) param;
		double c = -(p->R) / ( 4*(p->Dl)*t );
		double d = 8*(p->theta)*ExudationModel2::to32(M_PI*p->Dl*t);

		double tl = p->r_->calcLength( p->age_-t ); // tip
		if (tl<l) { // if root smaller l
			return 0.;
		}
		double agel = p->r_->calcAge(tl-l);
		Vector3d tipLS = p->ExudationModel2::pointAtAge(p->r_, agel);
		Vector3d z = p->x_.minus(tipLS);

		double l_ = z.times(z);
		l_ = std::max(p->eps,l_);

		return ((p->Q)*sqrt(p->R))/d *exp(c*l_ - p->k/p->R * t); // Eqn (11)
	}

	// moving point source, root is represented by a straight segments
	static double integrandMPS(double t, void* param) {
		ExudationModel2* p = (ExudationModel2*) param;
		double c = -p->R / ( 4*p->Dl*t );
		double d = 8*(p->theta)*ExudationModel2::to32(M_PI*p->Dl*t);

		Vector3d xtip = ExudationModel2::pointAtAge(p->r_, p->age_-t);
		Vector3d z = p->x_.minus(xtip);

		double l_ = z.times(z);
		l_ = std::max(p->eps,l_);

		return ((p->Q)*sqrt(p->R))/d *exp(c*l_ - p->k/p->R * t); // Eqn (11)
	}

	// Root system
	std::vector<std::shared_ptr<Root>> roots;
	std::vector<double> stopTime; // time when root stopped growing, 0 if it has not
	std::vector<Vector3d> tip;
	std::vector<Vector3d> v; // direction from tip towards root base
	double dx3 = 1;
	std::vector<SDF_RootSystem> sdfs; // direction from tip towards root base
	std::map<int, std::vector<int>> voxelList; // voxel list per root
	int i0_ = 0;
	int iend_ = 0;

	// Set before integrating
	Vector3d x_ = Vector3d(); // integration point
	int n_ = 0;
	std::shared_ptr<Root> r_ = nullptr; // current root
	double age_ = 0;
	Vector3d tip_ = Vector3d();
	Vector3d v_ = Vector3d();
	double st_ = 0; // stop time (eqn 13)

	int n_size = 0;
	//	std::vector<double> g_;
	//	std::vector<double> c_;
	std::map<int, double> g_;  // eqn 13
	std::map<int, double> c_;

};


}

#endif
