// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "PlantHydraulicModel.h"

#include <algorithm>
#include <set>
//#ifdef _OPENMP
//#include <omp.h>
//#endif

namespace CPlantBox {

PlantHydraulicModel::PlantHydraulicModel(std::shared_ptr<CPlantBox::MappedSegments> ms, std::shared_ptr<CPlantBox::PlantHydraulicParameters> params):
    ms(ms), params(params) 
    { params->ms = ms;
    //#ifdef _OPENMP
    //std::cout << "OpenMP enabled, max threads: " << omp_get_max_threads()<<" "<<omp_get_num_procs() <<std::flush<< std::endl;
    //#else
    //    std::cout << "OpenMP not enabled" <<std::flush<< std::endl;
    //#endif
     //Eigen::initParallel();
     //nthreads = Eigen::nbThreads( );
     //std::cout<<"PlantHydraulicModel::PlantHydraulicModel number of threads for openmp "<< nthreads <<std::flush<<std::endl;
	 
	 //#pragma omp parallel
//{
  //  int tid = omp_get_thread_num();
    //#pragma omp critical
    //std::cout << "Thread " << tid << " running on core " << sched_getcpu() << "\n";
//}
     //Eigen::setNbThreads(1);
    }

/** prescribes a Dirichlet boundary conditions for the system Qx=b
    @param n0         list of node indices, where the Dirichlet bc is applied
    @param d [cm]     list of Dirichlet values   
    @return Q, b, the updated matrix, and rhs vector 
 */
void PlantHydraulicModel::bc_dirichlet(Eigen::SparseMatrix<double>& mat, const std::vector<int>& n0, const std::vector<double>& d)
{
    for (size_t c = 0; c < n0.size(); ++c)
    {
        int i = n0[c];
        double di = d[c];

        // subtract column contribution from RHS
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it)
        {
            int j = it.row();
            if (j != i)
                b[j] -= it.value() * di;
        }

        // zero column i
        for (int k = 0; k < mat.outerSize(); ++k)
            for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
                if (it.col() == i)
                    it.valueRef() = 0.0;

        // zero row i
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, i); it; ++it)
            it.valueRef() = 0.0;

        mat.coeffRef(i,i) = 1.0;
        b[i] = di;
    }
}

/**
 * Solves the linear system filled by @see XylemFlux::linearSystem
 *
 * @param simTime[day]  	current simulation time, needed for age dependent conductivities,
 *                  		to calculate the age from the creation times (age = sim_time - segment creation time).
 * @param sx [cm]			soil matric potential in the cells or around the segments, given per cell or per segment
 * @param cells 			sx per cell (true), or segments (false)
 * @param soil_k [day-1]    optionally, soil conductivities can be prescribed per segment,
 *                          conductivity at the root surface will be limited by the value, i.e. kr = min(kr_root, k_soil)
 */
void PlantHydraulicModel::linearSystemMeunierSolve(double simTime, const std::vector<double> sx, 
													bool cells, const std::vector<double> soil_k, 
                                                   const std::vector<int> n0, const std::vector<double> d,
												   bool verbose)
{
	dovector = false;
	int Ns = ms->segments.size(); // number of segments
    int N = ms->nodes.size(); // number of nodes
    if (Ns != N - 1)
    {
        throw std::runtime_error("mismatch in the number of nodes and segments");
    }
    psiXyl.clear();
	tripletList.clear();
	tripletList.reserve(Ns*4);
	b = Eigen::VectorXd::Zero(N);
	//get "tripletList" and "b"
	auto start = std::chrono::high_resolution_clock::now();
	linearSystemMeunier(simTime, sx, cells, soil_k); 
    for (int i = 0; i < N; ++i){  b(i) = aB[i];}
	auto end = std::chrono::high_resolution_clock::now();
	if(verbose)
	{
		std::cout<< "time spent in linearSystemMeunier : "<<  std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << " seconds\n"<<std::flush;
	}
    if (!b.allFinite()) {
        throw std::runtime_error("RHS vector b contains NaN or Inf");
    }
    
	Eigen::SparseMatrix<double> mat(N,N);
	mat.reserve(Eigen::VectorXi::Constant(N,2));
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
    bc_dirichlet(mat, n0, d);
	mat.makeCompressed();
	//Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solverBiCGSTAB;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solverLU;
    
        Eigen::VectorXd v2;
    //Eigen::setNbThreads(nthreads);
	
        start = std::chrono::high_resolution_clock::now();
        solverLU.compute(mat);
    
	if(solverLU.info() != Eigen::Success){
		std::cout << "PlantHydraulicModel::linearSystemMeunierSolve  matrix Compute with Eigen failed: " << solverLU.info() << std::endl;
		throw std::runtime_error("Photosynthesis::linearSystemSolve  matrix Compute with Eigen failed" );
	}
		end = std::chrono::high_resolution_clock::now();
		if(verbose)
		{
			std::cout<< "time spent in solverLU.compute(mat) : "<<  std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count() << " seconds\n"<<std::flush;
		}
        if(solverLU.info() != Eigen::Success){
            std::cout << "PlantHydraulicModel::linearSystemMeunierSolve  matrix Compute with Eigen failed: " << solverLU.info() << std::endl;
            throw std::runtime_error("PlantHydraulicModel::linearSystemMeunierSolve  matrix Compute with Eigen failed" );
        }

        try{
            v2= solverLU.solve(b);
        }catch(...){
             throw std::runtime_error("PlantHydraulicModel::linearSystemMeunierSolve error when solving wat. pot. xylem with Eigen ");
        }
    if(verbose)
		{
        std::cout << "Matrix norm: " << mat.norm()<<  "Matrix rows: " <<mat.rows()  <<" minDelta "<<minDelta << "\n";
        double normMat = mat.norm();            // Frobenius norm
        double normInvApprox = v2.norm() / b.norm(); // crude estimate of ||A^{-1}|| via solution
        double condEstimate = normMat * normInvApprox;
        std::cout << "Estimated condition number: " << condEstimate << "\n";
    }
    //Eigen::setNbThreads(1);
        psiXyl.assign(v2.data(), v2.data() + v2.size());
	if(solverLU.info() != Eigen::Success){
		std::cout << "PlantHydraulicModel::linearSystemSolve error when solving wat. pot. xylem with Eigen: " << solverLU.info() << std::endl;
		 throw std::runtime_error("Photosynthesis::linearSystemSolve error when solving wat. pot. xylem with Eigen ");
	}
    for (int si = 0; si<ms->segments.size(); si++) {//todo: remove afterwards maybe
        if (!std::isfinite(psiXyl.at(si)))
        {
            std::cout<<"PlantHydraulicModel::linearSystemSolve nan of Inf psiXyl.at(si): " <<si <<" "<< psiXyl.at(si) << std::endl;
             throw std::runtime_error("Photosynthesis::linearSystemSolve nan of Inf ");
        }
    }
}

void PlantHydraulicModel::linearSystemMeunier(double simTime, const std::vector<double> sx, bool cells,
			const std::vector<double> soil_k)
{
    linearSystemMeunier_(simTime, sx, cells, soil_k);
}
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
void PlantHydraulicModel::linearSystemMeunier_(double simTime, const std::vector<double>& sx, bool cells,
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
    minDelta = 1e6;
    size_t k=0;
    //# pragma omp for schedule(static)
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
            kr = params->kr_f_wrapped(si, age, subType, organType, cells);
        } catch(...) {
            std::cout << "\n PlantHydraulicModel::linearSystem: conductivities failed" << std::flush;
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
            std::cout << "PlantHydraulicModel::linearSystem: warning segment "<< si <<" length smaller 1.e-5 \n";
            //l = 1.e-5; // valid quick fix? (also in segFluxes)
        }
		double perimeter = ms->getPerimeter(si, l);//perimeter of exchange surface
        double vz = v.z / l; // normed direction

        double cii, cij, bi;
        double tau = -1;
        double exp_tau_l = -1;
        double exp_minus_tau_l = -1;
        double delta = -1;
        double idelta = -1;

        if ((perimeter * kr>1.e-16) && (l > 1e-5)) {
            tau = std::sqrt(perimeter * kr / kx); // Eqn (6)
            
            exp_tau_l = std::exp(tau * l);
            exp_minus_tau_l = std::exp(-tau * l);
            delta = exp_minus_tau_l - exp_tau_l;// Eqn (12)
            idelta = 1.0 / delta;
            if (minDelta > delta){minDelta = delta;}
            cii = -kx * idelta * tau * (exp_minus_tau_l + exp_tau_l); // Eqn (23)
            cij = 2 * kx * idelta * tau;  // Eqn 24
            bi = kx * vz; //  # Eqn 25
            if (!std::isfinite(bi) || !std::isfinite(cij) || !std::isfinite(cii) ) {
            	std::cout << "PlantHydraulicModel::linearSystemMeunier_: nan or Inf bi, cii, or cij. segIdx "<<si<<" organType "<<organType<<" subType "<<subType;
				std::cout <<" bi " << bi << ", cii " << cii << ", cij "<< cij <<" perimeter "<<perimeter<<" kr "<<kr<<" kx "<<kx;
				std::cout<< ", delta "<<delta << ", tau "<< tau << ", psi_s " << psi_s << ", exp_tau_l " << exp_tau_l
                    << ", exp_minus_tau_l " << exp_minus_tau_l<< "\n";
				throw std::runtime_error("PlantHydraulicModel::linearSystemMeunier_: nan or Inf bi, cii, or cij");
			}
            
        } else { // solution for a=0, or kr = 0
            cii = kx/l;
            cij = -kx/l;
            bi = kx * vz;
            psi_s = 0;//
            if (!std::isfinite(bi) || !std::isfinite(cij) || !std::isfinite(cii) ) {
            	std::cout << "PlantHydraulicModel::linearSystemMeunier_: nan or Inf bi, cii, or cij for small seg. segIdx "<<si<<" organType "<<organType<<" subType "<<subType << ", psi_s " << psi_s ;
				std::cout <<" bi " << bi << ", cii " << cii << ", cij "<< cij <<" perimeter "<<perimeter<<" kr "<<kr<<" kx "<<kx<< "\n";
				throw std::runtime_error("PlantHydraulicModel::linearSystemMeunier_: nan or Inf bi, cii, or cij");
			}
        }
        if(dovector)
        {
            k = fillVectors(k, i, j, bi, cii, cij, psi_s);
        }else{
            //k = si * 4;
         try{
            k = fillTripletList(k, i, j, bi, cii, cij, psi_s);
        }catch(...){
             std::cout << "PlantHydraulicModel::linearSystemMeunier_ error when fillTripletList "<<si<<" organType "<<organType<<" subType "<<subType;
				std::cout <<" bi " << bi << ", cii " << cii << ", cij "<< cij <<" perimeter "<<perimeter<<" kr "<<kr<<" kx "<<kx;
				std::cout<< ", delta "<<delta << ", tau "<< tau << ", psi_s " << psi_s << ", exp_tau_l " << exp_tau_l
                    << ", exp_minus_tau_l " << exp_minus_tau_l << " minDelta " <<minDelta << " l "<<l <<" perimeter * kr " <<perimeter * kr<< "\n";
             throw std::runtime_error("PlantHydraulicModel::linearSystemMeunier_ error when fillTripletList");
        }
            
        }
    }/*end omp parallel*/ 
}


/**
 * fill the matrices to be solved. overloads @see XylemFlux::fillVectors
 * @param k				index for the row- and column-index vectors
 * @param i, j			indexes of the non-zero elements of the sparse matrix
 * @param psi_s 		outer water potential [cm]
 * @param bi			value of variable b at row i [cm3/d]
 * @param cii			value of variable c at row i col i [cm2/d]
 * @param cij			value of variable c at row i col j [cm2/d]
 * @return k			next index for the row- and column-index vectors
 */
size_t PlantHydraulicModel::fillTripletList(size_t k, int i, int j, double bi, double cii, double cij, double psi_s) 
{
	typedef Eigen::Triplet<double> Tri;
	aB[i] += ( bi + cii * psi_s +cij * psi_s) ;

	//b(i) = aB[i];
	tripletList.push_back(Tri(i,i,cii));
	k += 1;
	tripletList.push_back(Tri(i,j,cij));
	k += 1;

	int ii = i;
	i = j;  j = ii; // edge ji
	aB[i] += ( -bi + cii * psi_s +cij * psi_s) ; // (-bi) Eqn (14) with changed sign

	//b(i) = aB[i];

            if (!std::isfinite(aB[i])) {
            	std::cout << "PlantHydraulicModel::fillTripletList: nan or Inf 2nd b(i). i "<<i<< " aB.size() " << aB.size() <<" k "<<k<<" j "<<j
                   <<" aB[i] "<<aB[i] <<" bi "<<bi<<" cii "<<cii;
				std::cout <<" psi_s " << psi_s << ", cij " << cij<<" cii * psi_s "<< cii * psi_s << " cij * psi_s "
                    << cij * psi_s<<" tot "<< -bi + cii * psi_s +cij * psi_s << "\n";
				throw std::runtime_error("PlantHydraulicModel::segFluxes: nan or Inf 2nd b(i)");
			}
    
	tripletList.push_back(Tri(i,i,cii));
	k += 1;

	tripletList.push_back(Tri(i,j,cij));
	k += 1;
	return k;
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
            kr = params->kr_f_wrapped(si, age, subType, organType, cells);
        } catch(...) {
            std::cout << "\n PlantHydraulicModel::segFluxes: conductivities failed" << std::flush;
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
            std::cout << "PlantHydraulicModel::segFluxes: warning segment length smaller 1.e-5 \n";
            //l = 1.e-5; // valid quick fix? (also in segFluxes)
        }

		double perimeter = ms->getPerimeter(si, l);//perimeter of exchange surface


        if ((perimeter * kr>1.e-16)&& (l > 1e-5)) { // only relevant for exact solution
            double f = -perimeter*kr; // flux is proportional to f // *rho*g
            double fApprox = f*l*(psi_s - rx[j]); // cm3 / day

            double tau = std::sqrt(perimeter*kr/kx); // sqrt(c) [cm-1]
            double d = std::exp(-tau*l)-std::exp(tau*l); // det
            double fExact = -f*(1./(tau*d))*(rx[i]-psi_s+rx[j]-psi_s)*(2.-std::exp(-tau*l)-std::exp(tau*l));
            if (!std::isfinite(fExact)) {
            	std::cout << "PlantHydraulicModel::segFluxes: nan or Inf fExact. segIdx "<<si<<" organType "<<organType<<" subType "<<subType;
				std::cout <<" tau " << tau << ", l " << l <<" perimeter "<<perimeter<<" kr "<<kr<<" kx "<<kx;
				std::cout<< ", d "<<d << ", rx "<< rx[i] << ", psi_s " << psi_s << ", f " << f << "\n";
				throw std::runtime_error("PlantHydraulicModel::segFluxes: nan or Inf fExact");
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
    if((psi_s >0) || (!std::isfinite(psi_s)))
    {
        std::cout << "PlantHydraulicModel::getPsiOut: nan or Inf or > 0 psi_s. segIdx "<<si<<" organType "<<organType;
				std::cout <<" psi_s "<<psi_s<<" cells " << cells ;
        if (cells) { 
				std::cout<< ", cellIndex "<<ms->seg2cell.at(si) << ", params->psi_air "<< params->psi_air ;
        }else{std::cout <<" sx_.at(si) " << sx_.at(si) << " " <<sx_.at(0);}
        std::cout << "\n";
				throw std::runtime_error("PlantHydraulicModel::getPsiOut: nan or Inf or > 0 psi_s");
    }
	return psi_s;
}


} // namespace
