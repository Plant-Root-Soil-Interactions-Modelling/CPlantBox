#include "Plant.h"
#include <memory>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>

namespace CPlantBox {

/*
 * Sets up the xml reader Organism::readParameters
 */
Plant::Plant(): Organism()
{
    initPrototypes(new SeedRandomParameter(this), new RootRandomParameter(this), new StemRandomParameter(this), new LeafRandomParameter(this));
}

/**
 * Copy Constructor
 *
 * deep copies the plant
 * does not deep copy geometry, elongation functions, and soil (all not owned by plant)
 *
 * @param p        plant that is copied
 */
Plant::Plant(const Plant& p): Organism(p), geometry(p.geometry), soil(p.soil)
{
    std::cout << "Copying plant";
}

/**
 *
 */
void Plant::initPrototypes(OrganRandomParameter* seed, OrganRandomParameter* root, OrganRandomParameter* stem,
    OrganRandomParameter* leaf)
{
    seed->subType = 0;
    setOrganRandomParameter(seed);
    root->subType = 0;
    setOrganRandomParameter(root);
    stem->subType = 0;
    setOrganRandomParameter(stem);
    leaf->subType = 0;
    setOrganRandomParameter(leaf);
}


/**
 * Resets the root system: deletes all roots, sets simulation time to 0.
 */
void Plant::reset()
{
	if (seed!=nullptr) {
		delete seed;
		seed = nullptr;
	}
	baseOrgans.clear();
    simtime = 0;
    organId = -1;
    nodeId = -1;
}

/**
 * Sets up the plant according to the plant parameters,
 * a confining geometry, the tropism functions, and the growth functions.
 *
 * Call this method before simulation and after setting geometry, plant and root parameters
 */
void Plant::initialize()
{
    reset(); // just in case

    // create seed
    Seed seed = Seed(this);
    seed.initialize();

    // todo check what this does
    oldNumberOfNodes = getNumberOfNodes();

    // further intializations
    initCallbacks();
}

/**
 * Called by RootSystem::initialize.
 * Sets up tropism and growth functions call backs using
 * RootSystem::createTropismFunction and RootSystem::createGrowthFunction
 */
void Plant::initCallbacks()
{
    // Create tropisms and growth functions per root type
    for (auto& p_otp :organParam[Organism::ot_root]) {
        RootRandomParameter* rtp = (RootRandomParameter*)p_otp.second;
        Tropism* tropism = this->createTropismFunction(rtp->tropismT, rtp->tropismN, rtp->tropismS);
        tropism->setGeometry(geometry);
        delete rtp->f_tf; // delete old tropism
        rtp->f_tf = tropism; // set new one
        GrowthFunction* gf_ = this->createGrowthFunction(rtp->gf);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (otherwise an exception is thrown)
        delete rtp->f_gf;
        rtp->f_gf  = gf_;
    }
    for (auto& p_otp :organParam[Organism::ot_leaf]) {
        LeafRandomParameter* rtp = (LeafRandomParameter*)p_otp.second;
        Tropism* tropism = this->createTropismFunction(rtp->tropismT, rtp->tropismN, rtp->tropismS);
        tropism->setGeometry(geometry);
        delete rtp->f_tf; // delete old tropism
        rtp->f_tf = tropism; // set new one
        GrowthFunction* gf_ = this->createGrowthFunction(rtp->gf);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (otherwise an exception is thrown)
        delete rtp->f_gf;
        rtp->f_gf  = gf_;
    }
    for (auto& p_otp :organParam[Organism::ot_stem]) {
        StemRandomParameter* rtp = (StemRandomParameter*)p_otp.second;
        Tropism* tropism = this->createTropismFunction(rtp->tropismT, rtp->tropismN, rtp->tropismS);
        tropism->setGeometry(geometry);
        delete rtp->f_tf; // delete old tropism
        rtp->f_tf = tropism; // set new one
        GrowthFunction* gf_ = this->createGrowthFunction(rtp->gf);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (otherwise an exception is thrown)
        delete rtp->f_gf;
        rtp->f_gf  = gf_;
    }

}

/**
 * Creates a specific tropism from the tropism type index.
 * the function must be extended or overwritten to add more tropisms.
 *
 * @param tt        the tropism type index, given in the root type parameters
 * @param N         tropism parameter (passsed to the tropism class)
 * @param sigma     tropism parameter (passsed to the tropism class)
 * @return          the tropism class containing with the callback functions
 */
Tropism* Plant::createTropismFunction(int tt, double N, double sigma) {
    // std::cout << "Creating (" << tt << ", " << N << ", " << sigma <<")\n";
    switch (tt) {
    case tt_plagio: return new Plagiotropism(this,N,sigma);
    case tt_gravi: return new Gravitropism(this,N,sigma);
    case tt_exo: return new Exotropism(this,N,sigma);
    case tt_hydro: {
        Tropism* gt =  new Gravitropism(this,N,sigma);
        Tropism* ht= new Hydrotropism(this,N,sigma,soil);
        Tropism* cht = new CombinedTropism(this,N,sigma,ht,10.,gt,1.); // does only use the objective functions from gravitropism and hydrotropism
        return cht;
    }
    case tt_twist:  return new TwistTropism(this,N,sigma);
    case tt_antigravi: return new AntiGravitropism(this,N,sigma);


    default: throw std::invalid_argument( "Plant::createTropismFunction() tropism type not implemented" );
    }
}

/**
 * Creates a growth functions from the growth function index.
 * This function must bee extended or overwritten to add more growth functions
 *
 * @param gft       the growth function index, given in the root type parameters
 * @return          the growth function class containing with the callback functions
 */
GrowthFunction* Plant::createGrowthFunction(int gft) {
    switch (gft) {
    case gft_negexp: return new ExponentialGrowth();
    case gft_linear: return new LinearGrowth();
    default: throw std::invalid_argument( "Plant::createGrowthFunction() growth function type not implemented" );
    }
}

/**
 * Simulates plant growth for the time span defined in the root system parameters
 */
void Plant::simulate()
{
	auto srp = (SeedRandomParameter*)organParam[Organism::ot_seed][0];
	Organism::simulate(srp->simtime);
}

/**
 * todo most important debug informations
 */
std::string Plant::toString() const
{
	return "todo";
}

/**
write VTP using tinyXML
 **/
void Plant::TiXMLwriteVTP(int otype, std::ostream & os) const // Write .VTP file by using TinyXML2 performance slowed by 0.5 seconds but precision increased
{
	tinyxml2::XMLPrinter printer( 0, false, 0 );

	auto organs = this->getOrgans(otype); // update roots (if necessary)
	auto nodes = getPolylines(otype);
	auto times = getPolylineCTs(otype);

	os << "<?xml version=\"1.0\"?>";
	printer.OpenElement("VTKFile"); printer.PushAttribute("type", "PolyData"); printer.PushAttribute("version", "0.1"); printer.PushAttribute("byte_order", "LittleEndian");
	printer.OpenElement("PolyData");
	int non = 0; // number of nodes
    for (const auto& r : organs) {
		non += r->getNumberOfNodes();
	}
	int nol=organs.size(); // number of lines
	printer.OpenElement("Piece"); printer.PushAttribute("NumberOfLines",  nol); printer.PushAttribute("NumberOfPoints", non);

	// POINTDATA
	printer.OpenElement("PointData"); printer.PushAttribute("Scalars", "Pointdata");
	printer.OpenElement("DataArray"); printer.PushAttribute("type", "Float32");  printer.PushAttribute("Name", "time"); printer.PushAttribute("NumberOfComponents", "1"); printer.PushAttribute("format", "ascii" );
	for (std::vector<double> r: times) {
		for (double t : r) {
			printer.PushText(t); printer.PushText(" ");
		}
	}
	printer.CloseElement();
	printer.CloseElement();

	// CELLDATA (live on the polylines)
	printer.OpenElement("CellData"); printer.PushAttribute("Scalars", "CellData" );
	std::vector<std::string> sTypeNames = { "organtype", "id",  "emergencetime", "creationtime", "age", "subtype", "order", "radius"}; //  , "order", "radius", "subtype" ,
	for (size_t i=0; i<sTypeNames.size(); i++) {
		std::string sType = sTypeNames[i];
        const char *schar = sType.c_str();
		printer.OpenElement("DataArray"); printer.PushAttribute("type", "Float32");  printer.PushAttribute("Name", schar); printer.PushAttribute("NumberOfComponents", "1"); printer.PushAttribute("format", "ascii" );
		std::vector<double> scalars = getParameter(sTypeNames[i], otype);
		for (double s : scalars) {
			printer.PushText(s); printer.PushText(" ");
		}
		printer.CloseElement();
	}
	printer.CloseElement();

	// POINTS (=nodes)
	printer.OpenElement("Points");
	printer.OpenElement("DataArray"); printer.PushAttribute("type", "Float32");  printer.PushAttribute("Name", "Coordinates"); printer.PushAttribute("NumberOfComponents", "3"); printer.PushAttribute("format", "ascii" );
    for (const auto& r : nodes) {
        for (const auto& n : r) {
			printer.PushText(n.x); printer.PushText(" "); printer.PushText(n.y); printer.PushText(" "); printer.PushText(n.z); printer.PushText(" ");
		}
	}
	printer.CloseElement();
	printer.CloseElement();

	// LINES (polylines)
	printer.OpenElement("Lines");
	printer.OpenElement("DataArray"); printer.PushAttribute("type", "Float32");  printer.PushAttribute("Name", "connectivity"); printer.PushAttribute("NumberOfComponents", "1"); printer.PushAttribute("format", "ascii" );
	int c=0;
    for (const auto& r : organs) {
		for (size_t i=0; i<r->getNumberOfNodes(); i++) {
			printer.PushText(c); printer.PushText(" ");
			c++;
		}
	}
	printer.CloseElement();

	printer.OpenElement("DataArray"); printer.PushAttribute("type", "Float32");  printer.PushAttribute("Name", "offsets"); printer.PushAttribute("NumberOfComponents", "1"); printer.PushAttribute("format", "ascii" );
	c = 0;
    for (const auto& r : organs) {
		c += r->getNumberOfNodes();
		printer.PushText(c); printer.PushText(" ");
	}
	printer.CloseElement();
	printer.CloseElement();
	printer.CloseElement();
	printer.CloseElement();
	printer.CloseElement();
	os << std::string(printer.CStr());
}





} // namespace CPlantBox
