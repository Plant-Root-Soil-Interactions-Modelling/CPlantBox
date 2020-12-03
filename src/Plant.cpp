#include "Plant.h"

#include <memory>
#include <iostream>
#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include <fstream>

namespace CPlantBox {

/*
 * todo doc me
 */
Plant::Plant(): Organism()
{ }

/**
 * todo doc me
 */
std::shared_ptr<Organism> Plant::copy()
{
    auto no = std::make_shared<Plant>(*this); // copy constructor
    for (int i=0; i<baseOrgans.size(); i++) {
        no->baseOrgans[i] = baseOrgans[i]->copy(no);
    }
    for (int ot = 0; ot < numberOfOrganTypes; ot++) { // copy organ type parameters
        for (auto& otp : organParam[ot]) {
            otp.second = otp.second->copy(no);
        }
    }
    return no;
}

/**
 * Returns the seed of the plant
 */
std::shared_ptr<Seed> Plant::getSeed()
{
	return std::static_pointer_cast<Seed>(baseOrgans.at(0));
}

/**
 * todo docme , this could be made unique? and probably should be protected
 */
void Plant::initializeReader()
{
    auto rrp = std::make_shared<RootRandomParameter>(shared_from_this());
    rrp->subType = 0;
    setOrganRandomParameter(rrp);
    auto srp = std::make_shared<SeedRandomParameter>(shared_from_this());
    srp->subType = 0;
    setOrganRandomParameter(srp);
    auto strp = std::make_shared<StemRandomParameter>(shared_from_this());
    strp->subType = 0;
    setOrganRandomParameter(strp);
    auto lrp = std::make_shared<LeafRandomParameter>(shared_from_this());
    lrp->subType = 0;
    setOrganRandomParameter(lrp);
}

/**
 * Resets the root system: deletes all roots, sets simulation time to 0.
 */
void Plant::reset()
{
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
void Plant::initialize(bool verbose)
{
    reset(); // just in case

    // create seed
    auto seed = std::make_shared<Seed>(shared_from_this());
    seed->initialize(verbose);
    baseOrgans.push_back(seed);

    oldNumberOfNodes = getNumberOfNodes(); // todo check what this does

    // further initializations
    initCallbacks();
}

/**
 * Called by RootSystem::initialize.
 * Sets up tropism and growth functions call backs using
 * RootSystem::createTropismFunction and RootSystem::createGrowthFunction
 */
void Plant::initCallbacks()
{
    // Create tropisms and growth functions per random root parameter
    for (auto& p_otp :organParam[Organism::ot_root]) {
        auto rp = std::static_pointer_cast<RootRandomParameter>(p_otp.second);
        auto tropism = this->createTropismFunction(rp->tropismT, rp->tropismN, rp->tropismS);
        tropism->setGeometry(geometry);
        rp->f_tf = tropism; // set new one
        auto gf_ = this->createGrowthFunction(rp->gf);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (otherwise an exception is thrown)
        rp->f_gf  = gf_;
    }
    // Create tropisms and growth functions per random leaf parameter
    for (auto& p_otp :organParam[Organism::ot_leaf]) {
        auto rp = std::static_pointer_cast<LeafRandomParameter>(p_otp.second);
        auto tropism = this->createTropismFunction(rp->tropismT, rp->tropismN, rp->tropismS);
        tropism->setGeometry(geometry);
        rp->f_tf = tropism; // set new one
        auto gf_ = this->createGrowthFunction(rp->gf);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (otherwise an exception is thrown)
        rp->f_gf  = gf_;
    }
    // Create tropisms and growth functions per random stem parameter
    for (auto& p_otp :organParam[Organism::ot_stem]) {
        auto rp = std::static_pointer_cast<StemRandomParameter>(p_otp.second);
        auto tropism = this->createTropismFunction(rp->tropismT, rp->tropismN, rp->tropismS);
        tropism->setGeometry(geometry);
        rp->f_tf = tropism; // set new one
        auto gf_ = this->createGrowthFunction(rp->gf);
        gf_->getAge(1,1,1,nullptr);  // check if getAge is implemented (otherwise an exception is thrown)
        rp->f_gf  = gf_;
    }

}

/**
 * Manually sets a tropism function for a specific or for all root types. TODO
 * Must be called after RootSystem::initialize(), otherwise its overwritten.
 *
 * @param tf_           a tropism
 * @param rt            root type, if rt = -1 all types are set to this tropism (default).
 */
void Plant::setTropism(std::shared_ptr<Tropism> tf, int organType, int subType) // todo
{
    if (organType==-1) { // call for all relevant
        setTropism(tf, Organism::ot_root, subType);
        setTropism(tf, Organism::ot_stem, subType);
        setTropism(tf, Organism::ot_leaf, subType);
        return;
    }
    // put it into a vector
    std::vector<std::shared_ptr<OrganRandomParameter>> orp;
    if (subType>-1) { // set for a specific root type
        orp.push_back(getOrganRandomParameter(organType, subType));
    } else { // set for all root types (default)
        for (auto& orp_ :organParam[organType]) {
            orp.push_back(orp_.second);
        }
    }
    // static cast since tropisms are not in the base class
    for (auto& orp_ : orp)
    switch(organType) {
    case Organism::ot_root: {
        auto rtp = std::static_pointer_cast<RootRandomParameter>(orp_);
        rtp->f_tf = tf;
    } break;
    case Organism::ot_stem: {
        auto rtp = std::static_pointer_cast<StemRandomParameter>(orp_);
        rtp->f_tf = tf;
    } break;
    case Organism::ot_leaf: {
        auto rtp = std::static_pointer_cast<LeafRandomParameter>(orp_);
        rtp->f_tf = tf;
    } break;
    default: throw std::invalid_argument( "Plant::setTropism() cannot set tropism for organ type "+ Organism::organTypeName(organType) );
    }
}

/**
 * Simulates plant growth for the time span defined in the root system parameters
 */
	void Plant::simulate(double dt, bool verbose)	
{	
		
    Organism::simulate(dt, verbose);	
}

/**
 * Simulates plant growth for the time span defined in the root system parameters
 */
void Plant::simulate()
{
    auto srp = std::static_pointer_cast<SeedRandomParameter>(organParam[Organism::ot_seed][0]);
    Organism::simulate(srp->simtime);
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
std::shared_ptr<Tropism> Plant::createTropismFunction(int tt, double N, double sigma) {
    // std::cout << "Creating (" << tt << ", " << N << ", " << sigma <<")\n";
    switch (tt) {
    case tt_plagio: return std::make_shared<Plagiotropism>(shared_from_this(),N,sigma);
    case tt_gravi: return std::make_shared<Gravitropism>(shared_from_this(),N,sigma);
    case tt_exo: return std::make_shared<Exotropism>(shared_from_this(),N,sigma);
    case tt_hydro: { // uses weighted objective functions from gravitropism and hydrotropism
        auto gt =  std::make_shared<Gravitropism>(shared_from_this(),N,sigma);
        auto ht= std::make_shared<Hydrotropism>(shared_from_this(),N, sigma, soil);
        return std::make_shared<CombinedTropism>(shared_from_this(),N,sigma,ht,10.,gt,1.);
    }
    case tt_twist:  return std::make_shared<TwistTropism>(shared_from_this(),N,sigma);
    case tt_antigravi: return std::make_shared<AntiGravitropism>(shared_from_this(),N,sigma);
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
std::shared_ptr<GrowthFunction>Plant::createGrowthFunction(int gft) {
    switch (gft) {
    case gft_negexp: return std::make_shared<ExponentialGrowth>();
    case gft_linear: return std::make_shared<LinearGrowth>();
    default: throw std::invalid_argument( "Plant::createGrowthFunction() growth function type not implemented" );
    }
}

/**
 * todo most important debug informations
 */
std::string Plant::toString() const
{
    return "Plant " + Organism::toString();
}

/**
 * Exports the simulation results with the type from the extension in name
 * (that must be lower case)
 *
 * todo move to Organism
 *
 * @param name      file name e.g. output.vtp
 */
void Plant::write(std::string name) const
{
    std::string ext = name.substr(name.size()-3,name.size()); // pick the right writer
    if (ext.compare("sml")==0) {
        std::cout << "writing RSML... "<< name.c_str() <<"\n";
        writeRSML(name); // use base class writer
    } else if (ext.compare("vtp")==0) {
        std::cout << "writing VTP... "<< name.c_str() <<"\n";
        std::ofstream fos;
        fos.open(name.c_str());
        writeVTP(-1, fos);
        fos.close();
//    } else if (ext.compare(".py")==0)  {
//        std::cout << "writing Geometry ... "<< name.c_str() <<"\n";
//        writeGeometry(fos);
    } else {
        throw std::invalid_argument("RootSystem::write(): Unkwown file type");
    }
}

/**
write VTP using tinyXML
todo move to Organism
 **/
void Plant::writeVTP(int otype, std::ostream & os) const // Write .VTP file by using TinyXML2 performance slowed by 0.5 seconds but precision increased
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
    std::vector<std::string> sTypeNames = { "organType", "id", "creationTime", "age", "subType", "order", "radius"}; //  , "order", "radius", "subtype" ,
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
