///**

namespace CPlantBox {

// * Example

namespace CPlantBox {

// *

namespace CPlantBox {

// * defines a minimal organ from the base organ class

namespace CPlantBox {

// *

namespace CPlantBox {

// */

namespace CPlantBox {

//

namespace CPlantBox {

//

namespace CPlantBox {

///**

namespace CPlantBox {

// * example organ, having two nodes, from (0,0,0) and (0,0,1)

namespace CPlantBox {

// */

namespace CPlantBox {

//class MyOrgan : public Organ {

namespace CPlantBox {

//

namespace CPlantBox {

//public:

namespace CPlantBox {

//

namespace CPlantBox {

//	MyOrgan(Plant* plant, Organ* parent, int subtype, double delay): Organ(plant, parent, subtype, delay) {

namespace CPlantBox {

//		r_nodes.push_back(Vector3d(0,0,0));

namespace CPlantBox {

//	    nodeIDs.push_back(0);

namespace CPlantBox {

//	    nctimes.push_back(0.);

namespace CPlantBox {

//		r_nodes.push_back(Vector3d(0,0,1));

namespace CPlantBox {

//	    nodeIDs.push_back(0);

namespace CPlantBox {

//	    nctimes.push_back(0.);

namespace CPlantBox {

//	}

namespace CPlantBox {

//

namespace CPlantBox {

//	virtual void setRelativeOrigin(const Vector3d& o) override { this->o = o; }

namespace CPlantBox {

//    virtual Vector3d getRelativeOrigin() const override { return o;  }

namespace CPlantBox {

//	virtual void setRelativeHeading(const Matrix3d& m) override { this->A = m; }

namespace CPlantBox {

//    virtual Matrix3d getRelativeHeading() const override { return A; }

namespace CPlantBox {

//

namespace CPlantBox {

//	Vector3d o;

namespace CPlantBox {

//	Matrix3d A;

namespace CPlantBox {

//

namespace CPlantBox {

//};

namespace CPlantBox {

//

namespace CPlantBox {

//

namespace CPlantBox {

//

namespace CPlantBox {

//

namespace CPlantBox {

//void example_scenegraph()

namespace CPlantBox {

//{

namespace CPlantBox {

//	Plant* plant = new Plant(); // we use the plant only for the vtp export

namespace CPlantBox {

//

namespace CPlantBox {

//	Seed* seed = new Seed(plant);

namespace CPlantBox {

//	MyOrgan* base = new MyOrgan(plant, nullptr, 0, 0);

namespace CPlantBox {

//	MyOrgan* branch1 = new MyOrgan(plant, nullptr, 0, 0); // at top of base

namespace CPlantBox {

//	MyOrgan* branch2 = new MyOrgan(plant, nullptr, 0, 0); // half way base

namespace CPlantBox {

//	MyOrgan* branch3 = new MyOrgan(plant, nullptr, 0, 0); // at top of branch 1

namespace CPlantBox {

//

namespace CPlantBox {

//	seed->setRelativeOrigin(Vector3d(0.,0.,2));

namespace CPlantBox {

//

namespace CPlantBox {

//	base->setRelativeOrigin(Vector3d(0.,0.,0.));

namespace CPlantBox {

//	base->setRelativeHeading(Matrix3d::rotX(30./180.*3.1415)); // rotated 30 degrees

namespace CPlantBox {

//

namespace CPlantBox {

//	branch1->setRelativeOrigin(Vector3d(0.,0.,1.)); // at top

namespace CPlantBox {

//	branch1->setRelativeHeading(Matrix3d::rotX(30./180.*3.1415)); // rotated 30 degrees

namespace CPlantBox {

//

namespace CPlantBox {

//	branch2->setRelativeOrigin(Vector3d(0.,0.,0.5)); // half way

namespace CPlantBox {

//	branch2->setRelativeHeading(Matrix3d::rotX(30./180.*3.1415)); // X rotated 30 degrees

namespace CPlantBox {

//

namespace CPlantBox {

//	branch3->setRelativeOrigin(Vector3d(0.,0.,1.)); // at top

namespace CPlantBox {

//	branch3->setRelativeHeading(Matrix3d::rotY(30./180.*3.1415)); // Y rotated 30 degrees

namespace CPlantBox {

//

namespace CPlantBox {

//	// make tree by hand

namespace CPlantBox {

//	seed->children.push_back(base);

namespace CPlantBox {

//	base->children.push_back(branch1);

namespace CPlantBox {

//	base->children.push_back(branch2);

namespace CPlantBox {

//	branch1->children.push_back(branch3);

namespace CPlantBox {

//	base->parent = seed;

namespace CPlantBox {

//	branch1->parent = base;

namespace CPlantBox {

//	branch2->parent = base;

namespace CPlantBox {

//	branch3->parent = branch1;

namespace CPlantBox {

//

namespace CPlantBox {

//	plant->seed = seed;

namespace CPlantBox {

//	plant->write("example_scenegraph.vtp", Organ::ot_organ);

namespace CPlantBox {

//

namespace CPlantBox {

//}

namespace CPlantBox {

