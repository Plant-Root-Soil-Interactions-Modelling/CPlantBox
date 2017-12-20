/**
 * Example
 *
 * defines a minimal organ from the base organ class
 *
 */


/**
 * example organ, having two nodes, from (0,0,0) and (0,0,1)
 */
class MyOrgan : public Organ {

public:

	MyOrgan(Plant* plant, Organ* parent, int subtype, double delay): Organ(plant, parent, subtype, delay) {
		r_nodes.push_back(Vector3d(0,0,0));
	    nodeIDs.push_back(0);
	    nctimes.push_back(0.);
		r_nodes.push_back(Vector3d(0,0,1));
	    nodeIDs.push_back(0);
	    nctimes.push_back(0.);
	}

	virtual void setRelativeOrigin(const Vector3d& o) override { this->o = o; }
    virtual Vector3d getRelativeOrigin() const override { return o;  }
	virtual void setRelativeHeading(const Matrix3d& m) override { this->A = m; }
    virtual Matrix3d getRelativeHeading() const override { return A; }

	Vector3d o;
	Matrix3d A;
};




void example_scenegraph()
{
	Plant* plant = new Plant(); // we use the plant only for the vtp export

	Seed* seed = new Seed(plant);
	MyOrgan* base = new MyOrgan(plant, nullptr, 0, 0);
	MyOrgan* branch1 = new MyOrgan(plant, nullptr, 0, 0); // at top of base
	MyOrgan* branch2 = new MyOrgan(plant, nullptr, 0, 0); // half way base
	MyOrgan* branch3 = new MyOrgan(plant, nullptr, 0, 0); // at top of branch 1

	base->setRelativeOrigin(Vector3d(0.,0.,0.));
	base->setRelativeHeading(Matrix3d());

	branch1->setRelativeOrigin(Vector3d(0.,0.,1.)); // at top
	branch1->setRelativeHeading(Matrix3d()); // rotated 30 degrees

	branch2->setRelativeOrigin(Vector3d(0.,0.,0.5)); // half way
	branch2->setRelativeHeading(Matrix3d::rotX(30./180.*3.1415)); // X rotated 30 degrees

	branch3->setRelativeOrigin(Vector3d(0.,0.,1.)); // at top
	branch3->setRelativeHeading(Matrix3d::rotY(30./180.*3.1415)); // Y rotated 30 degrees

	// make tree by hand
	seed->children.push_back(base);
	base->children.push_back(branch1);
	base->children.push_back(branch2);
	branch1->children.push_back(branch3);
	base->parent = seed;
	branch1->parent = base;
	branch2->parent = base;
	branch3->parent = branch1;

	plant->seed = seed;
	plant->write("example_scenegraph.vtp", Organ::ot_organ);

}
