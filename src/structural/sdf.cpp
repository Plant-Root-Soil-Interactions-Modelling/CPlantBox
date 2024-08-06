// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
#include "sdf.h"

namespace CPlantBox {

std::string SignedDistanceFunction::writePVPScript() const
{
    std::stringstream str;
    this->writePVPScript(str,1);
    return str.str();
}



/**
 * Returns the signed distance to the next boundary of the box
 *
 * @param v     spatial position [cm]
 * \return      signed distance [cm], a minus sign means inside, plus outside
 */
double SDF_PlantBox::getDist(const Vector3d& v) const
{
    double z = v.z+dim.z; //  translate
    return -std::min(std::min(std::min(std::min(std::min(dim.z+z,dim.z-z),dim.y+v.y),dim.y-v.y),dim.x+v.x),dim.x-v.x);
}

/**
 * Writes a ParaView Phython script explicitly representing the implicit geometry
 *
 * @param cout      e.g. a file output stream
 * @param c         python object counter for the script (to avoid duplicate names)
 * \return          object counter
 */
int SDF_PlantBox::writePVPScript(std::ostream & cout, int c) const
{
    std::string name = "obj";
    name.append(std::to_string(c));
    cout << name <<" = Box()\n\n" <<
        name << ".XLength = " << 2.*dim.x << "\n" << name << ".YLength = " << 2.*dim.y << "\n" << name << ".ZLength = "<< 2.*dim.z << "\n" <<
        name << ".Center = [0.,0., "<< -dim.z << "]\n\n";
    cout << name << "Display = Show(" << name << ",renderView1)\n" << name << "Display.Opacity = 0.1\n" <<
        name<< "Display.DiffuseColor = [0., 0., 1.0]\n" << "renderView1.ResetCamera()\n";
    // created 1 object and corresponding Display
    c++;
    return c;
}



/**
 * Returns the signed distance to the next boundary of the cuboid
 *
 * @param v     spatial position [cm]
 * \return      signed distance [cm], a minus sign means inside, plus outside
 */
double SDF_Cuboid::getDist(const Vector3d& v) const
{
    // d=-min(min(min(min(min(-z1+p(:,3),z2-p(:,3)),-y1+p(:,2)),y2-p(:,2)),-x1+p(:,1)),x2-p(:,1));
    double d =   std::min(   -min.z+v.z,  max.z-v.z);
    d = std::min(std::min(d, -min.y+v.y), max.y-v.y);
    d = std::min(std::min(d, -min.x+v.x), max.x-v.x);
    return -d;
}



/**
 * Creates a cylindrical or square container
 *
 * @param r1_   top radius [cm]
 * @param r2_   bottom radius [cm]
 * @param h_    height of the container [cm]
 * @param sq    square (true) or cylindrical (false), default=false
 */
SDF_PlantContainer::SDF_PlantContainer(double r1_, double r2_, double h_, double sq) : r1(r1_), r2(r2_), h(h_), square(sq)
{ }

/**
 * Returns the signed distance to the next boundary
 *
 * @param v     spatial position [cm]
 * \return      signed distance [cm], a minus sign means inside, plus outside
 */
double SDF_PlantContainer::getDist(const Vector3d& v) const
{
    double z = v.z/h; // 0 .. -1
    double r =  (1+z)*r1 - z*r2;
    double d;
    if (square) { // rectangular pot
        d = std::max(std::abs(v.x),std::abs(v.y))-r;
    } else { // round pot
        d = sqrt(v.x*v.x+v.y*v.y)-r;
    }
    return std::max(d,-std::min(h+v.z,0.-v.z));
}

/**
 * Writes a ParaView Phython script explicitly representing the implicit geometry
 *
 * @param cout      e.g. a file output stream
 * @param c         python object counter for the script (to avoid duplicate names)
 * \return          object counter
 */
int SDF_PlantContainer::writePVPScript(std::ostream & cout, int c) const
{
    if (r1!=r2) {
        throw std::invalid_argument( "PlantContainer::writePVPScript() not implemented yet, only the case r1==r2 is implented" );
        // TODO could be done by a cone with a clipping plane
    }

    std::string name = "obj";
    name.append(std::to_string(c));

    if (square) {
        cout << name <<" = Box()\n\n" <<
            name << ".XLength = " << 2.*r1 << "\n" << name << ".YLength = " << 2.*r1 << "\n" << name << ".ZLength = "<< h << "\n" <<
            name << ".Center = [0.,0., "<< -h/2. << "]\n\n";
        cout << name << "Display = Show(" << name << ",renderView1)\n";
        c++; // created 1 object and corresponding Display
    } else {
        cout << name <<" = Cylinder()\n\n" <<
            name << ".Resolution = 50\n" <<
            name << ".Height = " << h << "\n" <<
            name << ".Radius = " << r1 << "\n" <<
            name << ".Center = [0., "<< -h/2. << ",0.]\n\n"; // y will be z after the rotation...
        c++;
        std::string name2 = "obj";
        name2.append(std::to_string(c));
        cout << name2 << " = Transform(Input="<< name <<")\n" <<
            name2 << ".Transform = 'Transform'\n" << // (?)
            name2 << ".Transform.Rotate = [90.0, 0.0, 0.0]\n\n";
        cout << name << "Display = Show(" << name2 << ",renderView1)\n";
        c++; // created 2 objects and 1 Display
    }
    cout << name << "Display.Opacity = 0.2\n" << name << "Display.DiffuseColor = [0., 0., 1.0]\n" << "renderView1.ResetCamera()\n";

    return c;
}



/**
 * Rotates the base geometry pot around a main axis and then translates it
 *
 * @param sdf_          points to the original geometry
 * @param a             rotation angle [degree]
 * @param axis_         rotation axis (SDF_Axes: sdf_xaxis, sdf_yaxis, or sdf_zaxis)
 * @param pos_          position of origin [cm], default=(0,0,0)
 */
SDF_RotateTranslate::SDF_RotateTranslate(std::shared_ptr<SignedDistanceFunction> sdf_, double a, int axis_, const Vector3d& pos_)
{
    sdf=sdf_;
    pos=pos_;
    axis =axis_;  // remember for the python script
    angle = a;
    a = a/180.*M_PI; // convert to radian
    assert(((axis>=0) && (axis<3)));
    switch (axis) {
    case xaxis:
        A = Matrix3d::rotX(-a); break;
    case yaxis:
        A = Matrix3d::rotY(-a); break;
    case zaxis:
        A = Matrix3d::rotZ(-a); break;
    default:
        throw std::exception();
    }
}

/**
 * Distance to the next boundary, after rotation and translation of the base geometry
 *
 * @param v    spatial position [cm]
 * \return      signed distance [cm], a minus sign6 means inside, plus outside
 */
double SDF_RotateTranslate::getDist(const Vector3d& v) const
{
    Vector3d p = (A.times(v.minus(pos))); // p= rot(a)v+pos, -> v = rot(-a) (p-pos)
    return sdf->getDist(p);
}

/**
 * Writes a ParaView Phython script explicitly representing the implicit geometry
 *
 * @param cout      e.g. a file output stream
 * @param c         python object counter for the script (to avoid duplicate names)
 * \return          object counter
 */
int SDF_RotateTranslate::writePVPScript(std::ostream & cout, int c) const
{
    c =  sdf->writePVPScript(cout,c); // write script for base geometry

    std::string base = "obj";
    base.append(std::to_string(c-1)); // name of phython object of the base geometry
    cout << "\nHide(" << base <<")\n";

    std::string name = "obj";
    name.append(std::to_string(c));
    cout << name << "= Transform(Input=" << base << ")\n" << name << ".Transform = 'Transform'\n" <<
        name <<".Transform.Translate = ["<< pos.x<<","<<pos.y<<","<<pos.z<<"]\n";
    cout << name <<".Transform.Rotate = [";
    switch (axis) {
    case xaxis:
        cout << angle << ",0,0"; break;
    case yaxis:
        cout << "0,"<< angle << ",0"; break;
    case zaxis:
        cout << "0,0,"<< angle; break;
    }
    cout << "]\n\n";
    cout << name << "Display = Show(" << name << ",renderView1)\n" << name << "Display.Opacity = 0.1\n" <<
        name<< "Display.DiffuseColor = [0., 0., 1.0]\n" << "renderView1.ResetCamera()\n";
    c++;
    return c;
}



/**
 * Computes the intersection of the original geometries
 *
 * @param v     spatial position [cm]
 * \return      signed distance [cm], a minus sign means inside, plus outside
 */
double SDF_Intersection::getDist(const Vector3d& v) const {
    double d = sdfs[0]->getDist(v);
    for (size_t i=1; i<sdfs.size(); i++) {
        d = std::max(d, sdfs[i]->getDist(v));
    }
    return d;
}

/**
 * Writes a ParaView Phython script explicitly representing the implicit geometry
 *
 * Intersection/Union/Difference is not calculated,
 * but all original geometry is represented and grouped together for vizualisation
 *
 * @param cout      e.g. a file output stream
 * @param c         python object counter for the script (to avoid duplicate names)
 * \return          object counter
 */
int SDF_Intersection::writePVPScript(std::ostream & cout, int c) const
{
    cout << "# GROUP \n";

    // write all
    std::vector<int> goi; // group object indices
    for (const auto& sdf : sdfs) {
        c = sdf->writePVPScript(cout, c);
        goi.push_back(c-1);
        cout << "#\n";
    }

    // TODO hide all

    // group all
    std::string name = "obj";
    name.append(std::to_string(c));
    cout << "\n"<< name << "= GroupDatasets(Input=[";
    cout << "obj" << goi.at(0);
    for (size_t i = 1; i<goi.size(); i++) {
        cout << "," << "obj"<< goi.at(i);
    }
    cout <<"]) \n";

    cout << name << "Display = Show(" << name <<", renderView1)\n";
    //cout << "ColorBy("<< name <<"Display, None)\n";
    cout << name << "Display.Opacity = 0.1\n" << name << "Display.DiffuseColor = [0., 0., 1.0]\n" << "renderView1.ResetCamera()\n";

    c++;
    return c;
}



/**
 * Computes the union of the original geometries
 *
 * @param v     spatial position [cm]
 * \return      signed distance [cm], a minus sign means inside, plus outside
 */
double SDF_Union::getDist(const Vector3d& v) const
{
    double d = sdfs[0]->getDist(v);
    for (size_t i=1; i<sdfs.size(); i++) {
        d = std::min(d, sdfs[i]->getDist(v));
    }
    return d;
}



/**
 * Computes the difference of the original geometries
 *
 * @param v     spatial position [cm]
 * \return      signed distance [cm], a minus sign means inside, plus outside
 */
double SDF_Difference::getDist(const Vector3d& v) const
{
    double d = sdfs[0]->getDist(v);
    for (size_t i=1; i<sdfs.size(); i++) {
        d = std::max(d, -sdfs[i]->getDist(v));
    }
    return d;
}



/**
 * Constructor
 */
SDF_HalfPlane::SDF_HalfPlane(const Vector3d& o, const Vector3d& n_) : o(o),n(n_)
{
    n.normalize();
    Matrix3d A = Matrix3d::ons(n);
    p1 = A.column(1);
    p2 = A.column(2);
    p1.times(10.); // rescale to 10 cm
    p2.times(10.);
}


SDF_HalfPlane::SDF_HalfPlane(const Vector3d& o, const Vector3d& p1, const Vector3d& p2): o(o), p1(p1), p2(p2)
{
    Vector3d v1 = p1.minus(o);
    Vector3d v2 = p2.minus(o);
    n = v1.cross(v2);
    n.normalize();
    //	std::cout << "SDF_HalfPlane normal:"<< n.toString() << "\n" ;
};

/**
 * Writes a ParaView Phython script explicitly representing the half plane,
 * the plane is given only by its normal, two orthogonal vectors are randomly chosen
 *
 * @param cout      e.g. a file output stream
 * @param c         python object counter for the script (to avoid duplicate names)
 * \return          object counter
 */
int SDF_HalfPlane::writePVPScript(std::ostream & cout, int c) const
{
    std::string name = "obj";
    name.append(std::to_string(c));

    cout << name <<" = Plane()\n" <<
        name << ".Origin = [" << o.x << ", " << o.y << ", " << o.z << "]\n" <<
        name << ".Point1 = [" << p1.x << ", " << p1.y << ", " << p1.z << "]\n" <<
        name << ".Point2 = [" << p2.x << ", " << p2.y << ", " << p2.z << "]\n\n";

    cout << name << "Display = Show(" << name << ",renderView1)\n" <<
        name << "Display.Opacity = 0.1\n";
    //name<< "Display.DiffuseColor = [0., 0., 1.0]\n" << "renderView1.ResetCamera()\n";

    // created 1 object and corresponding Display
    c++;
    return c;
}

} // end namespace CPlantBox
