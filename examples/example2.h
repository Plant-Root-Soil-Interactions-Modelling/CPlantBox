/**



 * Example 2



 *



 * Same as Example 1, but with a plant container confining root growth



 *



 * Example container are:



 * 1. Cylindrical soil core,



 * 2. Flat square rhizotron,



 * 3. Rotated square rhizotron



 *



 * Additionally, exports the confining geometry as paraview pyhton script



 * (run file in Paraview by Tools->Python Shell, Run Script)



 */

namespace CPlantBox {

using namespace std;







void example2()



{



    Plant rootsystem;







    string name = "Zea_mays_5_Leitner_2014";







    /*



     * Plant and root parameter from a file



     */



    rootsystem.openFile(name);







    /*



     * Set geometry



     */







    // 1.creates a cylindrical soil core with top radius 5 cm, bot radius 5 cm, height 50cm, not square but circular



    SDF_PlantContainer soilcore(5,5,40,false);







    //2. creates a square 27*27 cm containter with height 1.5 cm (used in parametrisation experiment)



    SDF_PlantBox rhizotron(1.4,27,27);







    //3. creates a square rhizotron r*r, with height h, rotated around the x-axis for angle alpha



    double r = 20;



    double h = 4;



    double alpha = 45.;



    SDF_PlantContainer rhizotron2(r,r,h,true);



    Vector3d posA = Vector3d(0,r,-h/2); // origin before rotation



    Matrix3d A = Matrix3d::rotX(alpha/180.*3.14);



    posA = A.times(posA); // origin after rotation



    SDF_RotateTranslate rotatedRhizotron(&rhizotron2,alpha,0,posA.times(-1));







    //rootsystem.setGeometry(&soilcore); // pick one of the geometries



	rootsystem.setGeometry(&rhizotron); // pick one of the geometries



    //rootsystem.setGeometry(&rotatedRhizotron); // pick one of the geometries







    /*



     * Initialize



     */



    rootsystem.initialize(); //it is important to call initialize() after setGeometry()







    /*



     * Simulate



     */



    double simtime = 60;



    double dt = 1;



    int N = round(simtime/dt);



    for (int i=0; i<N; i++) {



        rootsystem.simulate(dt);



    }







    /*



     * Export results (as vtp)



     */



    rootsystem.write(name+std::string(".vtp"));







    /*



     * Export container geometry as Paraview Python script (run file in Paraview by Tools->Python Shell, Run Script)



     */



    rootsystem.write(name+ std::string(".py"));







    cout << "Finished with a total of " << rootsystem.getNumberOfNodes()<< " nodes\n";







}
}


