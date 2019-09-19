// Copyright (C) 2016 Daniel Leitner and Andrea Schnepf. See //license.txt for details.

#include <iostream>
#include <fstream>
#include <unistd.h>

#include "../src/RootSystem.h"
#include "../src/analysis.h"

#include "example1.h"
#include "example1_wb_dgf.h"
#include "example2.h"
#include "example3.h"
#include "example4.h"
#include "example5.h"

#include "benchmarks.h"

#include "shehan_SoilCore.h"
#include "shehan_RhizoTubes.h"
#include "shehan_Trenches.h"
// #include "shehan_ScaleElongation.h"
#include "shehan_WC.h"
#include "shehan_ScaleElongation_CL.h"
#include "example_volume.h"
#include "example_exudation.h"

/**
 * Starts an examples (from the examples folder)
 */
int main(int argc, char* argv[])
{
    using namespace CRootBox;

    std::string name="";

    if (argc>1) {
        name= argv[1];
    }

    // example1(); // open parameter file, and output VTP
    // example1_wb_dgf(); // root growth inside a big box to simulate soil surface, open parameter file, and output VTP
    // example2(); // like example 1, but with put geometry
    // example3(); // more than 1 plant
    // example4(); // rhizotubes (an example for a more complex geomety)
    // example5(); // hydrotropism

    // benchmarks();

    //    if (argc>1) {
    //        cout<<"starting simulation: "<< name <<"\n";
    //        shehan_SoilCore(name, false); // put true here to export geometry
    //    } else {
    //        shehan_SoilCore(); // with default values
    //    }
    //
    // shehan_SoilCore("Anagallis_femina_Leitner_2010",true);

    // shehan_RhizoTubes("wheat",true);

    // shehan_Trenches("wheat",true);

    // example_dumux(); // tests the suggested dumux coupling

    example_exudation();
    // shehan_ScaleElongation();
    // shehan_ScaleElongation_CL();
    // shehan_WC();

    // example_volume();
    //
    return 0;

}
