// *** ADDED BY HEADER FIXUP ***
#include <istream>
// *** END ***
// Copyright (C) 2016 Daniel Leitner and Andrea Schnepf. See //license.txt for details.
///This is an under developed version with a simple stem growth.

#include "Plant.h"
#include "analysis.h"

#include <iostream>
#include <fstream>
//#include <unistd.h>
#include "tinyxml2.h"


#include "examples/example1.h"
#include "examples/example2.h"
#include "examples/example3.h"



/**
 * Starts an examples (from the examples folder)
 */
int main(int argc, char* argv[])
{
     string name = "CPlantBox_test_Xiaoran";

    if (argc>1) {
        name= argv[1];
    }

    example1(); // open parameter file, and output VTP
    // example2(); // like example 1, but with put geometry
    // example3(); // more than 1 plant



    /*
     * Open plant and root parameter from a file
     */



//    dxml.Parse(c); // only write if defined
//    dxml.SaveFile("xmltest.xml");


    return(0);

}



