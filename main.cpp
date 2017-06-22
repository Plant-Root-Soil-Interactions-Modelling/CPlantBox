// Copyright (C) 2016 Daniel Leitner and Andrea Schnepf. See //license.txt for details.

#include "Plant.h"
#include "analysis.h"

#include <iostream>
#include <fstream>
#include <unistd.h>

#include "examples/example1.h"
#include "examples/example2.h"
#include "examples/example3.h"



/**
 * Starts an examples (from the examples folder)
 */
int main(int argc, char* argv[])
{
    string name="";

    if (argc>1) {
        name= argv[1];
    }

    example1(); // open parameter file, and output VTP
    // example2(); // like example 1, but with put geometry
    // example3(); // more than 1 plant



    return(0);

}



