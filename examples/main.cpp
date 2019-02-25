// Copyright (C) 2016 Daniel Leitner and Andrea Schnepf. See //license.txt for details.
///This is an under developed version with a simple stem growth.
#include <istream>
#include "Plant.h"
#include "Organ.h"
#include "analysis.h"

#include <iostream>
#include <fstream>
//#include <unistd.h>
#include "tinyxml2.h"


#include "example1.h"
#include "example2.h"
//#include "example3.h"
//#include "example_scenegraph.h"





/**
 * Starts an examples (from the examples folder)
 */
int main(int argc, char* argv[])
{
	string name = "PMA2018";

	if (argc>1) {
		name= argv[1];
	}

	CPlantBox::example1(); // open parameter file, and output VTP
	// example2(); // like example 1, but with put geometry
	// example3(); // more than 1 plant

	// example_scenegraph();

	return(0);

}





