#include <vector>
#include <string>
#include <iostream>

//#include "MappedOrganism.h"
#include "Plant.h"


int testcpb()
{

  auto plant = std::make_shared<CPlantBox::Plant>();
  std::string filename = "fspm2023.xml"; //THIS IS JUST TO TEST WHERE ERRORS COME FROM
  plant->readParameters("../../../modelparameter/structural/plant/fspm2023.xml");

  plant->initialize();
  plant->writeParameters(filename);
  
	return 0;
}