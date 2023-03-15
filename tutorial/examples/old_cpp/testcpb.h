#include <vector>
#include <string>
#include <iostream>

#include "MappedOrganism.h"


int testcpb()
{

  auto plant = std::make_shared<CPlantBox::MappedPlant>();

  plant->readParameters("/mnt/c/work/CPlantBox/modelparameter/plant/Triticum_aestivum_adapted_2021.xml");
  plant->initialize();
  plant->simulate(30, true);
  std::cout << "Computing Geometry" << std::endl;
  plant->ComputeGeometry();
  auto points = plant->GetGeometry();

  std::cout << "Generated points are "  << points.size() << std::endl;
}