#include "cpb.h"

int main ( int argc, char *argv[] ) {
    std::string para = argv[1];
    std::string time = argv[2];
    double time2 = atof(time.c_str());
	CPlantBox::cpb(para, time2);



}