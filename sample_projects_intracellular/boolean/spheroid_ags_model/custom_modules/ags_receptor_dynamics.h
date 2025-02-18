#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

void ags_receptor_model_setup();

void ags_receptor_model( Cell* pCell, Phenotype& phenotype, double dt );

void ags_receptor_model_main( double dt );
