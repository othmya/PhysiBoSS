/*
    header file for primary_active.c 
*/

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 
#include "./custom.h"


using namespace BioFVM; 
using namespace PhysiCell;

// #include "./submodel_data_structures.h" 

void primary_active_model_setup();

void primary_active_model( Cell* pCell, Phenotype& phenotype, double dt );

void ligand_gated_channel_model_main( double dt );
