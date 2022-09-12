/*
    voltage_gated_channel.cpp

    Facilitated diffusion of voltage-gated channels. A small submodel. 
*/

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 
#include "./custom.h"


using namespace BioFVM; 
using namespace PhysiCell;

// #include "./submodel_data_structures.h" 

void ligand_gated_channel_model_setup();

void ligand_gated_channel_model( Cell* pCell, Phenotype& phenotype, double dt );

void ligand_gated_channel_model_main( double dt );
