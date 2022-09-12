/*
    porin_like_channel.cpp

    Facilitated diffusion of ligand-gated channels. A small submodel. 
*/

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h"
#include "./custom.h"

using namespace BioFVM; 
using namespace PhysiCell;

// #include "./submodel_data_structures.h" 
void add_cell_defaults_porin_like_channel( void ); // Fetch variables here
void setup_tissue_porin_like_channel( void ); // Setup tissue here
void porin_like_channel_model( Cell* pCell, Phenotype& phenotype, double dt );
void porin_like_channel_model_main( double dt ); // Called in main.cpp, assigns transport model to agent