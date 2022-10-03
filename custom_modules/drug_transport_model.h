#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

#include "./submodel_data_structures.h" 

float calculate_diffusion_flux(Cell* pCell, int density_idx, double permeability);

void drug_transport_model_setup();

void drug_transport_model_update( Cell* pCell, Phenotype& phenotype, double dt );

void drug_transport_model_main( double dt );
