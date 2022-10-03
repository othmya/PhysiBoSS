/*
	
	ags_receptor_dynamics.cpp

	each drug that enters does so with its own mechanism
*/

#include "./drug_transport_model.h" 
#include "./boolean_model_interface.h" 

using namespace PhysiCell; 

Submodel_Information ags_receptor_info;


float calculate_diffusion_flux(Cell* pCell, int density_idx, double permeability){
	
	// Cell and voxel geometries
	static double voxel_volume = microenvironment.mesh.dV;
	
	double cell_volume = pCell->phenotype.volume.total;
	double cell_radius = pCell->phenotype.geometry.radius;
	double cell_surface = 4 * M_PI * std::pow(cell_radius, 2);

	double density_ext = pCell->nearest_density_vector()[density_idx];	 // A density (mM)
	double density_int = pCell->phenotype.molecular.internalized_total_substrates[density_idx];
	density_int /= cell_volume; // divide int tot substrate to get concentration

	float flux = permeability * (density_int - density_ext) * cell_surface;

	/*
	TODO: Fix this -- write correctly 


	if ( flux < 0.0 ) { // Uptake
		flux = - std::min( std::abs(flux), std::abs(density_ext * voxel_volume) );
	}
	else if ( flux > 0.0 ) { // Secretion
		flux = std::min(flux, std::abs(density_int * cell_volume ) );
	}
	*/

	return flux;
}

void drug_transport_model_setup()
{
    ags_receptor_info.name = "AGS model"; 
	ags_receptor_info.version = "0.0.1";
	
    ags_receptor_info.main_function = drug_transport_model_update; 

	// what custom data do I need?
	ags_receptor_info.cell_variables.push_back( "activation_threshold" );

	// drug_X variables
	ags_receptor_info.cell_variables.push_back( "drug_X_permeability" );

	// drug_Y variables
	ags_receptor_info.cell_variables.push_back( "drug_Y_permeability" );

	ags_receptor_info.register_model();

	return;
}

void drug_transport_model_update( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.death.dead == true )
	{ return; }

	// Fetch microenvironment variables
	static int drug_X_idx = microenvironment.find_density_index("drug_X");
	static int drug_Y_idx = microenvironment.find_density_index("drug_Y");

	double drug_X_permeability = parameters.doubles("drug_X_permeability");
	double drug_Y_permeability = parameters.doubles("drug_Y_permeability");

	double drug_X_flux = calculate_diffusion_flux(pCell, drug_X_idx, drug_X_permeability);
	double drug_Y_flux = calculate_diffusion_flux(pCell, drug_Y_idx, drug_Y_permeability);
	
	// Mapping to net export rate of agent
	pCell->phenotype.secretion.net_export_rates[drug_X_idx] = drug_X_flux;
	pCell->phenotype.secretion.net_export_rates[drug_Y_idx] = drug_Y_flux;

	return;
}

void drug_transport_model_main( double dt )
{
	#pragma omp parallel for 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		Cell* pCell = (*all_cells)[i]; 
		if( pCell->phenotype.death.dead == true )
		{ continue; }
		drug_transport_model_update( pCell, pCell->phenotype , dt );
	}
	return;
}
