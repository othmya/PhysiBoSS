/*
	
	ags_receptor_dynamics.cpp

	each drug that enters does so with its own mechanism
*/

#include "./drug_transport_model.h" 
#include "./boolean_model_interface.h" 

using namespace PhysiCell; 

Submodel_Information ags_receptor_info;


float calculate_diffusion_flux(Cell* pCell, int density_idx, double permeability, std::string drug_name){
	
	// Cell and voxel geometries
	static double voxel_volume = microenvironment.mesh.dV;

	double cell_volume = pCell->phenotype.volume.total;
	double cell_radius = pCell->phenotype.geometry.radius;
	double cell_surface = 4 * M_PI * std::pow(cell_radius, 2);

	double density_ext = pCell->nearest_density_vector()[density_idx];	 // A density (mM)

	double density_int = pCell->phenotype.molecular.internalized_total_substrates[density_idx];
	density_int /= cell_volume; // divide int tot substrate to get concentration

	// float flux = permeability * (density_int - density_ext) * cell_surface;

	// if (flux < 0.0 && std::abs(flux) > std::abs(density_ext * voxel_volume))
	// {
	// 	flux = - std::abs(density_ext * voxel_volume);
	// 	flux /= 2;
	// }
	// else if (flux > 0.0 && flux > density_int * cell_volume)
	// {
	// 	flux = std::abs(density_int * cell_volume);
	// 	flux /= 2;
	// }

	// Then map to the custom data for post-simulation analysis
	// std::string substrate_external_idx = drug_name + "_external_density";
	// std::string substrate_internal_idx = drug_name + "_internal_density";

	// static int external_density_custom_data_idx = pCell->custom_data.find_variable_index(substrate_external_idx);
	// static int internal_density_custom_data_idx = pCell->custom_data.find_variable_index(substrate_internal_idx);

	// pCell->custom_data.variables[external_density_custom_data_idx].value = density_ext;
	// pCell->custom_data.variables[internal_density_custom_data_idx].value = density_int;

	// std::cout << substrate_external_idx << " : " << density_ext << std::endl; 
	// std::cout << substrate_internal_idx << " : " << density_int << std::endl; 


	// Previous adjusting 
	// if ( flux < 0.0 ) { // Uptake
	// 	flux = - std::min( std::abs(flux), std::abs(density_ext * voxel_volume) );
	// }
	// else if ( flux > 0.0 ) { // Secretion
	// 	flux = std::min(flux, std::abs(density_int * cell_volume ) );
	// }


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
	ags_receptor_info.cell_variables.push_back( "drug_X_external_density" );
	ags_receptor_info.cell_variables.push_back( "drug_X_internal_density" );

	// drug_Y variables
	ags_receptor_info.cell_variables.push_back( "drug_Y_permeability" );
	ags_receptor_info.cell_variables.push_back( "drug_Y_external_density" );
	ags_receptor_info.cell_variables.push_back( "drug_Y_internal_density" );

	// oxygen
	ags_receptor_info.cell_variables.push_back( "oxygen_external_density" );
	ags_receptor_info.cell_variables.push_back( "oxygen_internal_density" );

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

	double drug_X_flux = calculate_diffusion_flux(pCell, drug_X_idx, drug_X_permeability, "drug_X");
	pCell->phenotype.secretion.net_export_rates[drug_X_idx] = drug_X_flux;
	
	double drug_Y_flux = calculate_diffusion_flux(pCell, drug_Y_idx, drug_Y_permeability, "drug_Y");
	pCell->phenotype.secretion.net_export_rates[drug_Y_idx] = drug_Y_flux;


	// std::cout << "drug_X_flux " << drug_X_flux << std::endl;
	// std::cout << "drug_Y_flux " << drug_Y_flux << std::endl;
	
	// map oxygen to custom data (test)
	// static double V_cell = pCell->phenotype.volume.total; 
	// static int oxy_idx = microenvironment.find_density_index("oxygen");
	// float oxy_int = pCell->phenotype.molecular.internalized_total_substrates[oxy_idx] / V_cell;
	// float oxy_ext = pCell->nearest_density_vector()[oxy_idx];

	// static int oxy_int_index = pCell->custom_data.find_variable_index("oxygen_internal_density");
	// static int oxy_ext_index = pCell->custom_data.find_variable_index("oxygen_external_density");
	// pCell->custom_data.variables[oxy_int_index].value = oxy_int;
	// pCell->custom_data.variables[oxy_ext_index].value = oxy_ext;


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

		// Adding values to custom data
		static int time_index = pCell->custom_data.find_variable_index("time");
		pCell->custom_data.variables[time_index].value = PhysiCell_globals.current_time;

	}

	return;
}
