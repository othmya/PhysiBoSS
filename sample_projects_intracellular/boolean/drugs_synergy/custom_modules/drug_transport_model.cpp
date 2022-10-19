/*
	
	ags_receptor_dynamics.cpp

	each drug that enters does so with its own mechanism
*/

#include "./drug_transport_model.h" 
#include "./boolean_model_interface.h" 

using namespace PhysiCell; 

Submodel_Information ags_receptor_info;


double calculate_diffusion_flux(Cell* pCell, int density_idx, double permeability, std::string drug_name){
	
	// Cell and voxel geometries
	static double voxel_volume = microenvironment.mesh.dV;

	double cell_volume = pCell->phenotype.volume.total;
	double cell_radius = pCell->phenotype.geometry.radius;
	double cell_surface = 4 * M_PI * std::pow(cell_radius, 2);

	double density_ext = pCell->nearest_density_vector()[density_idx];	 // A density (mM)

	double density_int = pCell->phenotype.molecular.internalized_total_substrates[density_idx];
	density_int /= cell_volume; // divide int tot substrate to get concentration

	double flux = 0.0;
	flux = permeability * (density_int - density_ext) * cell_surface;

	// Adjusting for extreme cases -- For now only applies to uptake


	// -- Option A 
	// First, check for the equilibirum concentration

	// double flux_net_rate = flux*diffusion_dt; // in amol
	// std::string drug_initial_external_conc = drug_name + "_pulse_concentration";

	// float initial_external_density = parameters.doubles(drug_initial_external_conc); // in mM
	// float initial_external_net_amount = initial_external_density * voxel_volume; // in mmol

	// float initial_internal_density = 0.0; // mM
	// float initial_internal_net_amount = initial_internal_density * cell_volume; // amol

	// float total_net_amount = initial_external_net_amount + initial_internal_net_amount;
	
	// // At equilibirum
	// float eq_external_density = (total_net_amount / 2) * voxel_volume;
	// float eq_internal_density = (total_net_amount / 2) * cell_volume;

	// // Check if net export rate step will pass equilibirum
	// float internal_density_at_t1 = density_int + (flux_net_rate / cell_volume);

	// if ( flux < 0.0 && internal_density_at_t1 > eq_internal_density ){
	// 	// Introduce the difference between the internal density and 
	// 	flux = density_int - eq_internal_density; // remember: negative is uptake
	// 	flux *= cell_volume;
	// }
	

	// -- Option B: Using only the delta of net amount, not concentration

	if( flux < 0.0 && std::abs(flux) >= density_ext * voxel_volume ){
		flux = (density_int * cell_volume) - (density_ext * voxel_volume); // amol
	} else if ( flux > 0.0 && flux >= density_int * cell_volume){
		flux = (density_int * cell_volume) - (density_ext * voxel_volume); // amol
	}

	// Then map to the custom data for post-simulation analysis
	std::string substrate_external_idx = drug_name + "_external_density";
	std::string substrate_internal_idx = drug_name + "_internal_density";

	static int external_density_custom_data_idx = pCell->custom_data.find_variable_index(substrate_external_idx);
	static int internal_density_custom_data_idx = pCell->custom_data.find_variable_index(substrate_internal_idx);

	pCell->custom_data.variables[external_density_custom_data_idx].value = density_ext;
	pCell->custom_data.variables[internal_density_custom_data_idx].value = density_int;

	// std::cout << substrate_external_idx << " : " << density_ext << std::endl; 
	// std::cout << substrate_internal_idx << " : " << density_int << std::endl; 

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
	double drug_Y_flux = calculate_diffusion_flux(pCell, drug_Y_idx, drug_Y_permeability, "drug_Y");

	pCell->phenotype.secretion.net_export_rates[drug_X_idx] = drug_X_flux;
	pCell->phenotype.secretion.net_export_rates[drug_Y_idx] = drug_Y_flux;
	
	// map oxygen to custom data -- TODO: Encapsulate this in a function
	static double V_cell = pCell->phenotype.volume.total; 
	static int oxy_idx = microenvironment.find_density_index("oxygen");
	float oxy_int = pCell->phenotype.molecular.internalized_total_substrates[oxy_idx] / V_cell;
	float oxy_ext = pCell->nearest_density_vector()[oxy_idx];

	static int oxy_int_index = pCell->custom_data.find_variable_index("oxygen_internal_density");
	static int oxy_ext_index = pCell->custom_data.find_variable_index("oxygen_external_density");
	pCell->custom_data.variables[oxy_int_index].value = oxy_int;
	pCell->custom_data.variables[oxy_ext_index].value = oxy_ext;

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

		int cell_microenv_index = pCell->get_current_voxel_index();
		bool DC_node = microenvironment.is_dirichlet_node(cell_microenv_index);
		bool agent_out_domain = pCell->is_out_of_domain;

		if (cell_microenv_index >= 0 && DC_node == false && agent_out_domain == false) // Avoid segfault of NER
			drug_transport_model_update( pCell, pCell->phenotype , dt );

		// Adding time to custom data
		static int time_index = pCell->custom_data.find_variable_index("time");
		pCell->custom_data.variables[time_index].value = PhysiCell_globals.current_time;

	}

	return;
}