
/*

    Facilitated diffusion channel: voltage-gated.

*/

#include "./ligand_gated_channel.h" 

void ligand_gated_channel_model_setup( void )
{
    static int E_prot_A_near_index = pCell->custom_data.find_variable_index("E_prot_A_near");
	static int I_prot_A_index = pCell->custom_data.find_variable_index("I_prot_A");
	static int E_prot_A_index = pCell->custom_data.find_variable_index("E_prot_A");
	
	static int D_prot_A_index = pCell->custom_data.find_variable_index("D_prot_A");
	static int total_E_prot_A_index = pCell->custom_data.find_variable_index("total_E_prot_A");
	static int total_I_prot_A_index = pCell->custom_data.find_variable_index("total_I_prot_A");
	static int total_prot_A_index = pCell->custom_data.find_variable_index("total_prot_A");

	static int prot_A_flux_index = pCell->custom_data.find_variable_index("prot_A_flux");
	static int adjusted_prot_A_flux_index = pCell->custom_data.find_variable_index("adjusted_prot_A_flux");


	// Constants

	static double V_cell = pCell->phenotype.volume.total; 
	static double R_cell = pCell->phenotype.geometry.radius; // Around 8, but I don't get the difference between both
	static double R_cell_2 = std::pow(V_cell * 3.0/4.0, 1.0/3.0); // 12.3, use this one for now, check the radius in PCUG
	static double A_cell = 4 * M_PI * std::pow(R_cell, 2);
	static double A_channel = 2e-07; // um², ranges between 20 and 26 

	double V_voxel = ( microenvironment.mesh.dV ) ;
	double total_num_voxels = microenvironment.number_of_voxels();


	static double prot_A_k = parameters.doubles("k_prot_A");
	static double total_net_receptors = parameters.ints("number_of_receptors");

    static int prot_A_index = microenvironment.find_density_index("prot_A");

	double E_prot_A_near = pCell->nearest_density_vector()[prot_A_index]; // A density
	double I_prot_A = pCell->phenotype.molecular.internalized_total_substrates[prot_A_index] / V_cell; // Convert to density
	// I_prot_A just refers to one cell, not useful when more than one cell in the simulation

	double total_E_prot_A = 0.0;
	double total_I_prot_A = 0.0;
	
	static double Initial_E_prot_A = default_microenvironment_options.initial_condition_vector[prot_A_index];
	static double Initial_I_prot_A = parameters.doubles("Initial_I_prot_A");



	// Get total External net amount of substrate
	std::vector<double> E_prot_A_vector(microenvironment.number_of_voxels());


	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{ 	
		E_prot_A_vector[n] = (microenvironment.density_vector(n)[prot_A_index]);
		total_E_prot_A += (microenvironment.density_vector(n)[prot_A_index]) * V_voxel; // A net amount of prot_A

	}

	double E_prot_A_sum = std::accumulate(E_prot_A_vector.begin(), E_prot_A_vector.end(), 0.0);
	double E_prot_A_mean = E_prot_A_sum / total_num_voxels;


	// Counting net internal amount of substrate
	std::vector<double> I_prot_A_vector((*all_cells).size());
	
	#pragma omp parallel for
	for( int i=0; i < (*all_cells).size(); i++ )
	{
		//std::cout << "This is cell number " << i << std::endl; DOES NOT ITERATE CORRECTLY 
		Cell* pC = (*all_cells)[i]; 
		I_prot_A_vector[i] = pC -> phenotype.molecular.internalized_total_substrates[prot_A_index] / V_cell;
		total_I_prot_A += pC->phenotype.molecular.internalized_total_substrates[prot_A_index];
	}

	double I_prot_A_sum = std::accumulate(I_prot_A_vector.begin(), I_prot_A_vector.end(), 0);


	double D_prot_A = I_prot_A - E_prot_A_near;


	// replace hard-coded initial values at each time-step with the new values

	pCell->custom_data.variables[E_prot_A_near_index].value = E_prot_A_near; // Mean of near_density_vector 
	pCell->custom_data.variables[E_prot_A_index].value = E_prot_A_mean; // Mean of all voxels
	pCell->custom_data.variables[I_prot_A_index].value = I_prot_A; // Density
	pCell->custom_data.variables[D_prot_A_index].value = D_prot_A; // Gradient
	pCell->custom_data.variables[total_E_prot_A_index].value = total_E_prot_A; 
	pCell->custom_data.variables[total_I_prot_A_index].value = total_I_prot_A; 
	pCell->custom_data.variables[total_prot_A_index].value = total_I_prot_A + total_E_prot_A; // Net amount of prot_A
    
    return;

}





void voltage_gated_channel_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	
	// The model is based on Fick's law, but here the effective surface through which molecules pass depends 
	// on the amount of channels within an agents' surface.

	double flux; 	 // (amol/min)
	double adjusted_flux; // (amol/min)

	if( D_prot_A != 0)
	{
		flux = A_channel * prot_A_k * (D_prot_A);

	}
	else if ( D_prot_A == 0)
	{
		flux = 0;
	}

	// Sanity checks
	if ( flux_prot_A < 0.0 ) { // Uptake
		adjusted_flux_prot_A = - std::min( std::abs(flux_prot_A), std::abs(E_prot_A_near * V_voxel) );
	}
	else if ( flux_prot_A > 0.0 ) { // Secretion
		adjusted_flux_prot_A = std::min(flux_prot_A, std::abs(I_prot_A * V_cell ) );
	}
		
	pCell->custom_data.variables[prot_A_flux_index].value = flux_prot_A * 0.01;
	pCell->custom_data.variables[adjusted_prot_A_flux_index].value = adjusted_flux_prot_A * 0.01;

	pCell->phenotype.secretion.net_export_rates[prot_A_index] = adjusted_flux_prot_A; 


	// Do not update phenotype if dead 
    if( pCell->phenotype.death.dead == true )
    { 
		pCell->functions.update_phenotype = NULL; 		
		return; 
	}
   
    // Call the standard function - 
	// DEACTIVATED FOR NOW. Not the main focus, once model is built, implement it for a more realistic model

    // update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);


	// this is refered only to the Boolean model inside each agent - Does not affect the model so far, as cells are static
	if (pCell->phenotype.intracellular->need_update())
	{	
		if (
			pCell->type == get_cell_definition("default").type
			&& PhysiCell::PhysiCell_globals.current_time >= 0 
			&& pCell->phenotype.intracellular->get_parameter_value("$time_scale") == 0.0
		)
			pCell->phenotype.intracellular->set_parameter_value("$time_scale", 0.0);

		set_input_nodes(pCell);

		pCell->phenotype.intracellular -> update();
		
		from_nodes_to_cell(pCell, phenotype, dt);
		color_node(pCell);
	}

	
	return;

}




void voltage_gated_channel_model_main( double dt )
{
	#pragma omp parallel for 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		Cell* pC = (*all_cells)[i]; 
		if( pC->phenotype.death.dead == false )
		{ ligand_gated_channel_model( pC, pC->phenotype , dt ); }
	}
	
	return;
}