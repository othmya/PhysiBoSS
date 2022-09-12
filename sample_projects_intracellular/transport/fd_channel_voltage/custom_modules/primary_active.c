/*

    Primary active transport: ATP-fueled pumping of substrates

*/


	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; // is already called by standard_update_cell_velocity
	cell_defaults.functions.update_phenotype = facilitated_carrier_model;
	cell_defaults.functions.custom_cell_rule = NULL;

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;


	static int sK_index = microenvironment.find_density_index("sK");
	static int cytB_index = microenvironment.find_density_index("cytB");

	Cell *pCell;

	// Generic variables:
	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0); // for paraview visualization

	// Densities - sK
	cell_defaults.custom_data.add_variable("I_sK", "amol / um^3", parameters.doubles("Initial_I_sK"));
	cell_defaults.custom_data.add_variable("E_sK_near", "amol / um^3", default_microenvironment_options.initial_condition_vector[sK_index]); // Mean nearest substrate conc.
	cell_defaults.custom_data.add_variable("E_sK", "amol/um^3", default_microenvironment_options.initial_condition_vector[sK_index]);		 // Mean external substrate conc.

	// Densities - ATP, ADP
	cell_defaults.custom_data.add_variable("I_ATP", "amol / um^3", parameters.doubles("Initial_I_ATP"));
	cell_defaults.custom_data.add_variable("I_ADP", "amol / um^3", parameters.doubles("Initial_I_ADP")); 

	// Receptors - sK (Fac. diff. carrier) & ATP, ADP
	cell_defaults.custom_data.add_variable("Rcb_sK", "amol / um^3", parameters.doubles("Initial_Rcb_sK"));
	cell_defaults.custom_data.add_variable("Rcb_ATP", "amol / um^3", parameters.doubles("Initial_Rcb_ATP"));
	cell_defaults.custom_data.add_variable("Rcb_total", "amol / um^3", parameters.doubles("Initial_Rcb_sK") + parameters.doubles("Initial_Rcb_ATP")); // Bound receptor density

	// Total receptor - sum of 
	cell_defaults.custom_data.add_variable("Rc_total", "amol/um^3", parameters.doubles("Initial_Rc_total"));


	// Total amounts - sK & ATP, ADP
	cell_defaults.custom_data.add_variable("total_I_sK", "amol", 0.0); // total external sK
	cell_defaults.custom_data.add_variable("total_E_sK", "amol", 0.0); // total external sK
	cell_defaults.custom_data.add_variable("total_I_ADP", "amol", 0.0); // total internal sK
	cell_defaults.custom_data.add_variable("total_I_ATP", "amol", 0.0); // total internal sK
	cell_defaults.custom_data.add_variable("total_sK", "amol", 0.0);   // total sK = I_sK + E_sK

	// Flux - sK
	cell_defaults.custom_data.add_variable("sK_flux", "amol/min", 0.0);

	// Diffusion coefficients, concentration gradients
	cell_defaults.custom_data.add_variable("DC_sK", "um^3/min", microenvironment.diffusion_coefficients[sK_index]); // Diffusion coefficient
	cell_defaults.custom_data.add_variable("D_sK", "mmol/um^3", 0.0);


	// Other constants (Num. of vocels, Initial conditions, etc.)

	cell_defaults.custom_data.add_variable("Initial_E_sK", "amol/um^3", default_microenvironment_options.initial_condition_vector[sK_index]); // Initial external substrate_Y
	cell_defaults.custom_data.add_variable("Initial_I_sK", "amol/um^3", parameters.doubles("Initial_I_sK"));	
	
	cell_defaults.custom_data.add_variable("Initial_I_ATP", "amol/um^3", parameters.doubles("Initial_I_ATP"));
	cell_defaults.custom_data.add_variable("Initial_I_ADP", "amol/um^3", parameters.doubles("Initial_I_ADP"));

	cell_defaults.custom_data.add_variable("total_healthy_cells", "dimensionless", 0.0);				
	cell_defaults.custom_data.add_variable("total_cancer_cells", "dimensionless", 0.0);				

	cell_defaults.custom_data.add_variable("number_of_voxels", "dimensionless", microenvironment.number_of_voxels());



#include "./primary_active.h" 

void primary_active_model(Cell *pCell, Phenotype &phenotype, double dt)
{
	Cell_Definition* pCD = find_cell_definition( pCell->type_name );

	// live.phases[0].entry_function = my_mutation_function;
	// live.phase_link(0,0).arrest_function = volume_arrest_function;
	// // live.display(std::cout);

	static int sK_index = microenvironment.find_density_index("sK"); // microenv substrate index

	static int I_sK_index = pCell->custom_data.find_variable_index("I_sK");
	static int E_sK_near_index = pCell->custom_data.find_variable_index("E_sK_near");
	static int E_sK_index = pCell->custom_data.find_variable_index("E_sK");

	static int I_ATP_index = pCell->custom_data.find_variable_index("I_ATP");
	static int I_ADP_index = pCell->custom_data.find_variable_index("I_ADP");

	double Rc_total_index = pCell->custom_data.find_variable_index("Rc_total");	// [E] Receptor alone (mM)
	double Rcb_sK_index = pCell->custom_data.find_variable_index("Rcb_sK"); // [ES] Bound receptor (mM)
	double Rcb_ATP_index = pCell->custom_data.find_variable_index("Rcb_ATP"); // [ES] Bound receptor (mM)
	double Rcb_total_index = pCell->custom_data.find_variable_index("Rcb_total"); // [ES] Bound receptor (mM)

	static int DC_sK_index = pCell->custom_data.find_variable_index("DC_sK"); // Difussion coefficient
	static int D_sK_index = pCell->custom_data.find_variable_index("D_sK");	  // Concentration gradient

	static int total_E_sK_index = pCell->custom_data.find_variable_index("total_E_sK");
	static int total_I_sK_index = pCell->custom_data.find_variable_index("total_I_sK");
	static int total_I_ATP_index = pCell->custom_data.find_variable_index("total_I_ATP");
	static int total_I_ADP_index = pCell->custom_data.find_variable_index("total_I_ADP");
	static int total_sK_index = pCell->custom_data.find_variable_index("total_sK");

	static int sK_flux_index = pCell->custom_data.find_variable_index("sK_flux");

	static int total_healthy_cells_index = pCell->custom_data.find_variable_index("total_healthy_cells");
	static int total_cancer_cells_index = pCell->custom_data.find_variable_index("total_cancer_cells");

	// Constants
	static double DC_sK = microenvironment.diffusion_coefficients[sK_index];


	double V_cell = pCell->phenotype.volume.total;
	double R_cell = pCell->phenotype.geometry.radius;
	double A_cell = 4 * M_PI * std::pow(R_cell, 2);
	double V_voxel = (microenvironment.mesh.dV);
	double total_num_voxels = microenvironment.number_of_voxels();

	// Kinetic parameters
	static double sK_binding_rate = parameters.doubles("k_1B");
	static double sK_recycling_rate = parameters.doubles("k_n1B");
	static double sk_movement_rate = parameters.doubles("k_2B");
	static double sK_internal_binding_rate = parameters.doubles("k_n2B"); // not employed

	static double ATP_binding_rate = parameters.doubles("k_1A");
	static double ATP_recycling_rate = parameters.doubles("k_n1A");
	static double ATP_hydrolisis_rate = parameters.doubles("k_2A");
	static double ATP_recycling_rate = parameters.doubles("k_n1A");

	static double sK_internal_degradation_rate = parameters.doubles("k_iB"); // Not employed yet

	// sK is externalised at the rate ATP is hydrolised, this is the way of connecting both systems
	sK_recycling_rate = ATP_hydrolisis_rate; // Connect both systems

	// Receptor concentrations
	double Rc_total = pCell->custom_data.variables[Rc_total_index].value; // [E] Receptor alone (mM), total
	double Rcb_total = pCell->custom_data.variables[Rcb_total_index].value; // [E] Receptor alone (mM), total
	double Rcb_sK = pCell->custom_data.variables[Rcb_sK_index].value; // [ES] Bound receptor (mM) to sK
	double Rcb_ATP = pCell->custom_data.variables[Rcb_ATP_index].value; // [ES] Bound receptor (mM)

	// Substrate concentrations
	double E_sK_near = pCell->nearest_density_vector()[sK_index];							   // A density (mM)
	double I_sK = pCell->phenotype.molecular.internalized_total_substrates[sK_index] / V_cell; // Convert to density (mM)
	static double Initial_E_sK = default_microenvironment_options.initial_condition_vector[sK_index];
	static double Initial_I_sK = parameters.doubles("Initial_I_sK");

	double I_ATP = parameters.doubles("Initial_I_ATP");
	double I_ADP = parameters.doubles("Initial_I_ATP");



	// Obtaining the total net internal and external amounts through iterating over the cells

	double total_E_sK = 0.0;
	double total_I_sK = 0.0;
	int total_healthy_cells = 0;
	int total_cancer_cells = 0;

	std::vector<double> I_sK_vector(microenvironment.number_of_voxels());
	std::vector<double> I_ATP_vector(microenvironment.number_of_voxels());
	std::vector<double> I_ADP_vector(microenvironment.number_of_voxels());

	// Counting net internal amount
	#pragma omp parallel for
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell *pC = (*all_cells)[i];
		I_sK_vector[i] = pC->phenotype.molecular.internalized_total_substrates[sK_index];
		total_I_sK += pC->phenotype.molecular.internalized_total_substrates[sK_index];

		I_ATP_vector[i] = pC->phenotype.molecular.internalized_total_substrates[cytB_index];
		I_ADP_vector[i] = pC->phenotype.molecular.internalized_total_substrates[cytB_index];

		// Get total healthy and cancer cells
		if(pCell->type == 1){ total_healthy_cells +=1; }
		else if(pCell->type == 2){ total_cancer_cells += 1; }

	}

	double I_sK_sum = std::accumulate(I_sK_vector.begin(), I_sK_vector.end(), 0.0);
	double I_ATP_sum = std::accumulate(I_ATP_vector.begin(), I_ATP_vector.end(), 0.0);
	double I_ADP_sum = std::accumulate(I_ADP_vector.begin(), I_ADP_vector.end(), 0.0);

	// I_sK_sum and total_I_sK provide the same exact result. They reflect the actual total amount of
	// internal net sK

	std::vector<double> E_sK_vector(microenvironment.number_of_voxels());
	std::vector<double> E_cytB_vector(microenvironment.number_of_voxels());

	#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		E_sK_vector[n] = (microenvironment.density_vector(n)[sK_index]);
		total_E_sK += (microenvironment.density_vector(n)[sK_index]) * V_voxel; // A net amount of sK

		E_cytB_vector[n] = (microenvironment.density_vector(n)[cytB_index]);
		total_E_cytB += (microenvironment.density_vector(n)[cytB_index]) * V_voxel; 
	}

	double E_sK_sum = std::accumulate(E_sK_vector.begin(), E_sK_vector.end(), 0.0);
	double E_sK_mean = E_sK_sum / total_num_voxels; // NOT CORRECT, the accumulate just adds up the number of voxels

	double D_sK = I_sK - E_sK_near; // This provides the direction of the gradient 
	

	// ATP-dependet sK pump 

	double binding_component;
	double recycling_component;
	double movement_component; // could be externalisation or internalisation of a substrate, always against gradient
	double internal_binding_component;
	double internal_degradation_component;

	double ATP_binding_component; 
	double ATP_unbinding_component;
	double ATP_hydrolisis_component;
	double ADP_hrecycling_component;

	// step 1: ATP binding to the receptor
	ATP_binding_component = diffusion_dt * ATP_binding_rate * I_ATP * Rc_total;
	if (ATP_binding_component > Rc_total * diffusion_dt || ATP_binding_component > I_ATP * diffusion_dt )
	{ ATP_binding_component = std::min(Rc_total * V_cell * diffusion_dt, I_ATP * V_cell * diffusion_dt); } // Get the limiting step 
	// double internal_excess_ATP = ATP_binding_component - std::min(Rc_total * V_cell * diffusion_dt, I_ATP * diffusion_dt);
	// { I_ATP += internal_excess_ATP/V_cell; }
	I_ATP -= ATP_binding_component;
	Rcb_ATP += ATP_binding_component;
	Rcb_total += Rcb_ATP;

	// step 2: Some of the ATP unbinds
	ATP_unbinding_component = diffusion_dt * ATP_unbinding_rate * Rcb_ATP;
	if(ATP_unbinding_component > Rcb_ATP * diffusion_dt){ ATP_unbinding_component = Rcb_ATP * V_cell * diffusion_dt; }
	I_ATP += ATP_unbinding_component;
	Rcb_ATP -= ATP_unbinding_component;
	Rcb_total -= Rcb_ATP; 

	// step 3: sK binds the ATP bomb (from inside)
	binding_component = diffusion_dt * sK_binding_rate * I_sK * Rc_total;
	if (binding_component > Rc_total * diffusion_dt || binding_component > I_sK * diffusion_dt )
	{ binding_component = std::min(Rc_total * V_cell * diffusion_dt, I_sK * V_cell * diffusion_dt); } // Get the limiting step 
	// double excess_sK = binding_component - std::min(Rc_total * V_cell * diffusion_dt, I_sK * V_cell * diffusion_dt);
	// { E_sK_near += excess_sK * 1/microenvironment.mesh.dV; }
	I_sK -= binding_component;
	Rcb_sK += binding_component;
	Rcb_total += Rcb_sK;

	// step 4: sK unbinds the ATP bomb (from inside)
	recycling_component = diffusion_dt * sK_recycling_rate * Rcb_sK;
	if (recycling_component > Rcb_sK * diffusion_dt){ recycling_component = Rcb_sK * V_cell * diffusion_dt; } 
	I_sK += recycling_component;
	Rcb_sK -= recycling_component;
	Rcb_total -= Rcb_sK;

	Rc_total -= Rcb_total; // Update total unbound Rc after binding sK and ATP

	// step 5: ATP hydrolisis & sK expulsion
	ATP_hydrolisis_component = diffusion_dt * ATP_hydrolisis_rate * Rcb_ATP;
	if( ATP_hydrolisis_component > Rcb_ATP * diffusion_dt ){ ATP_hydrolisis_component = Rcb_ATP * V_cell * diffusion_dt; }

	movement_component = diffusion_dt * ATP_hydrolisis_rate * Rcb_sK;
	if(movement_component > Rcb_sK * diffusion_dt){ movement_component = Rcb_sK * V_cell * diffusion_dt; }
	
	I_ATP -= ATP_hydrolisis_component;
	I_ADP += ATP_hydrolisis_component;

	I_sK -= movement_component;
	E_sK_near += movement_component;

	Rcb_sK -= ATP_hydrolisis_component;
	Rcb_ATP -= ATP_hydrolisis_component;

	Rcb_total -= Rcb_ATP;
	Rcb_total -= Rcb_sK;
	Rc_total += Rcb_total;

	// step 6: ADP recycling
	ADP_recycling_component = diffusion_dt * ADP_recycling_rate * I_ADP;
	if (ADP_hrecycling_component > I_ADP * diffusion_dt){ ADP_recycling_component = I_ADP * V_cell * diffusion_dt; }
	I_ADP -= ADP_recycling_component;
	I_ATP += ADP_recycling_component;

	
	double v_max = sk_movement_rate * (Rc_total + Rcb_total); 
	double Km = (sK_recycling_rate + sk_movement_rate) / sK_binding_rate; // k_n1 + k_2 / k_1
	double sK_monod = v_max * (I_sK / (Km + I_sK));
	sK_monod *= diffusion_dt; 


	double v_formation = - (v_max * I_sK) / ( Km + I_sK); // Option A
	v_formation = v_formation * V_cell; // Convert mM/min to amol/min

	// double v_formation = - binding_component + recycling_component - endocytosis_component + internal_binding_component;  // Option B
	// double v_formation = - endocytosis_component;																	 // Option C

	
	// Internal degradation of sK - Set to kill HC in 1000 min if Initial_I_sK is 0.3 mM

	double I_sK_degradation;
	if(pCell->type == 1){ I_sK_degradation = parameters.doubles("k_3") * I_sK * diffusion_dt; }
	if(pCell->type == 2){ I_sK_degradation = (5 * parameters.doubles("k_3")) * I_sK * diffusion_dt; } // CC  die twice as fast - more requirement of Glucose
	
	I_sK -= I_sK_degradation + sK_monod;


	// replace hard-coded initial values at each time-step with the new values

	// Substrate concentrations

	// For sK
	pCell->custom_data.variables[E_sK_near_index].value = E_sK_near; // Mean of near_density_vector
	pCell->custom_data.variables[E_sK_index].value = E_sK_mean;		 // Mean of all voxels
	pCell->custom_data.variables[I_sK_index].value = I_sK;			 // Density

	pCell->custom_data.variables[DC_sK_index].value = DC_sK; // Substrate diffusion coefficient
	pCell->custom_data.variables[D_sK_index].value = D_sK;	 // Substrate diffusion coefficient

	// For cytB
	pCell->custom_data.variables[E_cytB_near_index].value = E_cytB_near; // Mean of near_density_vector
	pCell->custom_data.variables[E_cytB_index].value = E_cytB_mean;		 // Mean of all voxels
	pCell->custom_data.variables[I_cytB_index].value = I_cytB;			 // Density

	pCell->custom_data.variables[DC_cytB_index].value = DC_cytB; // Substrate diffusion coefficient
	pCell->custom_data.variables[D_cytB_index].value = D_cytB;	 // Substrate diffusion coefficient

	
	// Fluxes
	pCell->custom_data.variables[sK_flux_index].value = v_formation;
	pCell->custom_data.variables[Rc_sK_index].value = Rc_sK;
	pCell->custom_data.variables[Rcb_sK_index].value = Rcb_sK;
	pCell->phenotype.secretion.net_export_rates[sK_index] = v_formation; // Avoid second multiplication by diffusion_dt

	// Receptor concentrations
	pCell->custom_data.variables[total_Rc_sK_index].value = Rcb_sK + Rc_sK;


	// Net total amounts
	pCell->custom_data.variables[total_E_sK_index].value = total_E_sK;
	pCell->custom_data.variables[total_I_sK_index].value = total_I_sK;
	pCell->custom_data.variables[total_sK_index].value = total_E_sK + total_I_sK;

	pCell->custom_data.variables[total_cancer_cells_index].value = total_cancer_cells;
	pCell->custom_data.variables[total_healthy_cells_index].value = total_healthy_cells;


	// Do not update phenotype if dead
	if (pCell->phenotype.death.dead == true)
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	// Call the standard function -
	// DEACTIVATED FOR NOW. Not the main focus, once model is built, implement it for a more realistic model

	// update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);

	static int necrosis_idx = pCell->phenotype.death.find_death_model_index("Necrosis");
	static int apoptosis_idx = pCell->phenotype.death.find_death_model_index("Apoptosis");

	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 

	int glucose_flux_threshold;



	// Connecting the transport model to the phenotype

	// Monod function that modifies the rate of duplication

	// if (I_sK < parameters.doubles("cytotoxic_threshold") && parameters.doubles("cytotoxic_threshold") != 0.0)
	// {
	// 	std::cout << "Not enough glucose... " << std::endl;
	// 	pCell->start_death( necrosis_idx ); // "Custom" apoptosis 
	// }

	if (I_sK > parameters.doubles("growth_threshold"))
	{
		pCell->phenotype.cycle.data.transition_rate(0,0) = sK_monod * V_cell;
	}

	// if (I_sK < parameters.doubles("growth_threshold"))
	// {
	// 	// pCell->flag_for_division();
	// 	// pCell->phenotype.cycle.advance_cycle(pCell, phenotype, diffusion_dt);
	// 	pCell->phenotype.cycle.data.transition_rate(0,0) = 1e-08;
	// }

	if (std::abs(Rc_sK) < 1.0e-5) 
	{
		glucose_flux_threshold += 1;

		if (glucose_flux_threshold > 3)
		{
			std::cout << "Not enough glucose flux..." << std::endl;
			pCell->start_death( necrosis_idx );
		}

	}


	// this is refered only to the Boolean model inside each agent - Does not affect the model so far, as cells are static
	if (pCell->phenotype.intracellular->need_update())
	{
		if (
			pCell->type == get_cell_definition("default").type && PhysiCell::PhysiCell_globals.current_time >= 0 && pCell->phenotype.intracellular->get_parameter_value("$time_scale") == 0.0)
			pCell->phenotype.intracellular->set_parameter_value("$time_scale", 0.0);

		set_input_nodes(pCell);

		pCell->phenotype.intracellular->update();

		from_nodes_to_cell(pCell, phenotype, dt);
		color_node(pCell);
	}

	return;
}




void primary_active_model_main( double dt )
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