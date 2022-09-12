/*
    Aquaporin-like channels

*/

#include "./porin_like_channel.h"
#include "./custom.h"

void add_cell_defaults_porin_like_channel( void ) // SUBSTRATE IS sJ
{
    Cell *pCell;

    // Microenvironment indices
	static int sJ_index = microenvironment.find_density_index("sJ");
	static int H2O_index = microenvironment.find_density_index("H2O");

    // Densities - sJ
	cell_defaults.custom_data.add_variable("I_sJ", "amol / um^3", parameters.doubles("Initial_I_sJ"));
	cell_defaults.custom_data.add_variable("E_sJ_near", "amol / um^3", default_microenvironment_options.initial_condition_vector[sJ_index]); // Mean nearest substrate conc.
	cell_defaults.custom_data.add_variable("E_sJ", "amol/um^3", default_microenvironment_options.initial_condition_vector[sJ_index]);   	 // Mean external substrate conc.

    // Total net amounts
	cell_defaults.custom_data.add_variable("total_I_sJ", "amol", 0.0);
	cell_defaults.custom_data.add_variable("total_E_sJ", "amol", 0.0);
	cell_defaults.custom_data.add_variable("total_net_sJ", "amol", 0.0);
    
    // Transporters - net amount, and size of each, their product is the total effective surface
	cell_defaults.custom_data.add_variable("n_Rc", "dimensionless", parameters.ints("n_Rc"));
	cell_defaults.custom_data.add_variable("radius_Rc", "um", parameters.doubles("radius_Rc"));

    // Kinetic parameters
	cell_defaults.custom_data.add_variable("sJ_k", "um / min", parameters.doubles("sJ_k")); // Permeability of sJ
	cell_defaults.custom_data.add_variable("sJ_D", "um² / min", parameters.doubles("sJ_D")); // Diffusion coefficient

    // Fluxes
	cell_defaults.custom_data.add_variable("sJ_flux", "amol/min", 0.0);
	cell_defaults.custom_data.add_variable("sJ_flux_adjusted", "amol/min", 0.0);

    // Concentration gradient
	cell_defaults.custom_data.add_variable("sJ_DC", "amol/um³", 0.0);
    
    // Misc.
    cell_defaults.custom_data.add_variable("number_of_voxels", "min",  microenvironment.number_of_voxels());
	cell_defaults.custom_data.add_variable("time", "min", PhysiCell_globals.current_time );


    return;
}

void setup_tissue_porin_like_channel(void)
{

    // So far, place just one cell 
    double Initial_I_sJ = parameters.doubles("Initial_I_sJ");
    static int sJ_index = microenvironment.find_density_index("sJ");	

    #pragma omp parallel for
	for (int i=0; i < parameters.ints("number_of_cells"); i++) 
	{
		Cell* pC;
		pC = create_cell(get_cell_definition("default")); 

		// Distribute cells randomly on Unit Circle
		std::vector<double> position = UniformOnUnitCircle();
		position *= 50.0; position += {0,0,0}; // circle: radius of 50.0, centered at 0,0,0
		// pC -> assign_position( position );
		pC -> assign_position( {0,0,0} );

		// Add initial density of substrate -- Input is a density, convert to a net value
		pC -> phenotype.molecular.internalized_total_substrates[sJ] = Initial_I_sJ * pC->phenotype.volume.total;
		
	}

    return;
}



/*
    The model is based on Fick's law, but here the effective surface through which molecules pass depends 
    on the amount of channels within an agents' surface.

    Step-function ODE, where the condition is given by the state of the channels.

*/

void porin_like_channel_model(Cell *pCell, Phenotype &phenotype, double dt)
{
    Cell_Definition* pCD = find_cell_definition( pCell->type_name );

    // Microenvironment indices
    static int sJ_index = microenvironment.find_density_index("sJ");
	static int H2O_index = microenvironment.find_density_index("H2O");


    // Fetch custom data:

    // Densities - sJ
    static int E_sJ_near_index = pCell->custom_data.find_variable_index("E_sJ_near");
	static int I_sJ_index = pCell->custom_data.find_variable_index("I_sJ");
	static int E_sJ_index = pCell->custom_data.find_variable_index("E_sJ");

    // Total net amounts - sJ
	static int total_E_sJ_index = pCell->custom_data.find_variable_index("total_E_sJ");
	static int total_I_sJ_index = pCell->custom_data.find_variable_index("total_I_sJ");
	static int total_sJ_index = pCell->custom_data.find_variable_index("total_net_sJ");

    // Transporters - net amount, and size of each, their product is the total effective surface
    static int n_Rc_index = pCell->custom_data.find_variable_index("n_Rc");
    static int radius_Rc_index = pCell->custom_data.find_variable_index("radius_Rc");

    // Kinetic parameters
    static int sJ_k_index = pCell->custom_data.find_variable_index("sJ_k");
    static int sJ_D_index = pCell->custom_data.find_variable_index("sJ_D");

    // Fluxes
    static int sJ_flux_index = pCell->custom_data.find_variable_index("sJ_flux");
	static int adjusted_sJ_flux_index = pCell->custom_data.find_variable_index("adjusted_sJ_flux");

    // Concentration gradient
    static int DC_sJ_index = pCell->custom_data.find_variable_index("DC_sJ");

    // Misc.
    static int time_index = pCell->custom_data.find_variable_index("time");
    pCell->custom_data.variables[time_index].value = PhysiCell_globals.current_time; // Don't forget to round

    // Constants
	static double V_cell = pCell->phenotype.volume.total; 
	static double R_cell = pCell->phenotype.geometry.radius; 
	static double A_cell = 4 * M_PI * std::pow(R_cell, 2);
    double V_voxel = ( microenvironment.mesh.dV ) ;
	double total_num_voxels = microenvironment.number_of_voxels();

	static double A_channel = parameters.doubles("radius_Rc"); // um², ranges between 20 and 26 

    // Define variables

    // Densities - sJ
    double I_sJ = pCell->phenotype.molecular.internalized_total_substrates[prot_A_index] / V_cell; // amol/um³
    double E_sJ_near = pCell->nearest_density_vector()[sJ_index]; // amol/um³
	double Initial_I_sJ = parameters.doubles("Initial_I_sJ");
	double Initial_E_sJ = default_microenvironment_options.initial_condition_vector[sJ_index];

	// Total net amounts - sJ
	double total_E_sJ = 0.0;
	double total_I_sJ = 0.0;
	std::vector<double> I_sJ_vector(microenvironment.number_of_voxels());
	std::vector<double> E_sJ_vector(microenvironment.number_of_voxels());

	// Transporters - net amount, and size of each, their product is the total effective surface
	double n_Rc = parameters.doubles("n_Rc"); // net amount of transporters
	double radius_Rc = parameters.doubles("radius_Rc"); // um

	// Kinetic parameters
    double sJ_k = parameters.doubles("sJ_k");
    double sJ_DC = parameters.doubles("sJ_DC");

	// Misc.
	int total_healthy_cells = 0.0;
	int total_cancer_cells = 0.0;


	// Fetch total Internal and External sJ

	#pragma omp parallel for
	for (int i = 0; i < (*all_cells).size(); i++) // Internal
	{
		Cell *pC = (*all_cells)[i];
		I_sJ_vector[i] = pC->phenotype.molecular.internalized_total_substrates[sJ_index];
		total_I_sJ += pC->phenotype.molecular.internalized_total_substrates[sJ_index];

		// Get total healthy and cancer cells
		if(pC->type == 1){ total_healthy_cells +=1; }
		if(pC->type == 2){ total_cancer_cells += 1; }
	}

	pCell->custom_data.variables[total_cancer_cells_index].value = total_cancer_cells;
	pCell->custom_data.variables[total_healthy_cells_index].value = total_healthy_cells;

	double I_sJ_sum = std::accumulate(I_sJ_vector.begin(), I_sJ_vector.end(), 0.0);

	#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++) // External
	{
		E_sJ_vector[n] = (microenvironment.density_vector(n)[sJ_index]);
		total_E_sJ += (microenvironment.density_vector(n)[sY_index]) * V_voxel; // A net amount of sY
	}

	double E_sJ_sum = std::accumulate(E_sJ_vector.begin(), E_sJ_vector.end(), 0.0);
	double E_sY_mean = E_sY_sum / total_num_voxels; 


	// Define concentration gradient
	double sJ_D = I_sJ - E_sJ_near; // amol/um³ (mM)


	// Simple diffusion of sJ, ONLY when there are open channels (i.e. n_Rc > 0)

	double flux; 	 // (amol/min)
	double adjusted_flux; // (amol/min)
	double portion_of_open_channels = 1.0; //= NormalRandom(0,1);
	if( portion_of_open_channels < 0.0) { portion_of_open_channels = 0.0; }
	n_Rc *= portion_of_open_channels; // Stochastically open channels at each diffusion_dt
	double A_effective = n_Rc * A_channel; // um²

	flux = A_channel * sJ_k * (D_sJ);

	// Sanity checks
	if ( flux_sJ < 0.0 ) { // Uptake
		adjusted_flux_sJ = - std::min( std::abs(flux_sJ), std::abs(E_sJ_near * V_voxel) ); }
	else if ( flux_sJ > 0.0 ) { // Secretion
		adjusted_flux_sJ = std::min(flux_sJ, std::abs(I_sJ * V_cell ) ); }
		
	pCell->custom_data.variables[sJ_flux_index].value = flux_sJ * diffusion_dt;
	pCell->custom_data.variables[adjusted_sJ_flux_index].value = adjusted_flux_sJ * diffusion_dt;
	pCell->phenotype.secretion.net_export_rates[sJ_index] = adjusted_flux_sJ; 


	// Do not update phenotype if dead 
    if( pCell->phenotype.death.dead == true )
    { 
		pCell->functions.update_phenotype = NULL; 		
		return; 
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




void porin_like_channel_model_main( double dt )
{
	#pragma omp parallel for 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		Cell* pC = (*all_cells)[i]; 
		if( pC->phenotype.death.dead == false )
		{ porin_like_channel_model( pC, pC->phenotype , dt ); }
	}
	
	return;
}

