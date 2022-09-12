/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "custom.h" 
#include "../BioFVM/BioFVM.h"  
using namespace BioFVM;

// declare cell definitions here 

std::vector<bool> nodes;

void create_cell_types( void )
{
	std::cout << "creating cell types" << std::endl;
	// set the random seed 
	SeedRandom( parameters.ints("random_seed") );  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 

	// Not including the basic phenotype functionalities, in order to have the simplest model possible

	// Generic MCF breast epithelium cell definition
	// initialize_default_cell_definition(); 
	//			 DO NOT USE, overrides the custom cell rule and the phenotype

	// cell_defaults.parameters.o2_proliferation_saturation = 38.0;  
	// cell_defaults.parameters.o2_reference = 38.0; 

	Cell * pCell;
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; // is already called by standard_update_cell_velocity
	cell_defaults.functions.update_phenotype = porin_like_channel_model; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 

	/*
		Define specific sets of custom variables depending on the transport model
	*/
	
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
	cell_defaults.custom_data.add_variable("sJ_DC", "um² / min", parameters.doubles("sJ_D")); // Diffusion coefficient

    // Fluxes
	cell_defaults.custom_data.add_variable("sJ_flux", "amol/min", 0.0);
	cell_defaults.custom_data.add_variable("sJ_flux_adjusted", "amol/min", 0.0);

    // Concentration gradient
	cell_defaults.custom_data.add_variable("D_sJ", "amol/um³", 0.0);
    
    // Misc.
    cell_defaults.custom_data.add_variable("number_of_voxels", "min",  microenvironment.number_of_voxels());
    cell_defaults.custom_data.add_variable("A_effective", "min", parameters.doubles("n_RC") * parameters.doubles("radius_RC"));
	cell_defaults.custom_data.add_variable("time", "min", PhysiCell_globals.current_time );
	cell_defaults.custom_data.add_variable("total_cancer_cells", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable("total_healthy_cells", "dimensionless", 0.0 );

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 

	static Cell_Definition* pCH = find_cell_definition("healthy");
	static Cell_Definition* pCC = find_cell_definition("cancer_cell");
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/

	build_cell_definitions_maps(); 
	
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	

	
	return; 
}

void setup_tissue(void)
{
	std::cout << "Setting up the tissue..." << std::endl;
	Cell *pCell = NULL;
	static Cell_Definition* pCH = find_cell_definition("healthy");
	static Cell_Definition* pCC = find_cell_definition("cancer_cell");

	// setup_tissue_porin_like_channel(); 

	// PORIN-LIKE CHANNEL: Just one cell in the center of the domain
    double Initial_I_sJ = parameters.doubles("Initial_I_sJ");
    static int sJ_index = microenvironment.find_density_index("sJ");	

    #pragma omp parallel for
	for (int i=0; i < parameters.ints("number_of_cells"); i++) 
	{
		Cell* pC;
		pC = create_cell(get_cell_definition("default")); 

		// Distribute cells randomly on Unit Circle
		// std::vector<double> position = UniformOnUnitCircle();
		// position *= 50.0; position += {0,0,0}; // circle: radius of 50.0, centered at 0,0,0
		// pC -> assign_position( position );
		pC -> assign_position( {0,0,0} );

		// Add initial density of substrate -- Input is a density, convert to a net value
		pC -> phenotype.molecular.internalized_total_substrates[sJ_index] = Initial_I_sJ * pC->phenotype.volume.total;
		
	}

	return;
}

// custom cell phenotype function to run PhysiBoSS when is needed
void porin_like_channel_model(Cell *pCell, Phenotype &phenotype, double dt)
{
	// std::cout << "Running the model..." << std::endl;

	Cell_Definition* pCD = find_cell_definition( pCell->type_name );

	if (phenotype.death.dead == true)
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	// update_cell_and_death_parameters_O2_based(pCell, phenotype, dt); // Not employed

	/*
		
		Run a specific transport model, or all of them in one simulation.
		Useful to include them here, as they could be easily used in a specific cell type.

	*/

	// porin_like_channel_model(pCell, phenotype, dt);
	// ligand_gated_channel_main(pCell, phenotype, dt);
	// voltage_gated_channel_main(pCell, phenotype, dt);
	// primary_active_main(pCell, phenotype, dt);
	// secondary_active_main(pCell, phenotype, dt);

	

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
    static int sJ_DC_index = pCell->custom_data.find_variable_index("sJ_DC");

    // Fluxes
    static int sJ_flux_index = pCell->custom_data.find_variable_index("sJ_flux");
	static int adjusted_sJ_flux_index = pCell->custom_data.find_variable_index("adjusted_sJ_flux");

    // Concentration gradient
    static int D_sJ_index = pCell->custom_data.find_variable_index("D_sJ");

    // Misc.
	static int total_cancer_cells_index = pCell->custom_data.find_variable_index("total_cancer_cells");
	static int total_healthy_cells_index = pCell->custom_data.find_variable_index("total_healthy_cells");
	static int A_effective_index = pCell->custom_data.find_variable_index("A_effective");
    static int time_index = pCell->custom_data.find_variable_index("time");
    pCell->custom_data.variables[time_index].value = PhysiCell_globals.current_time; // Don't forget to round
	
    // Constants
	static double V_cell = pCell->phenotype.volume.total; 
	static double R_cell = pCell->phenotype.geometry.radius; 
	static double A_cell = 4 * M_PI * std::pow(R_cell, 2);
    double V_voxel = ( microenvironment.mesh.dV ) ;
	double total_num_voxels = microenvironment.number_of_voxels();

	double A_channel = parameters.doubles("radius_Rc"); // um², ranges between 20 and 26 

    // Define variables

    // Densities - sJ
    double I_sJ = pCell->phenotype.molecular.internalized_total_substrates[sJ_index] / V_cell; // amol/um³
    double E_sJ_near = pCell->nearest_density_vector()[sJ_index]; // amol/um³
	double Initial_I_sJ = parameters.doubles("Initial_I_sJ");
	double Initial_E_sJ = default_microenvironment_options.initial_condition_vector[sJ_index];

	// Total net amounts - sJ
	double total_E_sJ = 0.0;
	double total_I_sJ = 0.0;
	std::vector<double> I_sJ_vector(microenvironment.number_of_voxels());
	std::vector<double> E_sJ_vector(microenvironment.number_of_voxels());

	// Transporters - net amount, and size of each, their product is the total effective surface
	double n_Rc = parameters.ints("n_Rc"); // net amount of transporters
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
		total_E_sJ += (microenvironment.density_vector(n)[sJ_index]) * V_voxel; 
	}

	double E_sJ_sum = std::accumulate(E_sJ_vector.begin(), E_sJ_vector.end(), 0.0);
	double E_sJ_mean = E_sJ_sum / total_num_voxels; 


	// Define concentration gradient
	double D_sJ = I_sJ - E_sJ_near; // amol/um³ (mM)


	// Simple diffusion of sJ, ONLY when there are open channels (i.e. n_Rc > 0)

	double flux_sJ; 	 // (amol/min)
	double adjusted_flux_sJ; // (amol/min)
	double portion_of_open_channels = NormalRandom(0,1);
	if( portion_of_open_channels < 0.0) { portion_of_open_channels = 0.0; }
	n_Rc *= portion_of_open_channels; // Stochastically open channels at each diffusion_dt
	double A_effective = n_Rc * A_channel; // um²

	flux_sJ = A_effective * sJ_k * (D_sJ);

	// Sanity checks
	if ( flux_sJ < 0.0 ) { // Uptake
		adjusted_flux_sJ = - std::min( std::abs(flux_sJ), std::abs(E_sJ_near * V_voxel) ); }
	else if ( flux_sJ > 0.0 ) { // Secretion
		adjusted_flux_sJ = std::min(flux_sJ, std::abs(I_sJ * V_cell ) ); }

	
	// Replacing hard-coded values with variables at each time-step

	// Densities - sJ
	pCell->custom_data.variables[E_sJ_near_index].value = E_sJ_near; // Mean of near_density_vector
	pCell->custom_data.variables[E_sJ_index].value = E_sJ_mean;		 // Mean of all voxels
	pCell->custom_data.variables[I_sJ_index].value = I_sJ;			 // Density

	// Total net amounts - sJ
	pCell->custom_data.variables[total_E_sJ_index].value = total_E_sJ;
	pCell->custom_data.variables[total_I_sJ_index].value = total_I_sJ;
	pCell->custom_data.variables[total_sJ_index].value = total_E_sJ + total_I_sJ;

	// Transporters - net amount, and size of each, their product is the total effective surface
	pCell->custom_data.variables[A_effective_index].value = A_effective;
	
	// Concentration gradient - sJ
	pCell->custom_data.variables[D_sJ_index].value = D_sJ;

	// Fluxes
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


// Specific coloring function

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector<std::string> output = false_cell_coloring_live_dead(pCell);
	return output;
}


// Functions related to the BoSS part of the simulation, leave for now
void set_input_nodes(Cell* pCell) {}

void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt) {}


void color_node(Cell* pCell){
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
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

