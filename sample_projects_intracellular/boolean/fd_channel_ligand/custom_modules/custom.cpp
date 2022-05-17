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
	cell_defaults.functions.update_phenotype = fd_channel_ligand_model; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 

		/*
		Define specific sets of custom variables depending on the transport model
	*/

	static int sJ_index = microenvironment.find_density_index("sJ");

	// Generic variables:
	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0); // for paraview visualization

	// Densities - sJ
	cell_defaults.custom_data.add_variable("I_sJ", "amol / um^3", parameters.doubles("Initial_I_sJ"));
	cell_defaults.custom_data.add_variable("E_sJ_near", "amol / um^3", default_microenvironment_options.initial_condition_vector[sJ_index]); // Mean nearest substrate conc.
	cell_defaults.custom_data.add_variable("E_sJ", "amol/um^3", default_microenvironment_options.initial_condition_vector[sJ_index]);		 // Mean external substrate conc.

	// Densities - Lg, Rc
	cell_defaults.custom_data.add_variable("I_Lg", "amol / um^3", parameters.doubles("Initial_I_Lg"));
	cell_defaults.custom_data.add_variable("I_Rc", "amol / um^3", parameters.doubles("Initial_I_Rc"));
	cell_defaults.custom_data.add_variable("I_RL", "amol / um^3", parameters.doubles("Initial_I_RL"));
	cell_defaults.custom_data.add_variable("I_total_Rc", "amol / um^3", parameters.doubles("Initial_I_RL") + parameters.doubles("Initial_I_Rc"));

	// Total amounts - sJ, Lg, Rc
	cell_defaults.custom_data.add_variable("total_I_sJ", "amol", 0.0); // total external sJ
	cell_defaults.custom_data.add_variable("total_E_sJ", "amol", 0.0); // total external sJ

	cell_defaults.custom_data.add_variable("total_I_Lg", "amol", 0.0); 
	cell_defaults.custom_data.add_variable("total_I_Rc", "amol", 0.0); 

	cell_defaults.custom_data.add_variable("total_sJ", "amol", 0.0);   // total sJ = I_sJ + E_sJ

	// Flux - sJ
	cell_defaults.custom_data.add_variable("sJ_flux", "amol/min", 0.0);
	cell_defaults.custom_data.add_variable("sJ_adjusted_flux", "amol/min", 0.0);

	// Diffusion coefficients, concentration gradients
	cell_defaults.custom_data.add_variable("DC_sJ", "um^3/min", microenvironment.diffusion_coefficients[sJ_index]); // Diffusion coefficient
	cell_defaults.custom_data.add_variable("D_sJ", "mmol/um^3", 0.0); // Concentration gradient

	// Kinetic parameters
	cell_defaults.custom_data.add_variable("RL_binding_component", "mmol/um^3", 0.0);
	cell_defaults.custom_data.add_variable("RL_unbinding_component", "mmol/um^3", 0.0);


	// Other constants (Num. of vocels, Initial conditions, etc.)
	cell_defaults.custom_data.add_variable("Initial_E_sJ", "amol/um^3", default_microenvironment_options.initial_condition_vector[sJ_index]); // Initial external substrate_Y
	cell_defaults.custom_data.add_variable("Initial_I_sJ", "amol/um^3", parameters.doubles("Initial_I_sJ"));	
	

	cell_defaults.custom_data.add_variable("total_healthy_cells", "dimensionless", 0.0);				
	cell_defaults.custom_data.add_variable("total_cancer_cells", "dimensionless", 0.0);				

	cell_defaults.custom_data.add_variable("number_of_voxels", "dimensionless", microenvironment.number_of_voxels());
	cell_defaults.custom_data.add_variable("number_of_cells", "dimensionless", (*all_cells).size());
	cell_defaults.custom_data.add_variable("time", "min", PhysiCell_globals.current_time );
	cell_defaults.custom_data.add_variable("hill_cooperativity", "min", 0.0 );
	cell_defaults.custom_data.add_variable("gompertz_function", "min", 0.0 );
	cell_defaults.custom_data.add_variable("hill_function", "min", 0.0 );


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
void fd_channel_ligand_model(Cell *pCell, Phenotype &phenotype, double dt)
{
	Cell_Definition* pCD = find_cell_definition( pCell->type_name );

	// live.phases[0].entry_function = my_mutation_function;
	// live.phase_link(0,0).arrest_function = volume_arrest_function;
	// // live.display(std::cout);

	static int sJ_index = microenvironment.find_density_index("sJ"); // microenv substrate index

	// Densities
	static int I_sJ_index = pCell->custom_data.find_variable_index("I_sJ");
	static int E_sJ_near_index = pCell->custom_data.find_variable_index("E_sJ_near");
	static int E_sJ_index = pCell->custom_data.find_variable_index("E_sJ");
	static int I_Lg_index = pCell->custom_data.find_variable_index("I_Lg");
	static int I_Rc_index = pCell->custom_data.find_variable_index("I_Rc");
	static int I_RL_index = pCell->custom_data.find_variable_index("I_RL");
	static int I_total_Rc_index = pCell->custom_data.find_variable_index("I_total_Rc");

	// Total amounts
	double total_I_sJ_index = pCell->custom_data.find_variable_index("total_I_sJ");	
	double total_E_sJ_index = pCell->custom_data.find_variable_index("total_E_sJ"); 
	double total_E_Lg_index = pCell->custom_data.find_variable_index("total_I_Lg"); 
	double total_E_Rc_index = pCell->custom_data.find_variable_index("total_I_Rc"); 
	double total_sJ_index = pCell->custom_data.find_variable_index("total_sJ"); 
	
	// Flux data
	static int sJ_flux_index = pCell->custom_data.find_variable_index("sJ_flux");
	static int sJ_adjusted_flux_index = pCell->custom_data.find_variable_index("sJ_adjusted_flux");

	// Kinetic parameters
	static int DC_sJ_index = pCell->custom_data.find_variable_index("DC_sJ"); // Difussion coefficient
	static int D_sJ_index = pCell->custom_data.find_variable_index("D_sJ");	  // Concentration gradient
	static int RL_binding_component_index = pCell->custom_data.find_variable_index("RL_binding_component");
	static int RL_unbinding_component_index = pCell->custom_data.find_variable_index("RL_unbinding_component");
	
	// Other parameters
	static int total_healthy_cells_index = pCell->custom_data.find_variable_index("total_healthy_cells");
	static int total_cancer_cells_index = pCell->custom_data.find_variable_index("total_cancer_cells");
	static int total_cells_index = pCell->custom_data.find_variable_index("total_cells");
	pCell->custom_data.variables[total_cells_index].value = (*all_cells).size();
	
	static int time_index = pCell->custom_data.find_variable_index("time");
	pCell->custom_data.variables[time_index].value = PhysiCell_globals.current_time;

	static int hill_cooperativity_index = pCell->custom_data.find_variable_index("hill_cooperativity");
	static int gompertz_function_index = pCell->custom_data.find_variable_index("gompertz_function");
	static int hill_function_index = pCell->custom_data.find_variable_index("hill_function");



	// Agent, microenvironment mechanical parameters

	double V_cell = pCell->phenotype.volume.total;
	double R_cell = pCell->phenotype.geometry.radius;
	double A_cell = 4 * M_PI * std::pow(R_cell, 2);

	double DC_sJ = microenvironment.diffusion_coefficients[sJ_index];
	double V_voxel = (microenvironment.mesh.dV);
	double total_num_voxels = microenvironment.number_of_voxels();

	// Substrate concentrations
	double E_sJ_near = pCell->nearest_density_vector()[sJ_index];							   // A density (mM)
	double I_sJ = pCell->phenotype.molecular.internalized_total_substrates[sJ_index] / V_cell; // Convert to density (mM)
	static double Initial_E_sJ = default_microenvironment_options.initial_condition_vector[sJ_index];
	static double Initial_I_sJ = parameters.doubles("Initial_I_sJ");

	// Densities - RL system
	double I_Rc = pCell->custom_data.variables[I_Rc_index].value;
	double I_Lg = pCell->custom_data.variables[I_Lg_index].value;
	double I_RL = pCell->custom_data.variables[I_RL_index].value;
	

	// Obtaining the total net internal and external amounts through iterating over the cells
	double total_E_sJ = 0.0;
	double total_I_sJ = 0.0;
	int total_healthy_cells = 0;
	int total_cancer_cells = 0;

	std::vector<double> I_sJ_vector(microenvironment.number_of_voxels());
	

	// Counting net internal amount
	#pragma omp parallel for
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell *pC = (*all_cells)[i];
		I_sJ_vector[i] = pC->phenotype.molecular.internalized_total_substrates[sJ_index];
		total_I_sJ += pC->phenotype.molecular.internalized_total_substrates[sJ_index];

		// Get total healthy and cancer cells
		if(pC->type == 1){ total_healthy_cells +=1; }
		if(pC->type == 2){ total_cancer_cells += 1; }

	}

	double I_sJ_sum = std::accumulate(I_sJ_vector.begin(), I_sJ_vector.end(), 0.0);

	std::vector<double> E_sJ_vector(microenvironment.number_of_voxels());

	// Counting net external amount
	#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		E_sJ_vector[n] = (microenvironment.density_vector(n)[sJ_index]);
		total_E_sJ += (microenvironment.density_vector(n)[sJ_index]) * V_voxel; 

	}

	double E_sJ_sum = std::accumulate(E_sJ_vector.begin(), E_sJ_vector.end(), 0.0);
	double E_sJ_mean = E_sJ_sum / total_num_voxels; // NOT CORRECT, the accumulate just adds up the number of voxels

	
	
	double D_sJ = I_sJ - E_sJ_near; // This provides the direction of the gradient, <0 means that it goes FROM the cell
	

	/*
		Ligand-gated channel transport:
			Flux depends on the amount of RL complex, and on amount of Rc per cell.

			J = k_sJ * A_effective * D_sJ (amol/min)

	*/

	double sJ_flux;
	double sJ_adjusted_flux = sJ_flux;
	double k_sJ = parameters.doubles("k_sJ");

	double A_channel = parameters.doubles("radius_Rc");
	double A_effective = I_Rc * V_cell * A_channel;

	double RL_binding_component;
	double RL_unbinding_component;
	double RL_binding_rate = 0.0; // parameters.doubles("k_L_1");
	double RL_unbinding_rate = parameters.doubles("k_L_n1");

	double hill_index = parameters.doubles("hill_index");
	double stoichiometry = parameters.ints("RL_stoichiometry");
	double RL_Kd = ( std::pow(I_Lg, hill_index) * I_Rc ) / I_RL;

	// Fraction of Rc bound by Lg

	double hill_cooperativity = ( std::pow(I_Lg, hill_index) ) / ( 0.5 + (std::pow(I_Lg, hill_index)) );

	double gompertz_function = RL_binding_rate * std::exp(- std::exp( hill_index - hill_index * PhysiCell_globals.current_time));
	// std::cout << "Hill cooperativity: " << hill_cooperativity << std::endl;
	// std::cout << "gompertz_function: " << gompertz_function << std::endl;

	// step 1: Ligand binding to Receptor
	
	RL_binding_rate += parameters.doubles("k_L_1") * hill_cooperativity;
	RL_binding_component = RL_binding_rate * I_Rc * std::pow(I_Lg, stoichiometry);

	I_Rc -= RL_binding_component;
	I_Lg -= RL_binding_component;
	I_RL += RL_binding_component;

	if(I_Rc < 0){ I_Rc = 0; }
	if(I_Lg < 0){ I_Lg = 0; }

	// step 2: Ligand unbinds from Receptor
	RL_unbinding_component = RL_unbinding_rate * I_RL;
	I_Rc += RL_unbinding_component;
	I_Lg += RL_unbinding_component;
	I_RL -= RL_unbinding_component;
	if(I_RL < 0.0){ I_RL = 0.0; }

	// step 3: Compute flux based on the amount of RL complex
	if (I_RL > parameters.doubles("Open_threshold"))
	{
		sJ_flux = k_sJ * A_effective * D_sJ;

		// Sanity checks
		if ( sJ_flux < 0.0 ) { // Uptake
			sJ_adjusted_flux = - std::min( std::abs(sJ_flux), std::abs(E_sJ_near * V_cell) );
		}
		else if ( sJ_flux > 0.0 ) { // Secretion
			sJ_adjusted_flux = std::min(sJ_flux, std::abs(I_sJ * V_cell ) );
		}
	}
	else if (I_RL == 0.0)
	{
		sJ_flux = 0.0; // Therefore adjusted is 0.0 too
	}

	pCell->custom_data.variables[sJ_flux_index].value = sJ_flux * diffusion_dt; 
	pCell->custom_data.variables[sJ_adjusted_flux_index].value = sJ_adjusted_flux * diffusion_dt;

	pCell->phenotype.secretion.net_export_rates[sJ_index] = sJ_adjusted_flux; 

	

	/* 
	
		replace hard-coded initial values at each time-step with the new values

	*/


	// Densities - sJ
	pCell->custom_data.variables[E_sJ_near_index].value = E_sJ_near; // Mean of near_density_vector
	pCell->custom_data.variables[E_sJ_index].value = E_sJ_mean;		 // Mean of all voxels
	pCell->custom_data.variables[I_sJ_index].value = I_sJ;

	// Densities - RL system
	pCell->custom_data.variables[I_Rc_index].value = I_Rc;
	pCell->custom_data.variables[I_Lg_index].value = I_Lg;
	pCell->custom_data.variables[I_RL_index].value = I_RL;
	pCell->custom_data.variables[I_total_Rc_index].value = I_RL + I_Rc;

	// Kinetics
	pCell->custom_data.variables[DC_sJ_index].value = DC_sJ; 
	pCell->custom_data.variables[D_sJ_index].value = D_sJ;

	pCell->custom_data.variables[RL_binding_component_index].value = RL_binding_component;
	pCell->custom_data.variables[RL_unbinding_component_index].value = RL_unbinding_component;
	

	// Net total amounts
	pCell->custom_data.variables[total_E_sJ_index].value = total_E_sJ;
	pCell->custom_data.variables[total_I_sJ_index].value = total_I_sJ;
	pCell->custom_data.variables[total_sJ_index].value = total_E_sJ + total_I_sJ;

	// Other constants
	pCell->custom_data.variables[total_cancer_cells_index].value = total_cancer_cells;
	pCell->custom_data.variables[total_healthy_cells_index].value = total_healthy_cells;
	
	pCell->custom_data.variables[hill_cooperativity_index].value = hill_cooperativity;
	pCell->custom_data.variables[gompertz_function_index].value = gompertz_function;


	// Do not update phenotype if dead
	if (pCell->phenotype.death.dead == true)
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


void fd_channel_ligand_model_main( double dt )
{
	#pragma omp parallel for 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		Cell* pC = (*all_cells)[i]; 
		if( pC->phenotype.death.dead == false )
		{ fd_channel_ligand_model( pC, pC->phenotype , dt ); }
	}
	
	return;
}

