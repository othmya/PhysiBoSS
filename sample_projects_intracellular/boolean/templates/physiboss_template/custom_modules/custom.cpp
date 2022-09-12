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

	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; // is already called by standard_update_cell_velocity
	cell_defaults.functions.update_phenotype = simple_diffusion_model; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 

	// add here custom variables (same as adding them on the XML file) - for analysis with pyMCDS
	// ??? Maybe better to add them within the XML file, for easier manipulation???

	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0 ); //for paraview visualization


	// Naming convention -  I: Internal; E: External
	// !!! Keep the same substrate names as defined in the microenvironment

	static int prot_A_index = microenvironment.find_density_index("prot_A");
	Cell * pCell;


	cell_defaults.custom_data.add_variable("I_prot_A", "amol / um^3", 0.0);
	cell_defaults.custom_data.add_variable("E_prot_A_near", "amol / um^3", 0.0); // Mean nearest substrate conc.
	cell_defaults.custom_data.add_variable("E_prot_A", "amol/um^3", 0.0 ); // Mean external substrate conc.
	cell_defaults.custom_data.add_variable("total_E_prot_A", "amol", 0.0 ); // total external prot_A
	cell_defaults.custom_data.add_variable("total_I_prot_A", "amol", 0.0 ); // total internal prot_A


	cell_defaults.custom_data.add_variable("prot_A_flux", "amol/min", 0.0 ); 
	cell_defaults.custom_data.add_variable("adjusted_prot_A_flux", "amol/min", 0.0 ); 

	cell_defaults.custom_data.add_variable("D_prot_A", "mmol/um^3", 0.0 ); // Concentration gradient
	cell_defaults.custom_data.add_variable("Initial_E_prot_A", "amol/um^3", default_microenvironment_options.initial_condition_vector[prot_A_index] ); // Initial external prot_A
	cell_defaults.custom_data.add_variable("Initial_I_prot_A", "amol/um^3", parameters.doubles("Initial_I_prot_A") ); // Initial external prot_A
	cell_defaults.custom_data.add_variable("DC_prot_A", "um^3/min", microenvironment.diffusion_coefficients[prot_A_index] ); // Diffusion coefficient
	cell_defaults.custom_data.add_variable("k_prot_A", "um^2/min", parameters.doubles("k_prot_A") );
	cell_defaults.custom_data.add_variable("total_net_E_amount_theoretical", "dimensionless", microenvironment.number_of_voxels() * microenvironment.mesh.dV );
	cell_defaults.custom_data.add_variable("total_net_I_amount_theoretical", "dimensionless", 0);
	cell_defaults.custom_data.add_variable("number_of_voxels", "dimensionless", microenvironment.number_of_voxels());





	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 


	
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

void setup_tissue( void )
{

	double Initial_I_prot_A = parameters.doubles("Initial_I_prot_A");
    static int prot_A_index = microenvironment.find_density_index("prot_A");

	// Start with only 1 cell, a simple example
	// If indicated, an internal amount of substrate can be placed inside the cells
	

	for (int i=0; i < 1; i++) 
	{
		Cell* pC;
		pC = create_cell(get_cell_definition("default")); 

		// Distribute cells randomly on Unit Circle
		std::vector<double> position = UniformOnUnitCircle();
		position *= 40.0; position += {0,0,0}; // circle: radius of 50.0, centered at 0,0,0
		pC -> assign_position( position );
		//pC -> assign_position( {0,0,0} );

		// Add initial density of substrate
		pC -> phenotype.molecular.internalized_total_substrates[prot_A_index] = Initial_I_prot_A * pC->phenotype.volume.total;
		// Input is a density, convert to a net value
	}
		
	return; 
}

/*
	Simple diffusion model 

	ODE model for simple diffusion based on Fick's Law. 
	The only experimental value that is needed is the diffusion coefficient, which is in (micron^2/min)

	TO DO:
		- Add Sumbodel information class object "glucose simple diffusion" to include the specific information

*/


void simple_diffusion_model( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	static int E_prot_A_near_index = pCell->custom_data.find_variable_index("E_prot_A_near");
	static int I_prot_A_index = pCell->custom_data.find_variable_index("I_prot_A");
	static int E_prot_A_index = pCell->custom_data.find_variable_index("E_prot_A");
	
	static int D_prot_A_index = pCell->custom_data.find_variable_index("D_prot_A");
	static int total_E_prot_A_index = pCell->custom_data.find_variable_index("total_E_prot_A");
	static int total_I_prot_A_index = pCell->custom_data.find_variable_index("total_I_prot_A");
	static int total_net_I_amount_theoretical_index = pCell->custom_data.find_variable_index("total_net_I_amount_theoretical");

	static int prot_A_flux_index = pCell->custom_data.find_variable_index("prot_A_flux");
	static int adjusted_prot_A_flux_index = pCell->custom_data.find_variable_index("adjusted_prot_A_flux");





	// Constants

	static double V_cell = pCell->phenotype.volume.total; 
	static double R_cell = pCell->phenotype.geometry.radius; // Around 8, but I don't get the difference between both
	static double R_cell_2 = std::pow(V_cell * 3.0/4.0, 1.0/3.0); // 12.3, use this one for now, check the radius in PCUG
	static double A_cell = 4 * M_PI * std::pow(R_cell, 2);

	double V_voxel = ( microenvironment.mesh.dV ) ;
	double total_num_voxels = microenvironment.number_of_voxels();

	pCell->custom_data.variables[total_net_I_amount_theoretical_index].value = V_cell * (*all_cells).size(); /// Get the theoretical maximum amount of prot_A inside the cells (this is still a density)

	
	static double prot_A_k = parameters.doubles("k_prot_A");

    static int prot_A_index = microenvironment.find_density_index("prot_A");

	double E_prot_A_near = pCell->nearest_density_vector()[prot_A_index]; // A density
	double I_prot_A = pCell->phenotype.molecular.internalized_total_substrates[prot_A_index] / V_cell; // Convert to density
	// I_prot_A just refers to one cell, not useful when more than one cell in the simulation

	double total_E_prot_A;
	double total_I_prot_A;
	
	static double Initial_E_prot_A = default_microenvironment_options.initial_condition_vector[prot_A_index];
	static double Initial_I_prot_A = parameters.doubles("Initial_I_prot_A");

	std::vector<double> I_prot_A_vector(microenvironment.number_of_voxels());

	// Counting net internal amount of prot_A
	#pragma omp parallel for 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		Cell* pC = (*all_cells)[i]; 
		I_prot_A_vector[i] = pC->phenotype.molecular.internalized_total_substrates[prot_A_index];
		total_I_prot_A += pC->phenotype.molecular.internalized_total_substrates[prot_A_index];
	}

	double I_prot_A_sum = std::accumulate(I_prot_A_vector.begin(), I_prot_A_vector.end(), 0.0);

	// I_prot_A_sum and total_I_prot_A provide the same exact result. They reflect the actual total amount of
	// internal net prot_A

	std::vector<double> E_prot_A_vector(microenvironment.number_of_voxels());


	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{ 	
		E_prot_A_vector[n] = (microenvironment.density_vector(n)[prot_A_index]);
		total_E_prot_A += (microenvironment.density_vector(n)[prot_A_index]) * V_voxel; // A net amount of prot_A

	}

	double E_prot_A_sum = std::accumulate(E_prot_A_vector.begin(), E_prot_A_vector.end(), 0.0);
	double E_prot_A_mean = E_prot_A_sum / total_num_voxels; // NOT CORRECT, the accumulate just adds up the number of voxels


	double D_prot_A = I_prot_A - E_prot_A_near;


	// replace hard-coded initial values at each time-step with the new values

	pCell->custom_data.variables[E_prot_A_near_index].value = E_prot_A_near; // Mean of near_density_vector 
	pCell->custom_data.variables[E_prot_A_index].value = E_prot_A_mean; // Mean of all voxels
	pCell->custom_data.variables[I_prot_A_index].value = I_prot_A; // Density
	pCell->custom_data.variables[D_prot_A_index].value = D_prot_A; // Gradient
	pCell->custom_data.variables[total_E_prot_A_index].value = total_E_prot_A; // Net amount of prot_A
	pCell->custom_data.variables[total_I_prot_A_index].value = total_I_prot_A; // Net amount of prot_A

	//std::cout << "E_prot_A_near = " << E_prot_A_near << std::endl;

	// Following Fick's Second Law of diffusion

	double flux_prot_A  = (A_cell) * prot_A_k * (D_prot_A); // (amol/min)
	double adjusted_flux_prot_A = flux_prot_A;

	if (prot_A_k > 4.0){  

		if ( flux_prot_A < 0.0 ) { // Uptake
			adjusted_flux_prot_A = - std::min( std::abs(flux_prot_A), std::abs(E_prot_A_mean * V_voxel * dt) );
			if (adjusted_flux_prot_A = std::abs(E_prot_A_mean * V_voxel * dt)) { 
				pCell->phenotype.secretion.net_export_rates[prot_A_index] = adjusted_flux_prot_A / dt;)
			}
			 

		}
		else if ( flux_prot_A > 0.0 ) { // Secretion
			adjusted_flux_prot_A = std::min(flux_prot_A, std::abs(I_prot_A * V_cell * dt) );
			if (adjusted_flux_prot_A = std::abs(I_prot_A * V_cell * dt)) { 
				pCell->phenotype.secretion.net_export_rates[prot_A_index] = adjusted_flux_prot_A / dt;
			} 

		}

	}

	// Additional constraint
	/*
	if (Initial_E_prot_A > Initial_I_prot_A) {
		if (I_prot_A > E_prot_A_near) { adjusted_flux_prot_A = 0; }
	} else if (Initial_I_prot_A > Initial_E_prot_A) {
		if (I_prot_A < E_prot_A_near) { adjusted_flux_prot_A = 0;  }
	}
	*/

		
	pCell->custom_data.variables[prot_A_flux_index].value = flux_prot_A * dt; // this is what happens in the solver, actually
	pCell->custom_data.variables[adjusted_prot_A_flux_index].value = adjusted_flux_prot_A * dt;

	pCell->phenotype.secretion.net_export_rates[prot_A_index] = adjusted_flux_prot_A / dt; 




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

// Both functions related to the BoSS part of the simulation, leave for now
void set_input_nodes(Cell* pCell) {}

void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt) {}


// Specific coloring function

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	std::vector< std::string > output( 4 , "black" );

	if (pCell -> type == 0) // "other" cell type
	{
		
		if ( pCell -> phenotype.death.dead == false) 
		{
			output[0] = "rgb(242,161,29)"; // 
			output[1] = "rgb(201,51,242)";
			output[2] = "rgb(255,0,0)";
			output[3] = "rgb(156, 34, 109)";

			return output;
		}

	}
	else if (pCell -> type == 1) // "default" cell type
	{
		static int o2_index = microenvironment.find_density_index( "prot_A" );
		double o2 = pCell->nearest_density_vector()[o2_index];
		double o2_int = pCell->phenotype.molecular.internalized_total_substrates[o2_index];		

		double delta = o2_int - o2; // Difference between inside and outside

		int color = delta;
		char szTempString [182];
		sprintf( szTempString , "rgb(100, %u, 100)", color);
		output[2] = "rgb(255,0,0)"; // nucleus
		output[0].assign( szTempString ); // Cytoplasm


		/*

		if (delta > 0) // Red cells indicate more Glucose inside than outside
		{
			output[0] = "rgb(242,161,29)"; // 
			output[1] = "rgb(201,51,242)";
			output[2] = "rgb(255,0,0)";
			output[3] = "rgb(156, 34, 109)";

		}
		else if ( delta < 0) // Green cells indicate more Glucose OUTSIDE than Inside
		{
			output[0] = "rgb(242,161,29)"; // 
			output[1] = "rgb(201,51,242)";
			output[2] = "rgb(0,255,0)"; // nucleus
			output[3] = "rgb(156, 34, 109)";
		}
		*/

		


		


	}
	
	return output;
}


void color_node(Cell* pCell){
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}

void diffusion_model_main( double dt )
{
	#pragma omp parallel for 
	for( int i=0; i < (*all_cells).size() ; i++ )
	{
		Cell* pC = (*all_cells)[i]; 
		if( pC->phenotype.death.dead == false )
		{ simple_diffusion_model( pC, pC->phenotype , dt ); }
	}
	
	return;
}



	// Adjusting the flux to employ any permeability coefficient value

	/*

	if (prot_A_k > 200.0){ // or std::abs(flux_prot_A) > 20e3?

		if ( flux_prot_A < 0.0 ) { // Uptake
			adjusted_flux_prot_A = - std::min( std::abs(flux_prot_A), std::abs(E_prot_A_near * V_cell * dt) );
		}
		else if ( flux_prot_A > 0.0 ) { // Secretion
			adjusted_flux_prot_A = std::min(flux_prot_A, std::abs(I_prot_A * V_cell * dt ) );
		}
	}

	else { adjusted_flux_prot_A = flux_prot_A; }
	pCell->custom_data.variables[adjusted_prot_A_flux_index].value = adjusted_flux_prot_A;

	*/



/*

########## REMINDER

cell_parameters in PCUG 9.4.2.
custom phenotype in 19.1.4.
oxygen-dependent phenotype (update_cells_and_death_parameters_O2_based) in PCUG 17.6


*/