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

void create_cell_types(void)
{
	// set the random seed
	SeedRandom(parameters.ints("random_seed"));

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
	cell_defaults.functions.update_phenotype = secondary_active_model;
	cell_defaults.functions.custom_cell_rule = NULL;

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;

	/*

		Custom variables.

	*/

	static int sY_index = microenvironment.find_density_index("sY");
	static int sK_index = microenvironment.find_density_index("sK");

	Cell *pCell;
	
	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0); // for paraview visualization

	// Densities - sY & sK
	cell_defaults.custom_data.add_variable("I_sY", "amol / um^3", parameters.doubles("Initial_I_sY"));
	cell_defaults.custom_data.add_variable("E_sY_near", "amol / um^3", default_microenvironment_options.initial_condition_vector[sY_index]); // Mean nearest substrate conc.
	cell_defaults.custom_data.add_variable("E_sY", "amol/um^3", default_microenvironment_options.initial_condition_vector[sY_index]);		 // Mean external substrate conc.

	cell_defaults.custom_data.add_variable("I_sK", "amol / um^3", parameters.doubles("Initial_I_sK"));
	cell_defaults.custom_data.add_variable("E_sK_near", "amol / um^3", default_microenvironment_options.initial_condition_vector[sK_index]); // Mean nearest substrate conc.
	cell_defaults.custom_data.add_variable("E_sK", "amol/um^3", default_microenvironment_options.initial_condition_vector[sK_index]);		 // Mean external substrate conc.


	// Receptor - just one for both sY and sK
	cell_defaults.custom_data.add_variable("Rcb_sY", "amol / um^3", parameters.doubles("Initial_Rcb_sY"));											  // Bound receptor density
	cell_defaults.custom_data.add_variable("Rcb_sK", "amol / um^3", parameters.doubles("Initial_Rcb_sK"));
	cell_defaults.custom_data.add_variable("total_Rcb", "amol / um^3", parameters.doubles("Initial_Rcb_sK") + parameters.doubles("Initial_Rcb_sY")); 
	cell_defaults.custom_data.add_variable("total_Rc", "amol / um^3", parameters.doubles("Initial_Rcb_sK") + parameters.doubles("Initial_Rcb_sY") + parameters.doubles("Initial_Rc")); 

	// Total amounts - sY & sK
	cell_defaults.custom_data.add_variable("total_E_sY", "amol", 0.0); // total external sY
	cell_defaults.custom_data.add_variable("total_I_sY", "amol", 0.0); // total internal sY
	cell_defaults.custom_data.add_variable("total_sY", "amol", 0.0);   // total internal sY

	cell_defaults.custom_data.add_variable("total_E_sK", "amol", 0.0); 
	cell_defaults.custom_data.add_variable("total_I_sK", "amol", 0.0);
	cell_defaults.custom_data.add_variable("total_sK", "amol", 0.0);  


	// Kinetics of SA model - sY goes along gradient, sK is moved in opposition to gradient direction
	cell_defaults.custom_data.add_variable("DC_sY", "um^3/min", microenvironment.diffusion_coefficients[sY_index]); // Diffusion coefficient
	cell_defaults.custom_data.add_variable("D_sY", "mmol/um^3", 0.0); // Concentration gradient

	cell_defaults.custom_data.add_variable("DC_sK", "um^3/min", microenvironment.diffusion_coefficients[sY_index]); // Diffusion coefficient
	cell_defaults.custom_data.add_variable("D_sK", "mmol/um^3", 0.0); // Concentration gradient

	cell_defaults.custom_data.add_variable("sY_binding_component", "dimensionless", 0.0); // Kinetic component of sY binding
	cell_defaults.custom_data.add_variable("sY_unbinding_component", "dimensionless", 0.0); // Kinetic component of sY binding
	cell_defaults.custom_data.add_variable("sY_internal_binding_component", "dimensionless", 0.0); // Kinetic component of sY binding
	cell_defaults.custom_data.add_variable("sY_internal_unbinding_component", "dimensionless", 0.0); // Kinetic component of sY binding

	cell_defaults.custom_data.add_variable("sK_binding_component", "dimensionless", 0.0); // Kinetic component of sK binding
	cell_defaults.custom_data.add_variable("sK_unbinding_component", "dimensionless", 0.0); // Kinetic component of sK binding
	cell_defaults.custom_data.add_variable("sK_internal_binding_component", "dimensionless", 0.0); // Kinetic component of sK binding
	cell_defaults.custom_data.add_variable("sK_internal_unbinding_component", "dimensionless", 0.0); // Kinetic component of sK binding

	// Fluxes - sY & sK
	cell_defaults.custom_data.add_variable("sY_flux", "amol/min", 0.0);
	cell_defaults.custom_data.add_variable("sK_flux", "amol/min", 0.0);

	// Initial parameters and other relevant info
	cell_defaults.custom_data.add_variable("Initial_E_sY", "amol/um^3", default_microenvironment_options.initial_condition_vector[sY_index]); // Initial external substrate_Y
	cell_defaults.custom_data.add_variable("Initial_I_sY", "amol/um^3", parameters.doubles("Initial_I_sY"));								  // Initial internal substrate_Y

	cell_defaults.custom_data.add_variable("num_cells", "dimensionless", (*all_cells).size() );
	cell_defaults.custom_data.add_variable("number_of_voxels", "dimensionless", microenvironment.number_of_voxels());
	cell_defaults.custom_data.add_variable("time", "min", PhysiCell_globals.current_time );

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

	display_cell_definitions(std::cout);

	return;
}

void setup_microenvironment(void)
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
	Cell *pCell = NULL;

	/*

		Initialize agents with specific amounts of substrate inside & outside

	*/

	static int sY_index = microenvironment.find_density_index("sY");
	static int sK_index = microenvironment.find_density_index("sK");

	double Initial_Rc = parameters.doubles("Initial_Rc");
	double Initial_Rcb_sY = parameters.doubles("Initial_Rcb_sY");
	double Initial_Rcb_sK = parameters.doubles("Initial_Rcb_sK");

	double Initial_I_sY = parameters.doubles("Initial_I_sY");
	double Initial_I_sK = parameters.doubles("Initial_I_sK");

	// Disposition of agents
	// int number_of_cells = parameters.ints("number_of_cells");
	// int just_one_cell = parameters.ints("single_cell");

	
	#pragma omp parallel for
	for (int i=0; i < 1; i++) 
	{
		Cell* pC;
		pC = create_cell(get_cell_definition("default")); 
		pC -> assign_position( {0,0,0} );

		// Add initial density of substrate -- Input is a density, convert to a net value
		pC -> phenotype.molecular.internalized_total_substrates[sY_index] = Initial_I_sY * pC->phenotype.volume.total;
		pC -> phenotype.molecular.internalized_total_substrates[sK_index] = Initial_I_sK * pC->phenotype.volume.total;

	}
	
	return;

}


/*
	Secondary active transport model.
	Built as a combination of 2 coupled FDC models, a total of 8 ODEs.

*/

void secondary_active_model(Cell *pCell, Phenotype &phenotype, double dt)
{

	static int sY_index = microenvironment.find_density_index("sY"); // microenv substrate index
	static int sK_index = microenvironment.find_density_index("sK"); // microenv substrate index

	/*

		First, declare indexes for custom data which we will overwrite afterwards.
	
	*/


	// Densities - sY & sK
	static int E_sY_near_index = pCell->custom_data.find_variable_index("E_sY_near");
	static int I_sY_index = pCell->custom_data.find_variable_index("I_sY");
	static int E_sY_index = pCell->custom_data.find_variable_index("E_sY");

	static int E_sK_near_index = pCell->custom_data.find_variable_index("E_sK_near");
	static int I_sK_index = pCell->custom_data.find_variable_index("I_sK");
	static int E_sK_index = pCell->custom_data.find_variable_index("E_sK");

	// Receptor densities 
	double total_Rc_index = pCell->custom_data.find_variable_index("total_Rc"); // [E] Receptor alone (mM)
	double total_Rcb_index = pCell->custom_data.find_variable_index("total_Rcb");
	double Rcb_sY_index = pCell->custom_data.find_variable_index("Rcb_sY");
	double Rcb_sK_index = pCell->custom_data.find_variable_index("Rcb_sK");

	// Total net amounts
	static int total_E_sY_index = pCell->custom_data.find_variable_index("total_E_sY");
	static int total_I_sY_index = pCell->custom_data.find_variable_index("total_I_sY");
	static int total_sY_index = pCell->custom_data.find_variable_index("total_sY");

	static int total_E_sK_index = pCell->custom_data.find_variable_index("total_E_sK");
	static int total_I_sK_index = pCell->custom_data.find_variable_index("total_I_sK");
	static int total_sK_index = pCell->custom_data.find_variable_index("total_sK");

	// Kinetics
	static int DC_sY_index = pCell->custom_data.find_variable_index("DC_sY"); // Difussion coefficient
	static int D_sY_index = pCell->custom_data.find_variable_index("D_sY");	  
	static int DC_sK_index = pCell->custom_data.find_variable_index("DC_sK"); // Difussion coefficient
	static int D_sK_index = pCell->custom_data.find_variable_index("D_sK");	

	static int sY_binding_component_index = pCell->custom_data.find_variable_index("sY_binding_component");
	static int sY_unbinding_component_index = pCell->custom_data.find_variable_index("sY_unbinding_component");
	static int sY_internal_binding_component_index = pCell->custom_data.find_variable_index("sY_internal_binding_component");
	static int sY_internal_unbinding_component_index = pCell->custom_data.find_variable_index("sY_internal_unbinding_component");

	static int sK_binding_component_index = pCell->custom_data.find_variable_index("sK_binding_component");
	static int sK_unbinding_component_index = pCell->custom_data.find_variable_index("sK_unbinding_component");
	static int sK_internal_binding_component_index = pCell->custom_data.find_variable_index("sK_internal_binding_component");
	static int sK_internal_unbinding_component_index = pCell->custom_data.find_variable_index("sK_internal_unbinding_component");

	// Fluxes
	static int sY_flux_index = pCell->custom_data.find_variable_index("sY_flux");
	static int sK_flux_index = pCell->custom_data.find_variable_index("sK_flux");


	// Other variables 
	static int num_cells_index = pCell->custom_data.find_variable_index("num_cells"); 
	pCell->custom_data.variables[num_cells_index].value = (*all_cells).size(); // Total number of cells
	static int time_index = pCell->custom_data.find_variable_index("time");
	pCell->custom_data.variables[time_index].value = PhysiCell_globals.current_time; // Current time


	/*

		2. Declare actual variables that will be employed within the model

	*/

	// Constants
	static double DC_sY = microenvironment.diffusion_coefficients[sY_index];
	static double DC_sK = microenvironment.diffusion_coefficients[sK_index];

	static double V_cell = pCell->phenotype.volume.total;
	static double R_cell = pCell->phenotype.geometry.radius;
	static double A_cell = 4 * M_PI * std::pow(R_cell, 2);

	double V_voxel = (microenvironment.mesh.dV);
	double total_num_voxels = microenvironment.number_of_voxels();

	// Substrate concentrations
	double E_sY_near = pCell->nearest_density_vector()[sY_index];							   // A density (mM)
	double I_sY = pCell->phenotype.molecular.internalized_total_substrates[sY_index] / V_cell; // Convert to density (mM)
	
	double E_sK_near = pCell->nearest_density_vector()[sK_index];							   // A density (mM)
	double I_sK = pCell->phenotype.molecular.internalized_total_substrates[sK_index] / V_cell; // Convert to density (mM)


	// Receptor concentrations
	double total_Rc = pCell->custom_data.variables[total_Rc_index].value;	  // [E] Receptor alone (mM)
	double Rcb_sY = pCell->custom_data.variables[Rcb_sY_index].value; // [ES] Bound receptor (mM)
	double Rcb_sK = pCell->custom_data.variables[Rcb_sK_index].value; // [ES] Bound receptor (mM)
	double total_Rcb = pCell->custom_data.variables[total_Rcb_index].value; // [ES] Bound receptor (mM)
	double total_Rc_alone = total_Rc - total_Rcb;

	// Initial amounts
	static double Initial_E_sY = default_microenvironment_options.initial_condition_vector[sY_index];
	static double Initial_I_sY = parameters.doubles("Initial_I_sY"); 


	// Obtaining the total net internal and external amounts through iterating over the cells
	double total_E_sY;
	double total_I_sY;
	double total_E_sK;
	double total_I_sK;
	std::vector<double> I_sY_vector(microenvironment.number_of_voxels());
	std::vector<double> I_sK_vector(microenvironment.number_of_voxels());
	std::vector<double> E_sY_vector(microenvironment.number_of_voxels());
	std::vector<double> E_sK_vector(microenvironment.number_of_voxels());


	#pragma omp parallel for
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell *pC = (*all_cells)[i];
		I_sY_vector[i] = pC->phenotype.molecular.internalized_total_substrates[sY_index];
		total_I_sY += pC->phenotype.molecular.internalized_total_substrates[sY_index];

		I_sK_vector[i] = pC->phenotype.molecular.internalized_total_substrates[sK_index];
		total_I_sK += pC->phenotype.molecular.internalized_total_substrates[sK_index];
	}
	#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		E_sY_vector[n] = (microenvironment.density_vector(n)[sY_index]);
		E_sK_vector[n] = (microenvironment.density_vector(n)[sK_index]);
		total_E_sY += (microenvironment.density_vector(n)[sY_index]) * V_voxel;
		total_E_sK += (microenvironment.density_vector(n)[sK_index]) * V_voxel;
	}

	double I_sY_sum = std::accumulate(I_sY_vector.begin(), I_sY_vector.end(), 0.0);
	double I_sK_sum = std::accumulate(I_sK_vector.begin(), I_sK_vector.end(), 0.0);

	double E_sY_sum = std::accumulate(E_sY_vector.begin(), E_sY_vector.end(), 0.0);
	double E_sK_sum = std::accumulate(E_sK_vector.begin(), E_sK_vector.end(), 0.0);

	double E_sY_mean = E_sY_sum / total_num_voxels; 
	double E_sK_mean = E_sK_sum / total_num_voxels; 

	double D_sY = I_sY - E_sY_near;
	double D_sK = I_sK - E_sK_near;

	/*

		3. Compute net export rate as a result of a Secondary Active transport model.

	*/

	// KKinetic rates for the SA model
	static double sY_binding_rate = parameters.doubles("sY_k_1A");
	static double sY_unbinding_rate = parameters.doubles("sY_k_n1A");
	static double sY_internal_binding_rate = parameters.doubles("sY_k_2A");
	static double sY_internal_unbinding_rate = parameters.doubles("sY_k_n2A");

	static double sK_binding_rate = parameters.doubles("sK_k_1B");
	static double sK_unbinding_rate = parameters.doubles("sK_k_n1B");
	static double sK_internal_binding_rate = parameters.doubles("sK_k_2B");
	static double sK_internal_unbinding_rate = parameters.doubles("sK_k_n2B");


	// Separated components

	double sY_binding_component; // k1A
	double sY_unbinding_component; // kn1A
	double sY_internal_binding_component; // k2A
	double sY_internal_unbinding_component; // kn2A

	double sK_binding_component; // k1B
	double sK_unbinding_component; // kn1B
	double sK_internal_binding_component; // k2B
	double sK_internal_unbinding_component; // kn2B

	std::cout << "I_sY: " << I_sY << std::endl;
	std::cout << "I_sK: " << I_sK << std::endl;


	// step 1: sY binding from outside
	sY_binding_component = diffusion_dt * sY_binding_rate * (total_Rc - total_Rcb) * E_sY_near;

	if (sY_binding_component > total_Rc * diffusion_dt || sY_binding_component > E_sY_near * diffusion_dt )
	{ sY_binding_component = std::min(total_Rc * diffusion_dt, E_sY_near * diffusion_dt); } // Get the limiting step 
	
	Rcb_sY += sY_binding_component;
	// total_Rc -= sY_binding_component;
	if(total_Rc < 0.0) { total_Rc = 0.0; }

	// step 2: some of the sY unbinds
	sY_unbinding_component = diffusion_dt * sY_unbinding_rate * Rcb_sY;
	if (sY_unbinding_component > Rcb_sY * diffusion_dt){ sY_unbinding_component = Rcb_sY * diffusion_dt; }

	// total_Rc += sY_unbinding_component;
	Rcb_sY -= sY_unbinding_component;
	if(Rcb_sY < 0.0) { Rcb_sY = 0.0; }
	if(total_Rc > parameters.doubles("Initial_Rc")) { total_Rc = parameters.doubles("Initial_Rc"); }

	// step 3: sK binds from inside
	sK_internal_binding_component = diffusion_dt * sK_internal_binding_rate *  (total_Rc - total_Rcb)  * I_sK;
	if (sK_internal_binding_component > total_Rc * diffusion_dt || sK_internal_binding_component > I_sK * diffusion_dt )
	{ sK_internal_binding_component = std::min(total_Rc  * diffusion_dt, I_sK * diffusion_dt); } // Get the limiting step 
	
	I_sK -= sK_internal_binding_component;
	Rcb_sK += sK_internal_binding_component;
	// total_Rc -= sK_internal_binding_component;

	if(total_Rc < 0.0) { total_Rc = 0.0; }
	if(I_sK < 0.0) { I_sK = 0.0; }

	// step 4: some of the sK unbinds from inside
	sK_internal_unbinding_component = diffusion_dt * sK_internal_unbinding_rate * Rcb_sK;
	if (sK_internal_unbinding_component > Rcb_sK * diffusion_dt){ sK_internal_unbinding_component = Rcb_sK * diffusion_dt; } 

	// total_Rc += sK_internal_unbinding_component;
	Rcb_sK -= sK_internal_unbinding_component;
	if(Rcb_sK < 0.0) { Rcb_sK = 0.0; }
	if(total_Rc > parameters.doubles("Initial_Rc")) { total_Rc = parameters.doubles("Initial_Rc"); }


	// step 5: sY enters the agent
	sY_internal_unbinding_component = diffusion_dt * sY_internal_unbinding_rate * Rcb_sY;
	if (sY_internal_unbinding_component > Rcb_sY * diffusion_dt){ sY_internal_unbinding_component = Rcb_sY * diffusion_dt; } 

	I_sY += sY_internal_unbinding_component;
	Rcb_sY -= sY_internal_unbinding_component;

	// Then, every sY that goes in, sK goes out
	
	// Antiporting
	// I_sK -= sY_internal_unbinding_component;
	// E_sK_near += sY_internal_unbinding_component;
	// Rcb_sK -= sK_unbinding_component;

	// Symporting
	I_sK += sY_internal_unbinding_component;
	E_sK_near -= sY_internal_unbinding_component;
	Rcb_sK += sK_unbinding_component;


	// total_Rc += sY_internal_unbinding_component;
	if(E_sK_near < 0.0) { E_sK_near = 0.0; }
	if(Rcb_sY < 0.0) { Rcb_sY = 0.0; }
	if(I_sK < 0.0) { I_sK = 0.0; }
	// if(total_Rc > parameters.doubles("Initial_Rc")) { total_Rc = parameters.doubles("Initial_Rc"); }

	// // step 6: sK is pumped outside the agent - against gradient
	// sK_unbinding_component = diffusion_dt * sY_internal_unbinding_rate * Rcb_sK; 
	// if (sK_unbinding_component > Rcb_sK * diffusion_dt){ sK_unbinding_component = Rcb_sK * diffusion_dt; }

	// I_sK -= sK_unbinding_component;
	// E_sK_near += sK_unbinding_component;
	// Rcb_sK -= sK_unbinding_component;
	// // total_Rc += sK_unbinding_component;
	// if(Rcb_sK < 0.0) { Rcb_sK = 0.0; }
	// if(I_sK < 0.0) { I_sK = 0.0; }
	// if(total_Rc > parameters.doubles("Initial_Rc")) { total_Rc = parameters.doubles("Initial_Rc"); }


	// step 7: some of sY binds inside
	sY_internal_binding_component = diffusion_dt * sY_internal_binding_rate *  (total_Rc - total_Rcb)  * I_sY;
	if (sY_internal_binding_component > Rcb_sY * diffusion_dt || sY_internal_binding_component > I_sY * diffusion_dt )
	{ sY_internal_binding_component = std::min(Rcb_sY  * diffusion_dt, I_sY * diffusion_dt); } // Get the limiting step 
	
	I_sY -= sY_internal_binding_component;
	Rcb_sY += sY_internal_binding_component;
	// total_Rcb -= sY_internal_binding_component;

	if(total_Rc < 0.0) { total_Rc = 0.0; }
	if(I_sY < 0.0) { I_sY = 0.0; }

	// step 8: some of sK binds outside
	sK_binding_component = diffusion_dt * sK_binding_rate *  (total_Rc - total_Rcb)  * E_sK_near;
	if (sK_binding_component > total_Rc * diffusion_dt || sK_binding_component > E_sK_near * diffusion_dt )
	{ sK_binding_component = std::min(total_Rc  * diffusion_dt, E_sK_near * diffusion_dt); }

	E_sK_near -= sK_binding_component;
	Rcb_sK += sK_binding_component;
	// total_Rcb -= sK_binding_component;
	if(total_Rc < 0.0) { total_Rc = 0.0; }
	if(E_sK_near < 0.0) { E_sK_near = 0.0; }


	// step 9: Compute reaction velocities for both substrates
	double sY_vf_max = sY_internal_unbinding_rate * (total_Rc); // both are constant given that total receptor is constant
	double sY_vb_max = sY_internal_binding_rate * (total_Rc);
	double sY_Km_1 = (sY_internal_binding_rate + sY_internal_unbinding_rate) / sY_binding_rate;
	double sY_Km_2 = (sY_internal_binding_rate + sY_internal_unbinding_rate) / sY_internal_binding_rate;

	double sY_monod = sY_vf_max * (I_sY / (sY_Km_1 + I_sY)); // Using the "forward" (entry) kinetic values
	sY_monod *= diffusion_dt;

	double sY_velocity = -((sY_vf_max * E_sY_near) / sY_Km_1 - (sY_vb_max * I_sY) / sY_Km_2) / (1 + (I_sY / sY_Km_2) + (E_sY_near / sY_Km_1));
	sY_velocity *= V_cell;

	double sK_vf_max = sK_internal_unbinding_rate * (total_Rc); // both are constant given that total receptor is constant
	double sK_vb_max = sK_internal_binding_rate * (total_Rc);
	double sK_Km_1 = (sK_internal_binding_rate + sK_internal_unbinding_rate) / sK_binding_rate;
	double sK_Km_2 = (sK_internal_binding_rate + sK_internal_unbinding_rate) / sK_internal_binding_rate;
	double sK_monod = sK_vf_max * (I_sK / (sK_Km_1 + I_sK));
	sK_monod *= diffusion_dt;



	double sK_velocity = -sY_velocity;
	
	// if(std::abs(sK_velocity / V_cell) > I_sY || Rcb_sY = 0.0 || std::abs(sK_velocity / V_cell) > E_sY_near){
	// 	sK_velocity = 0.0;
	// 	sY_velocity = 0.0;
	// }



	/*

		4. Replace hard-coded initial values at each time-step with the new values

	*/

	// Substrate concentrations
	
	pCell->custom_data.variables[E_sY_near_index].value = E_sY_near; // Mean of near_density_vector
	pCell->custom_data.variables[E_sY_index].value = E_sY_mean;		 // Mean of all voxels
	pCell->custom_data.variables[I_sY_index].value = I_sY;			 // Density
	pCell->custom_data.variables[E_sK_near_index].value = E_sK_near; // Mean of near_density_vector
	pCell->custom_data.variables[E_sK_index].value = E_sK_mean;		 // Mean of all voxels
	pCell->custom_data.variables[I_sK_index].value = I_sK;			 // Density

	// Receptor concentrations

	pCell->custom_data.variables[total_Rc_index].value = total_Rc;
	pCell->custom_data.variables[total_Rcb_index].value = Rcb_sY + Rcb_sK;
	pCell->custom_data.variables[Rcb_sK_index].value = Rcb_sK;
	pCell->custom_data.variables[Rcb_sY_index].value = Rcb_sY;

	// Kinetics
	pCell->custom_data.variables[DC_sY_index].value = DC_sY; 
	pCell->custom_data.variables[D_sY_index].value = D_sY;	 
	pCell->custom_data.variables[D_sK_index].value = D_sK;	 

	pCell->custom_data.variables[sY_binding_component_index].value = sY_binding_component;
	pCell->custom_data.variables[sY_internal_binding_component_index].value = sY_internal_binding_component;
	pCell->custom_data.variables[sY_internal_unbinding_component_index].value = sY_internal_unbinding_component;
	pCell->custom_data.variables[sY_unbinding_component_index].value = sY_unbinding_component;

	pCell->custom_data.variables[sK_binding_component_index].value = sK_binding_component;
	pCell->custom_data.variables[sK_internal_binding_component_index].value = sK_internal_binding_component;
	pCell->custom_data.variables[sK_internal_unbinding_component_index].value = sK_internal_unbinding_component;
	pCell->custom_data.variables[sK_unbinding_component_index].value = sK_unbinding_component;

	// Fluxes
	pCell->custom_data.variables[sY_flux_index].value = sY_velocity;
	pCell->custom_data.variables[sK_flux_index].value = sK_velocity;
	
	// Net total amounts
	pCell->custom_data.variables[total_E_sY_index].value = E_sY_sum;			 
	pCell->custom_data.variables[total_I_sY_index].value = I_sY_sum;			  
	pCell->custom_data.variables[total_sY_index].value = total_E_sY + total_I_sY; 
	pCell->custom_data.variables[total_E_sK_index].value = E_sK_sum;			 
	pCell->custom_data.variables[total_I_sK_index].value = I_sK_sum;			  
	pCell->custom_data.variables[total_sK_index].value = total_E_sK + total_I_sK; 


	/*

		5. Connecting the model with different aspects of the phenotype

	*/

	pCell->phenotype.secretion.net_export_rates[sY_index] = sY_velocity;
	pCell->phenotype.secretion.net_export_rates[sK_index] = sK_velocity;

	static int necrosis_idx = pCell->phenotype.death.find_death_model_index("Necrosis");
	static int apoptosis_idx = pCell->phenotype.death.find_death_model_index("Apoptosis");

	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 

	// pCell->phenotype.cycle.data.transition_rate(0,0) = sY_monod * V_cell;


	
	// Do not update phenotype if dead
	if (pCell->phenotype.death.dead == true)
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
			pCell->type == get_cell_definition("default").type && PhysiCell::PhysiCell_globals.current_time >= 0 && pCell->phenotype.intracellular->get_parameter_value("$time_scale") == 0.0)
			pCell->phenotype.intracellular->set_parameter_value("$time_scale", 0.0);

		set_input_nodes(pCell);

		pCell->phenotype.intracellular->update();

		from_nodes_to_cell(pCell, phenotype, dt);
		color_node(pCell);
	}

	return;
}

// Both functions related to the BoSS part of the simulation, leave for now
void set_input_nodes(Cell *pCell) {}

void from_nodes_to_cell(Cell *pCell, Phenotype &phenotype, double dt) {}

// Specific coloring function

std::vector<std::string> my_coloring_function(Cell *pCell)
{
	std::vector<std::string> output = false_cell_coloring_live_dead(pCell);

	if (pCell->type == 1) // "default" cell type
	{
		static int sY_index = microenvironment.find_density_index("sY");
		double E_sY_near = pCell->nearest_density_vector()[sY_index];
		double I_sY = pCell->phenotype.molecular.internalized_total_substrates[sY_index];

		double delta = I_sY - E_sY_near; // Difference between inside and outside

		int color = delta;
		char szTempString[182];
		sprintf(szTempString, "rgb(100, %u, 100)", color);
		output[2] = "rgb(255,0,0)";		// nucleus
		output[0].assign(szTempString); // Cytoplasm
	}

	return output;
}

void color_node(Cell *pCell)
{
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}

void secondary_active_model_main(double dt)
{
#pragma omp parallel for
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell *pC = (*all_cells)[i];
		if (pC->phenotype.death.dead == false)
		{
			secondary_active_model(pC, pC->phenotype, dt);
		}
	}

	return;
}


// for (int i = 0; i < parameters.ints("cells_number"); i++)
// {
// 	Cell *pC;
// 	pC = create_cell(get_cell_definition("default"));

// 	// Distribute cells randomly on Unit Circle
// 	std::vector<double> position = UniformOnUnitCircle();
// 	position *= 80.0;
// 	position += {0, 0, 0}; // circle: radius of 50.0, centered at 0,0,0
// 	pC->assign_position(position);
// 	// pC -> assign_position( {0,0,0} );

// 	// Add initial density of substrate
// 	pC->phenotype.molecular.internalized_total_substrates[sY_index] = Initial_I_sY * pC->phenotype.volume.total;
// 	// Input is a density, convert to a net value
// }