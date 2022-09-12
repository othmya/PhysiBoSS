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
	cell_defaults.functions.update_phenotype = drug_resistance_model;
	cell_defaults.functions.custom_cell_rule = NULL;

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;

	/*

		Custom variables.

	*/

	static int sY_index = microenvironment.find_density_index("sY");
	static int cytB_index = microenvironment.find_density_index("cytB");

	Cell *pCell;
	
	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0); // for paraview visualization

	// Densities - sY, ATP-ADP, cytB
	cell_defaults.custom_data.add_variable("I_sY", "amol / um^3", parameters.doubles("Initial_I_sY"));
	cell_defaults.custom_data.add_variable("E_sY_near", "amol / um^3", default_microenvironment_options.initial_condition_vector[sY_index]); // Mean nearest substrate conc.
	cell_defaults.custom_data.add_variable("E_sY", "amol/um^3", default_microenvironment_options.initial_condition_vector[sY_index]);		 // Mean external substrate conc.

	cell_defaults.custom_data.add_variable("I_cytB", "amol / um^3", parameters.doubles("Initial_I_cytB"));
	cell_defaults.custom_data.add_variable("E_cytB_near", "amol / um^3", default_microenvironment_options.initial_condition_vector[cytB_index]); // Mean nearest substrate conc.
	cell_defaults.custom_data.add_variable("E_cytB", "amol/um^3", default_microenvironment_options.initial_condition_vector[cytB_index]);		 // Mean external substrate conc.

	cell_defaults.custom_data.add_variable("I_ATP", "amol / um^3", parameters.doubles("Initial_I_ATP"));
	cell_defaults.custom_data.add_variable("I_ADP", "amol / um^3", parameters.doubles("Initial_I_ADP"));

	// Receptors - sY FDC, then PA for cytB 
	cell_defaults.custom_data.add_variable("Rcb_sY", "amol / um^3", parameters.doubles("Initial_Rcb_sY"));											  // Bound receptor density
	cell_defaults.custom_data.add_variable("Rc_sY", "amol / um^3", parameters.doubles("Initial_Rc_sY")); 
	cell_defaults.custom_data.add_variable("total_Rc_sY", "amol / um^3", parameters.doubles("Initial_Rc_sY") + parameters.doubles("Initial_Rcb_sY")); 

	cell_defaults.custom_data.add_variable("Rcb_cytB", "amol / um^3", parameters.doubles("Initial_Rcb_cytB"));											  // Bound receptor density
	cell_defaults.custom_data.add_variable("Rc_cytB", "amol / um^3", parameters.doubles("Initial_Rc_cytB"));											  // Bound receptor density
	cell_defaults.custom_data.add_variable("total_Rc_cytB", "amol / um^3", parameters.doubles("Initial_Rcb_cytB") + parameters.doubles("Initial_Rc_cytB") );  // Bound receptor density

	// Total amounts - sY & cytB
	cell_defaults.custom_data.add_variable("total_E_sY", "amol", 0.0); // total external sY
	cell_defaults.custom_data.add_variable("total_I_sY", "amol", 0.0); // total internal sY
	cell_defaults.custom_data.add_variable("total_sY", "amol", 0.0);   // total internal sY

	cell_defaults.custom_data.add_variable("total_E_cytB", "amol", 0.0); 
	cell_defaults.custom_data.add_variable("total_I_cytB", "amol", 0.0);
	cell_defaults.custom_data.add_variable("total_cytB", "amol", 0.0);  

	// Kinetics of sY -- FDC model
	cell_defaults.custom_data.add_variable("DC_sY", "um^3/min", microenvironment.diffusion_coefficients[sY_index]); // Diffusion coefficient
	cell_defaults.custom_data.add_variable("D_sY", "mmol/um^3", 0.0); // Concentration gradient

	cell_defaults.custom_data.add_variable("sY_binding_component", "mM/min", 0.0); 
	cell_defaults.custom_data.add_variable("sY_unbinding_component", "mM/min", 0.0); 
	cell_defaults.custom_data.add_variable("sY_internal_binding_component", "mM/min", 0.0); 
	cell_defaults.custom_data.add_variable("sY_movement_component", "mM/min", 0.0); 


	// Kinetics of cytB -- SD model
	cell_defaults.custom_data.add_variable("DC_cytB", "um^3/min", microenvironment.diffusion_coefficients[sY_index]); // Diffusion coefficient
	cell_defaults.custom_data.add_variable("D_cytB", "mmol/um^3", 0.0); // Concentration gradient
	cell_defaults.custom_data.add_variable("k_cytB", "mmol/um^3", 0.0); // Permeability coefficient


	// Kinetics of PA bomb of cytB	
	cell_defaults.custom_data.add_variable("cytB_binding_component", "dimensionless", 0.0); 
	cell_defaults.custom_data.add_variable("cytB_unbinding_component", "dimensionless", 0.0); 
	cell_defaults.custom_data.add_variable("cytB_internal_binding_component", "dimensionless", 0.0);
	cell_defaults.custom_data.add_variable("cytB_movement_component", "dimensionless", 0.0);

	cell_defaults.custom_data.add_variable("ATP_binding_component", "dimensionless", 0.0);
	cell_defaults.custom_data.add_variable("ATP_hydrolisis_component", "dimensionless", 0.0);
	cell_defaults.custom_data.add_variable("ATP_regeneration_component", "dimensionless", 0.0);

	// Net fluxes - sY & cytB
	cell_defaults.custom_data.add_variable("sY_flux", "amol/min", 0.0);
	cell_defaults.custom_data.add_variable("cytB_flux", "amol/min", 0.0);

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
	static int cytB_index = microenvironment.find_density_index("cytB");

	// double Initial_Rc_sY = parameters.doubles("Initial_Rc_sY");
	// double Initial_Rc_cytB = parameters.doubles("Initial_Rc_cytB");

	// double Initial_Rcb_sY = parameters.doubles("Initial_Rcb_sY");
	// double Initial_Rcb_cytB = parameters.doubles("Initial_Rcb_cytB");

	double Initial_I_sY = parameters.doubles("Initial_I_sY");
	double Initial_I_cytB = parameters.doubles("Initial_I_cytB");
	// double Initial_I_ATP = parameters.doubles("Initial_I_ATP");
	// double Initial_I_ADP = parameters.doubles("Initial_I_ADP");

	// Disposition of agents
	// int number_of_cells = parameters.ints("number_of_cells");
	// int just_one_cell = parameters.ints("single_cell");

	
	double cell_radius = cell_defaults.phenotype.geometry.radius;
	double cell_spacing = 0.95 * 2.0 * cell_radius;
	double tumor_radius = parameters.doubles("tumor_radius");  		

	// Parameter<double> temp;

	double x = 0.0;
	double x_outer = tumor_radius;
	double y = 0.0;


	int n = 0;
	while (y < tumor_radius)
	{
		x = 0.0;
		if (n % 2 == 1)
		{
			x = 0.5 * cell_spacing;
		}
		x_outer = sqrt(tumor_radius * tumor_radius - y * y);

		while (x < x_outer)
		{
			pCell = create_cell(); // tumor cell
			pCell->assign_position(x, y, 0.0);
			pCell->phenotype.molecular.internalized_total_substrates[sY_index] = Initial_I_sY * pCell->phenotype.volume.total;
			pCell->phenotype.molecular.internalized_total_substrates[cytB_index] = Initial_I_cytB * pCell->phenotype.volume.total;

			if (fabs(y) > 0.01)
			{
				pCell = create_cell(); // tumor cell
				pCell->assign_position(x, -y, 0.0);
				pCell->phenotype.molecular.internalized_total_substrates[sY_index] = Initial_I_sY * pCell->phenotype.volume.total;
				pCell->phenotype.molecular.internalized_total_substrates[cytB_index] = Initial_I_cytB * pCell->phenotype.volume.total;

			}

			if (fabs(x) > 0.01)
			{
				pCell = create_cell(); // tumor cell
				pCell->assign_position(-x, y, 0.0);
				pCell->phenotype.molecular.internalized_total_substrates[sY_index] = Initial_I_sY * pCell->phenotype.volume.total;
				pCell->phenotype.molecular.internalized_total_substrates[cytB_index] = Initial_I_cytB * pCell->phenotype.volume.total;


				if (fabs(y) > 0.01)
				{
					pCell = create_cell(); // tumor cell
					pCell->assign_position(-x, -y, 0.0);
					pCell->phenotype.molecular.internalized_total_substrates[sY_index] = Initial_I_sY * pCell->phenotype.volume.total;
					pCell->phenotype.molecular.internalized_total_substrates[cytB_index] = Initial_I_cytB * pCell->phenotype.volume.total;
				
				}
			}
			x += cell_spacing;
		}

		y += cell_spacing * sqrt(3.0) / 2.0;
		n++;
	}
	
	return;

}


/*
	Secondary active transport model.
	Built as a combination of 2 coupled FDC models, a total of 8 ODEs.

*/

void drug_resistance_model(Cell *pCell, Phenotype &phenotype, double dt)
{

	static int sY_index = microenvironment.find_density_index("sY"); // microenv substrate index
	static int cytB_index = microenvironment.find_density_index("cytB"); // microenv substrate index

	/*

		First, declare indexes for custom data which we will overwrite afterwards.
	
	*/


	// Densities - sY & cytB
	static int E_sY_near_index = pCell->custom_data.find_variable_index("E_sY_near");
	static int I_sY_index = pCell->custom_data.find_variable_index("I_sY");
	static int E_sY_index = pCell->custom_data.find_variable_index("E_sY");

	static int E_cytB_near_index = pCell->custom_data.find_variable_index("E_cytB_near");
	static int I_cytB_index = pCell->custom_data.find_variable_index("I_cytB");
	static int E_cytB_index = pCell->custom_data.find_variable_index("E_cytB");

	// Receptor densities -- sY (FDC)
	double total_Rc_sY_index = pCell->custom_data.find_variable_index("total_Rc_sY"); // [E] Receptor alone (mM)
	double Rc_sY_index = pCell->custom_data.find_variable_index("Rc_sY");
	double Rcb_sY_index = pCell->custom_data.find_variable_index("Rcb_sY");

	// Receptor densities -- cytB (PA)

	double total_Rc_cytB_index = pCell->custom_data.find_variable_index("total_Rc_cytB");
	double total_Rcb_cytB_index = pCell->custom_data.find_variable_index("total_Rcb_cytB");
	double Rcb_cytB_index = pCell->custom_data.find_variable_index("Rcb_cytB");
	double Rc_cytB_index = pCell->custom_data.find_variable_index("Rc_cytB");

	// Total net amounts
	static int total_E_sY_index = pCell->custom_data.find_variable_index("total_E_sY");
	static int total_I_sY_index = pCell->custom_data.find_variable_index("total_I_sY");
	static int total_sY_index = pCell->custom_data.find_variable_index("total_sY");

	static int total_E_cytB_index = pCell->custom_data.find_variable_index("total_E_cytB");
	static int total_I_cytB_index = pCell->custom_data.find_variable_index("total_I_cytB");
	static int total_cytB_index = pCell->custom_data.find_variable_index("total_cytB");

	// Kinetics -- sY (FDC)
	static int DC_sY_index = pCell->custom_data.find_variable_index("DC_sY"); // Difussion coefficient
	static int D_sY_index = pCell->custom_data.find_variable_index("D_sY");	  

	static int sY_binding_component_index = pCell->custom_data.find_variable_index("sY_binding_component");
	static int sY_unbinding_component_index = pCell->custom_data.find_variable_index("sY_unbinding_component");
	static int sY_internal_binding_component_index = pCell->custom_data.find_variable_index("sY_internal_binding_component");
	static int sY_internal_unbinding_component_index = pCell->custom_data.find_variable_index("sY_internal_unbinding_component");


	// Kinetics -- cytB (PA)
	static int DC_cytB_index = pCell->custom_data.find_variable_index("DC_cytB"); // Difussion coefficient
	static int D_cytB_index = pCell->custom_data.find_variable_index("D_cytB");	


	static int cytB_binding_component_index = pCell->custom_data.find_variable_index("cytB_binding_component");
	static int cytB_unbinding_component_index = pCell->custom_data.find_variable_index("cytB_unbinding_component");
	static int cytB_internal_binding_component_index = pCell->custom_data.find_variable_index("cytB_internal_binding_component");
	static int cytB_movement_component_index = pCell->custom_data.find_variable_index("cytB_movement_component");

	static int ATP_binding_component_index = pCell->custom_data.find_variable_index("ATP_binding_component");
	static int ATP_hydrolisis_component_index = pCell->custom_data.find_variable_index("ATP_hydrolisis_component");
	static int ATP_regeneration_component_index = pCell->custom_data.find_variable_index("ATP_regeneration_component");

	// Fluxes
	static int sY_flux_index = pCell->custom_data.find_variable_index("sY_flux");
	static int cytB_flux_index = pCell->custom_data.find_variable_index("cytB_flux");


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
	static double DC_cytB = microenvironment.diffusion_coefficients[cytB_index];

	static double V_cell = pCell->phenotype.volume.total;
	static double R_cell = pCell->phenotype.geometry.radius;
	static double A_cell = 4 * M_PI * std::pow(R_cell, 2);

	double V_voxel = (microenvironment.mesh.dV);
	double total_num_voxels = microenvironment.number_of_voxels();

	// Substrate concentrations
	double E_sY_near = pCell->nearest_density_vector()[sY_index];							   // A density (mM)
	double I_sY = pCell->phenotype.molecular.internalized_total_substrates[sY_index] / V_cell; // Convert to density (mM)
	
	double E_cytB_near = pCell->nearest_density_vector()[cytB_index];							   // A density (mM)
	double I_cytB = pCell->phenotype.molecular.internalized_total_substrates[cytB_index] / V_cell; // Convert to density (mM)xx


	// Receptor concentrations
	double total_Rc_sY = pCell->custom_data.variables[total_Rc_sY_index].value;	  // [E] Receptor alone (mM)
	double Rcb_sY = pCell->custom_data.variables[Rcb_sY_index].value; // [ES] Bound receptor (mM)
	double Rc_sY = pCell->custom_data.variables[Rc_sY_index].value; // [ES] Bound receptor (mM)


	double total_Rc_cytB = pCell->custom_data.variables[total_Rc_cytB_index].value; // [ES] Bound receptor (mM)
	double Rcb_cytB = pCell->custom_data.variables[Rcb_cytB_index].value; // [ES] Bound receptor (mM)
	double Rc_cytB = pCell->custom_data.variables[Rc_cytB_index].value; // [ES] Bound receptor (mM)

	// Initial amounts
	static double Initial_E_sY = default_microenvironment_options.initial_condition_vector[sY_index];
	static double Initial_I_sY = parameters.doubles("Initial_I_sY"); 

	static double Initial_E_cytB = default_microenvironment_options.initial_condition_vector[cytB_index];
	static double Initial_I_cytB = parameters.doubles("Initial_I_cytB"); 


	// Obtaining the total net internal and external amounts through iterating over the cells
	double total_E_sY;
	double total_I_sY;
	double total_E_cytB;
	double total_I_cytB;
	std::vector<double> I_sY_vector(microenvironment.number_of_voxels());
	std::vector<double> I_cytB_vector(microenvironment.number_of_voxels());
	std::vector<double> E_sY_vector(microenvironment.number_of_voxels());
	std::vector<double> E_cytB_vector(microenvironment.number_of_voxels());


	#pragma omp parallel for
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell *pC = (*all_cells)[i];
		I_sY_vector[i] = pC->phenotype.molecular.internalized_total_substrates[sY_index];
		total_I_sY += pC->phenotype.molecular.internalized_total_substrates[sY_index];

		I_cytB_vector[i] = pC->phenotype.molecular.internalized_total_substrates[cytB_index];
		total_I_cytB += pC->phenotype.molecular.internalized_total_substrates[cytB_index];
	}
	#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		E_sY_vector[n] = (microenvironment.density_vector(n)[sY_index]);
		E_cytB_vector[n] = (microenvironment.density_vector(n)[cytB_index]);
		total_E_sY += (microenvironment.density_vector(n)[sY_index]) * V_voxel;
		total_E_cytB += (microenvironment.density_vector(n)[cytB_index]) * V_voxel;
	}

	double I_sY_sum = std::accumulate(I_sY_vector.begin(), I_sY_vector.end(), 0.0);
	double I_cytB_sum = std::accumulate(I_cytB_vector.begin(), I_cytB_vector.end(), 0.0);

	double E_sY_sum = std::accumulate(E_sY_vector.begin(), E_sY_vector.end(), 0.0);
	double E_cytB_sum = std::accumulate(E_cytB_vector.begin(), E_cytB_vector.end(), 0.0);

	double E_sY_mean = E_sY_sum / total_num_voxels; 
	double E_cytB_mean = E_cytB_sum / total_num_voxels; 

	double D_sY = I_sY - E_sY_near;
	double D_cytB = I_cytB - E_cytB_near;

	/*

		3. Compute net export rate as a result of a combination of SD, FDC and PA models

	*/

	// Kinetic rates for the FDC model
	static double sY_binding_rate = parameters.doubles("k_1");
	static double sY_unbinding_rate = parameters.doubles("k_n1");
	static double sY_movement_rate = parameters.doubles("k_2");
	static double sY_internal_binding_rate = parameters.doubles("k_n2");


	// Kinetic rates for the PA model
	static double cytB_binding_rate = parameters.doubles("cytB_k_1B");
	static double cytB_unbinding_rate = parameters.doubles("cytB_k_n1B");
	static double cytB_movement_rate = parameters.doubles("cytB_k_2B");
	static double cytB_internal_binding_rate = parameters.doubles("cytB_k_n2B");

	static double ATP_binding_rate = parameters.doubles("k_1");
	static double ATP_hydrolisis_rate = parameters.doubles("k_n1");
	static double ATP_regeneration_rate = parameters.doubles("k_n2");


	// Separated components

	double sY_binding_component; // k1A
	double sY_unbinding_component; // kn1A
	double sY_movement_component; // k2A
	double sY_internal_binding_component; // kn2A

	double cytB_binding_component; // k1B
	double cytB_unbinding_component; // kn1B
	double cytB_movement_component; // k2B
	double cytB_internal_binding_component; // kn2B

	double ATP_hydrolisis_component;
	double ADP_recycling_component;

	/* 
	
	STEP 1: Movement of cytB through simple diffusion

	*/

	double k_cytB = parameters.doubles("k_cytB");

	double cytB_flux  = k_cytB * D_cytB * A_cell; // (amol/min)
	double adjusted_cytB_flux = cytB_flux;
	

	if ( cytB_flux < 0.0 ) { // Uptake
		adjusted_cytB_flux = - std::min( std::abs(cytB_flux), std::abs(E_cytB_near * V_voxel) );
		
	}
	else if ( cytB_flux > 0.0 ) { // Secretion
		adjusted_cytB_flux = std::min(cytB_flux, std::abs(I_cytB * V_cell ) );
	}

	pCell->phenotype.secretion.net_export_rates[cytB_index] = adjusted_cytB_flux;
	// std::cout << "adjusted cytB flux" << adjusted_cytB_flux << std::endl;


	/* 
	
	STEP 2: PA secretion of cytB

	*/

	double I_ATP = parameters.doubles("Initial_I_ATP");
	double I_ADP = parameters.doubles("Initial_I_ADP");

	ATP_hydrolisis_component = diffusion_dt * ATP_hydrolisis_rate * Rcb_cytB;
	if( ATP_hydrolisis_component > Rcb_cytB * diffusion_dt ){ ATP_hydrolisis_component = Rcb_cytB  * diffusion_dt; }
	
	I_ATP -= ATP_hydrolisis_component;
	I_ADP += ATP_hydrolisis_component;
	// I_sK -= movement_component;
	if(I_ATP < 0.0){ I_ATP = 0.0; }
	// if(I_sK < 0.0){ I_sK = 0.0; }

	// step 6: ADP recycling
	ADP_recycling_component = diffusion_dt * ATP_regeneration_rate * I_ADP;
	if (ADP_recycling_component > I_ADP * diffusion_dt){ ADP_recycling_component = I_ADP   * diffusion_dt; }

	I_ADP -= ADP_recycling_component;
	if(I_ADP < 0.0){ I_ADP = 0.0; }
	I_ATP += ADP_recycling_component;


	double cytB_v_max = cytB_movement_rate * (Rc_cytB + Rcb_cytB); 
	double cytB_Km = (cytB_unbinding_rate + ATP_hydrolisis_component) / cytB_binding_rate; // k_n1 + k_2 / k_1

	double cytB_monod = cytB_v_max * (I_cytB / (cytB_Km + I_cytB));
	cytB_monod *= diffusion_dt; 


	double cytB_pump_flux = (cytB_v_max * I_cytB) / ( cytB_Km + I_cytB); // Option A
	cytB_pump_flux *= V_cell  ; // Convert mM/min to amol/min

	pCell->phenotype.secretion.net_export_rates[cytB_index] = adjusted_cytB_flux + cytB_pump_flux;

	

	// // STEP 3: FDC Of sY -- Using just the simple model

	double inhibition_fraction = - (adjusted_cytB_flux + cytB_pump_flux) * diffusion_dt; // Net molar amount
	inhibition_fraction /= V_cell;

	// std::cout << inhibition_fraction << std::endl;

	Rc_sY -= inhibition_fraction * 10;
	Rcb_sY += inhibition_fraction * 10;
	if(Rc_sY < 0.0){ Rc_sY = 0.0;}
	// if(Rcb_sY > total_Rc_sY){Rcb_sY}


	double sY_vf_max = sY_movement_rate * (Rc_sY + Rcb_sY);
	double sY_vb_max = sY_internal_binding_rate * (Rc_sY + Rcb_sY);
	double sY_Km_1 = (sY_internal_binding_rate + sY_movement_rate) / sY_binding_rate;
	double sY_Km_2 = (sY_internal_binding_rate + sY_movement_rate) / sY_internal_binding_rate;

	double sY_velocity = -((sY_vf_max * E_sY_near) / sY_Km_1 - (sY_vb_max * I_sY) / sY_Km_2) / (1 + (I_sY / sY_Km_2) + (E_sY_near / sY_Km_1));
	sY_velocity *= V_cell;

	// if(sY_velocity < 0.0){sY_velocity = 0.0;}

	// std::cout << sY_velocity << std::endl;
	
	pCell->phenotype.secretion.net_export_rates[sY_index] = sY_velocity;


	double sY_monod = sY_vf_max * (I_sY / (sY_Km_1 + I_sY));
	sY_monod *= diffusion_dt; 


	/*

		4. Replace hard-coded initial values at each time-step with the new values

	*/

	// Substrate concentrations
	
	pCell->custom_data.variables[E_sY_near_index].value = E_sY_near; // Mean of near_density_vector
	pCell->custom_data.variables[E_sY_index].value = E_sY_mean;		 // Mean of all voxels
	pCell->custom_data.variables[I_sY_index].value = I_sY;			 // Density

	pCell->custom_data.variables[E_cytB_near_index].value = E_cytB_near; // Mean of near_density_vector
	pCell->custom_data.variables[E_cytB_index].value = E_cytB_mean;		 // Mean of all voxels
	pCell->custom_data.variables[I_cytB_index].value = I_cytB;			 // Density


	// Receptor concentrations

	// pCell->custom_data.variables[total_Rc_index].value = total_Rc;
	// pCell->custom_data.variables[total_Rcb_index].value = Rcb_sY + Rcb_cytB;
	// pCell->custom_data.variables[Rcb_cytB_index].value = Rcb_cytB;
	// pCell->custom_data.variables[Rcb_sY_index].value = Rcb_sY;

	// Kinetics
	pCell->custom_data.variables[DC_sY_index].value = DC_sY; 
	pCell->custom_data.variables[D_sY_index].value = D_sY;	 
	pCell->custom_data.variables[D_cytB_index].value = D_cytB;	 

	pCell->custom_data.variables[sY_binding_component_index].value = sY_binding_component;
	pCell->custom_data.variables[sY_internal_binding_component_index].value = sY_internal_binding_component;
	pCell->custom_data.variables[sY_internal_unbinding_component_index].value = sY_movement_component;
	pCell->custom_data.variables[sY_unbinding_component_index].value = sY_unbinding_component;

	pCell->custom_data.variables[cytB_binding_component_index].value = cytB_binding_component;
	pCell->custom_data.variables[cytB_internal_binding_component_index].value = cytB_internal_binding_component;
	pCell->custom_data.variables[cytB_movement_component_index].value = cytB_movement_component;
	pCell->custom_data.variables[cytB_unbinding_component_index].value = cytB_unbinding_component;
	
	// Net total amounts
	pCell->custom_data.variables[total_E_sY_index].value = E_sY_sum;			 
	pCell->custom_data.variables[total_I_sY_index].value = I_sY_sum;			  
	pCell->custom_data.variables[total_sY_index].value = total_E_sY + total_I_sY; 
	pCell->custom_data.variables[total_E_cytB_index].value = E_cytB_sum;			 
	pCell->custom_data.variables[total_I_cytB_index].value = I_cytB_sum;			  
	pCell->custom_data.variables[total_cytB_index].value = total_E_cytB + total_I_cytB; 


	/*

		5. Connecting the model with different aspects of the phenotype

	*/

	// pCell->phenotype.secretion.net_export_rates[sY_index] = sY_velocity;
	// pCell->phenotype.secretion.net_export_rates[cytB_index] = cytB_velocity;

	static int necrosis_idx = pCell->phenotype.death.find_death_model_index("Necrosis");
	static int apoptosis_idx = pCell->phenotype.death.find_death_model_index("Apoptosis");

	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 

	pCell->phenotype.cycle.data.transition_rate(cycle_start_index, cycle_end_index) = I_sY * 1e-02;


	
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
	static int sY_index = microenvironment.find_density_index("sY");
	double E_sY_near = pCell->nearest_density_vector()[sY_index];
	double I_sY = pCell->phenotype.molecular.internalized_total_substrates[sY_index];

	double delta = I_sY - E_sY_near; // Difference between inside and outside


	std::vector<std::string> output = false_cell_coloring_live_dead(pCell);

	if (pCell->type == 0) // "default" cell type
	{
		int color = I_sY;
		char szTempString[182];
		sprintf(szTempString, "rgb(100, %u, 100)", color);
		output[0].assign(szTempString); // Cytoplasm

		sprintf(szTempString, "rgb(100, %u, 100)", color);
		output[2] = "rgb(255,0,0)";		// nucleus

	}

	return output;
}

void color_node(Cell *pCell)
{
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}

void drug_resistance_model_main(double dt)
{
#pragma omp parallel for
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell *pC = (*all_cells)[i];
		if (pC->phenotype.death.dead == false)
		{
			drug_resistance_model(pC, pC->phenotype, dt);
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