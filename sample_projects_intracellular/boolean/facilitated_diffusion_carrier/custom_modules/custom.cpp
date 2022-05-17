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
	cell_defaults.functions.update_phenotype = facilitated_carrier_model;
	cell_defaults.functions.custom_cell_rule = NULL;

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;

	// add here custom variables (same as adding them on the XML file) - for analysis with pyMCDS
	// ??? Maybe better to add them within the XML file, for easier manipulation???

	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0); // for paraview visualization

	// Naming convention -  I: Internal; E: External
	// !!! Keep the same substrate names as defined in the microenvironment

	static int sY_index = microenvironment.find_density_index("substrate_Y");
	Cell *pCell;

	// Generic variables:
	cell_defaults.custom_data.add_variable("I_sY", "amol / um^3", parameters.doubles("Initial_I_sY"));
	cell_defaults.custom_data.add_variable("E_sY_near", "amol / um^3", default_microenvironment_options.initial_condition_vector[sY_index]); // Mean nearest substrate conc.
	cell_defaults.custom_data.add_variable("E_sY", "amol/um^3", default_microenvironment_options.initial_condition_vector[sY_index]);		 // Mean external substrate conc.

	cell_defaults.custom_data.add_variable("Rc_sY", "amol / um^3", parameters.doubles("Initial_Rc_sY"));											  // Receptor density
	cell_defaults.custom_data.add_variable("Rcb_sY", "amol / um^3", parameters.doubles("Initial_Rcb_sY"));											  // Bound receptor density
	cell_defaults.custom_data.add_variable("total_Rc_sY", "amol / um^3", parameters.doubles("Initial_Rcb_sY") + parameters.doubles("Initial_Rc_sY")); // Bound receptor density

	cell_defaults.custom_data.add_variable("total_E_sY", "amol", 0.0); // total external sY
	cell_defaults.custom_data.add_variable("total_I_sY", "amol", 0.0); // total internal sY
	cell_defaults.custom_data.add_variable("total_sY", "amol", 0.0);   // total internal sY
	cell_defaults.custom_data.add_variable("sY_flux", "amol/min", 0.0);
	cell_defaults.custom_data.add_variable("sY_flux_explicit", "amol/min", 0.0);
	cell_defaults.custom_data.add_variable("adjusted_sY_flux", "amol/min", 0.0);

	cell_defaults.custom_data.add_variable("Initial_E_sY", "amol/um^3", default_microenvironment_options.initial_condition_vector[sY_index]); // Initial external substrate_Y
	cell_defaults.custom_data.add_variable("Initial_I_sY", "amol/um^3", parameters.doubles("Initial_I_sY"));								  // Initial internal substrate_Y

	cell_defaults.custom_data.add_variable("DC_sY", "um^3/min", microenvironment.diffusion_coefficients[sY_index]); // Diffusion coefficient
	cell_defaults.custom_data.add_variable("D_sY", "mmol/um^3", 0.0);

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

	static int sY_index = microenvironment.find_density_index("sY");

	double Initial_Rc_sY = parameters.doubles("Initial_Rc_sY");
	double Initial_Rcb_sY = parameters.doubles("Initial_Rcb_sY");

	double Initial_I_sY = parameters.doubles("Initial_I_sY");
	int cells_number = parameters.ints("cells_number");
	int just_one_cell = parameters.ints("one_cell");
	int tumor_radius = parameters.ints("tumor_radius");

	double center_x;
	double center_y;
	double center_z;

	// double Rc_sY_index = pCell->custom_data.find_variable_index("Rcb_sY"); // [E] Receptor alone (mM)

	// for(int k=0; k<parameters.ints("cells_number"); k++)
	// {

	// 	std::cout << "Placing cells.. " << k << std::endl;

	// 	std::vector<double> position = {0.0, 0.0, 0.0};
	// 	double r = NormalRandom(0, 1) * tumor_radius;
	// 	double theta = UniformRandom() * M_PI * 2;

	// 	position[0] = center_x + r*std::cos(theta);
	// 	position[1] = center_y + r*std::sin(theta);
	// 	position[2] = center_z;

	// 	pCell->assign_position( position );
      
	// 	pCell->phenotype.molecular.internalized_total_substrates[sY_index] = parameters.doubles("Initial_I_sY");
	
	// }




	if (just_one_cell == 1)
	{
		#pragma omp parallel for
		for (int i=0; i < parameters.ints("cells_number"); i++) 
		{
			Cell* pC;
			pC = create_cell(get_cell_definition("default")); 

			// Distribute cells randomly on Unit Circle
			std::vector<double> position = UniformOnUnitCircle();
			position *= 50.0; position += {0,0,0}; // circle: radius of 50.0, centered at 0,0,0
			pC -> assign_position( position );
			pC -> assign_position( {0,0,0} );

			// Add randomly within tumor radius
			// std::vector<double> position = {0.0, 0.0, 0.0};

			// double r = NormalRandom(0, 1) * tumor_radius;
			
			// double theta = UniformRandom() * M_PI * 2;

			// position[0] = center_x + r*std::cos(theta);
			// position[1] = center_y + r*std::sin(theta);
			// position[2] = center_z;

			// pC->assign_position( position );

			// Add initial density of substrate -- Input is a density, convert to a net value
			pC -> phenotype.molecular.internalized_total_substrates[sY_index] = Initial_I_sY * pC->phenotype.volume.total;
			
		}

	}
	
	if(just_one_cell == 0)
	{
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

				if (fabs(y) > 0.01)
				{
					pCell = create_cell(); // tumor cell
					pCell->assign_position(x, -y, 0.0);
					pCell->phenotype.molecular.internalized_total_substrates[sY_index] = Initial_I_sY * pCell->phenotype.volume.total;
				}

				if (fabs(x) > 0.01)
				{
					pCell = create_cell(); // tumor cell
					pCell->assign_position(-x, y, 0.0);
					pCell->phenotype.molecular.internalized_total_substrates[sY_index] = Initial_I_sY * pCell->phenotype.volume.total;

					if (fabs(y) > 0.01)
					{
						pCell = create_cell(); // tumor cell
						pCell->assign_position(-x, -y, 0.0);
						pCell->phenotype.molecular.internalized_total_substrates[sY_index] = Initial_I_sY * pCell->phenotype.volume.total;
					
					}
				}
				x += cell_spacing;
			}

			y += cell_spacing * sqrt(3.0) / 2.0;
			n++;
		}

		// load cells from your CSV file (if enabled)
		// load_cells_from_pugixml();


	}

	
	
	return;

}


/*
	Facilitated diffusion model

	ODE model for simple diffusion based on a set of differential equations.

	Parameters:
	k_1: Binding rate
	k_n1: Recycling rate
	k_2: Endocytosis rate
	k_n2: Internal binding rate

	Species:
	[S]_out: External substrate
	[S]_in: Internal substrate
	[E]: Receptor
	[ES]: Bound receptor

	ODE-system: 4 equations

*/

void facilitated_carrier_model(Cell *pCell, Phenotype &phenotype, double dt)
{

	static int sY_index = microenvironment.find_density_index("sY"); // microenv substrate index

	static int E_sY_near_index = pCell->custom_data.find_variable_index("E_sY_near");
	static int I_sY_index = pCell->custom_data.find_variable_index("I_sY");
	static int E_sY_index = pCell->custom_data.find_variable_index("E_sY");

	double Rc_sY_index = pCell->custom_data.find_variable_index("Rc_sY"); // [E] Receptor alone (mM)
	double Rcb_sY_index = pCell->custom_data.find_variable_index("Rcb_sY");
	double total_Rc_sY_index = pCell->custom_data.find_variable_index("total_Rc_sY");

	static int DC_sY_index = pCell->custom_data.find_variable_index("DC_sY"); // Difussion coefficient
	static int D_sY_index = pCell->custom_data.find_variable_index("D_sY");	  // Difussion coefficient
	static int total_E_sY_index = pCell->custom_data.find_variable_index("total_E_sY");
	static int total_I_sY_index = pCell->custom_data.find_variable_index("total_I_sY");
	static int total_sY_index = pCell->custom_data.find_variable_index("total_sY");

	static int sY_flux_index = pCell->custom_data.find_variable_index("sY_flux");
	static int sY_flux_explicit_index = pCell->custom_data.find_variable_index("sY_flux_explicit");
	static int adjusted_sY_flux_index = pCell->custom_data.find_variable_index("adjusted_sY_flux");
	static int num_cells_idx = pCell->custom_data.find_variable_index("num_cells"); 

	static int time_index = pCell->custom_data.find_variable_index("time");
	pCell->custom_data.variables[time_index].value = PhysiCell_globals.current_time;

	// Constants

	static double DC_sY = microenvironment.diffusion_coefficients[sY_index];
	static double V_cell = pCell->phenotype.volume.total;
	static double R_cell = pCell->phenotype.geometry.radius;
	static double A_cell = 4 * M_PI * std::pow(R_cell, 2);

	double V_voxel = (microenvironment.mesh.dV);
	double total_num_voxels = microenvironment.number_of_voxels();

	// Kinetic parameters
	static double sY_binding_rate = parameters.doubles("k_1");
	static double sY_recycling_rate = parameters.doubles("k_n1");
	static double sY_endocytosis_rate = parameters.doubles("k_2");
	static double sY_internal_binding_rate = parameters.doubles("k_n2");

	// Receptor concentrations
	double Rc_sY = pCell->custom_data.variables[Rc_sY_index].value;	  // [E] Receptor alone (mM)
	double Rcb_sY = pCell->custom_data.variables[Rcb_sY_index].value; // [ES] Bound receptor (mM)

	// Substrate concentrations
	double E_sY_near = pCell->nearest_density_vector()[sY_index];							   // A density (mM)
	double I_sY = pCell->phenotype.molecular.internalized_total_substrates[sY_index] / V_cell; // Convert to density (mM)

	static double Initial_E_sY = default_microenvironment_options.initial_condition_vector[sY_index];
	static double Initial_I_sY = parameters.doubles("Initial_I_sY"); 

	pCell->custom_data.variables[num_cells_idx].value = (*all_cells).size(); // Total number of cells


	// Obtaining the total net internal and external amounts through iterating over the cells

	double total_E_sY;
	double total_I_sY;

	std::vector<double> I_sY_vector(microenvironment.number_of_voxels());

	// Counting net internal amount
	#pragma omp parallel for
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell *pC = (*all_cells)[i];
		I_sY_vector[i] = pC->phenotype.molecular.internalized_total_substrates[sY_index];
		total_I_sY += pC->phenotype.molecular.internalized_total_substrates[sY_index];
	}

	double I_sY_sum = std::accumulate(I_sY_vector.begin(), I_sY_vector.end(), 0.0);

	// I_sY_sum and total_I_sY provide the same exact result. They reflect the actual total amount of
	// internal net sY

	std::vector<double> E_sY_vector(microenvironment.number_of_voxels());
	
	#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		E_sY_vector[n] = (microenvironment.density_vector(n)[sY_index]);
		total_E_sY += (microenvironment.density_vector(n)[sY_index]) * V_voxel; // A net amount of sY
	}

	double E_sY_sum = std::accumulate(E_sY_vector.begin(), E_sY_vector.end(), 0.0);
	double E_sY_mean = E_sY_sum / total_num_voxels; // NOT CORRECT, the accumulate just adds up the number of voxels

	double D_sY = I_sY - E_sY_near;

	// std::cout << "I_sY: " << I_sY << std::endl;
	// std::cout << "E_sY_near: " << E_sY_near << std::endl;
	// std::cout << "D_sY: " << D_sY << std::endl;

	// replace hard-coded initial values at each time-step with the new values

	// Substrate concentrations
	
	pCell->custom_data.variables[E_sY_near_index].value = E_sY_near; // Mean of near_density_vector
	pCell->custom_data.variables[E_sY_index].value = E_sY_mean;		 // Mean of all voxels
	pCell->custom_data.variables[I_sY_index].value = I_sY;			 // Density

	pCell->custom_data.variables[DC_sY_index].value = DC_sY; // Substrate diffusion coefficient
	pCell->custom_data.variables[D_sY_index].value = D_sY;	 // Substrate diffusion coefficient


	// Receptor concentrations

	pCell->custom_data.variables[total_Rc_sY_index].value = Rcb_sY + Rc_sY;

	/*

	The modification is in the net_export_rate, but it depends on the direction.
	We assume 1:1 stoichiometry between receptor and substrate, thus the reaction will be of order 1.

	*/

	// Separated components

	double binding_component;
	double recycling_component;
	double endocytosis_component;
	double internal_binding_component;

	binding_component = diffusion_dt * sY_binding_rate * Rc_sY * E_sY_near;
	if (binding_component > Rc_sY || binding_component > E_sY_near )
	{ binding_component = std::min(Rc_sY, E_sY_near); } // Get the limiting step 
	
	Rcb_sY += binding_component;
	Rc_sY -= binding_component;
	// E_sY_near -= binding_component;
	if(Rc_sY < 0.0) { Rc_sY = 0.0; }
	if(E_sY_near < 0.0) { E_sY_near = 0.0; }

	recycling_component = diffusion_dt * sY_recycling_rate * Rcb_sY;
	if (recycling_component > Rcb_sY){ recycling_component = Rcb_sY; } 
	Rc_sY += recycling_component;
	Rcb_sY -= recycling_component;
	// E_sY_near += recycling_component;
	if(Rcb_sY < 0.0) { Rcb_sY = 0.0; }

	endocytosis_component = diffusion_dt * sY_endocytosis_rate * Rcb_sY;
	if (endocytosis_component > Rcb_sY){ endocytosis_component = Rcb_sY; } 
	I_sY += endocytosis_component;
	Rcb_sY -= endocytosis_component;
	Rc_sY += endocytosis_component;
	if(Rcb_sY < 0.0) { Rcb_sY = 0.0; }
	// if(I_sY > E_sY_near) { endocytosis_component = 0.0; } // hard-code gradient

	internal_binding_component = diffusion_dt * sY_internal_binding_rate * Rc_sY * I_sY;
	if (internal_binding_component > Rcb_sY || internal_binding_component > I_sY )
	{ internal_binding_component = std::min(Rcb_sY, I_sY); } // Get the limiting step 
	
	I_sY -= internal_binding_component;
	Rcb_sY += internal_binding_component;
	Rc_sY -= internal_binding_component;

	if(Rc_sY < 0.0) { Rc_sY = 0.0; }
	if(I_sY < 0.0) { I_sY = 0.0; }


	// double vf_max = sY_endocytosis_rate * (Rc_sY + Rcb_sY); // both are constant given that total receptor is constant
	// double vb_max = sY_recycling_rate * (Rc_sY + Rcb_sY);
	// double Km_1 = (sY_recycling_rate + sY_endocytosis_rate) / sY_binding_rate;
	// double Km_2 = (sY_recycling_rate + sY_endocytosis_rate) / sY_internal_binding_rate;

	double vf_max = parameters.doubles("v_max_f");
	double vb_max = parameters.doubles("v_max_r");
	double Km_1 = parameters.doubles("K_M_1");
	double Km_2 = parameters.doubles("K_M_2");


	double sY_monod = vf_max * (I_sY / (Km_1 + I_sY)); // Using the "forward" (entry) kinetic values
	sY_monod *= diffusion_dt;


	double v_formation = -((vf_max * E_sY_near) / Km_1 - (vb_max * I_sY) / Km_2) / (1 + (I_sY / Km_2) + (E_sY_near / Km_1));
	v_formation *= V_cell;
	// std::cout << "v_formation: " << v_formation << std::endl;

	double v_formation_explicit = - endocytosis_component + internal_binding_component;
	v_formation_explicit *= V_cell;
	// std::cout  << "v_formation_explicit: " << v_formation_explicit << std::endl;

	pCell->custom_data.variables[sY_flux_index].value = v_formation;
	pCell->custom_data.variables[sY_flux_explicit_index].value = v_formation_explicit / diffusion_dt;
	pCell->custom_data.variables[Rc_sY_index].value = Rc_sY;
	pCell->custom_data.variables[Rcb_sY_index].value = Rcb_sY;

	if(vf_max == 0 && vb_max == 0 && Km_1 == 0 && Km_2 == 0)
	{
		// If the explicit system is employed, account for the dt multiplications made when obtaining the deltas.
		pCell->phenotype.secretion.net_export_rates[sY_index] = v_formation_explicit / diffusion_dt;
	}
	else
	{
		// If the velocity equation is employed
		pCell->phenotype.secretion.net_export_rates[sY_index] = v_formation;
	}


	// pCell->phenotype.secretion.net_export_rates[sY_index] = v_formation;
	
	// Net total amounts
	pCell->custom_data.variables[total_E_sY_index].value = E_sY_sum;			 
	pCell->custom_data.variables[total_I_sY_index].value = I_sY_sum;			  
	pCell->custom_data.variables[total_sY_index].value = total_E_sY + total_I_sY; 

	// Do not update phenotype if dead
	if (pCell->phenotype.death.dead == true)
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	/*

		Connecting the model with different aspects of the phenotype

	*/

	static int necrosis_idx = pCell->phenotype.death.find_death_model_index("Necrosis");
	static int apoptosis_idx = pCell->phenotype.death.find_death_model_index("Apoptosis");

	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 

	// pCell->phenotype.cycle.data.transition_rate(0,0) = sY_monod * V_cell;




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

void facilitated_carrier_main(double dt)
{
#pragma omp parallel for
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell *pC = (*all_cells)[i];
		if (pC->phenotype.death.dead == false)
		{
			facilitated_carrier_model(pC, pC->phenotype, dt);
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

	// v_formation = v_formation * V_cell; // Is this necessary? Copasi does it
	// double v_formation = - binding_component + recycling_component - endocytosis_component + internal_binding_component;
	// double v_formation = - endocytosis_component;