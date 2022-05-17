
///////////// USEFUL CHUNKS OF CODE


/*


============= FROM setup_tissue, different ways to initially place the cells

	FIRST OPTION: USE THE SAME AS THE 

	double Xmin = microenvironment.mesh.bounding_box[0]; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = microenvironment.mesh.bounding_box[2]; 

	double Xmax = microenvironment.mesh.bounding_box[3]; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = microenvironment.mesh.bounding_box[5]; 
	
	if( default_microenvironment_options.simulate_2D == true )
	{
		Zmin = 0.0; 
		Zmax = 0.0; 
	}
	
	double Xrange = Xmax - Xmin; 
	double Yrange = Ymax - Ymin; 
	double Zrange = Zmax - Zmin; 

	for( int k=0; k < cell_definitions_by_index.size() ; k++ )
	{
		Cell_Definition* pCD = cell_definitions_by_index[k]; 
		std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl; 
		for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
		{
			std::vector<double> position = {0,0,0}; 
			position[0] = Xmin + UniformRandom()*Xrange; 
			position[1] = Ymin + UniformRandom()*Yrange; 
			position[2] = Zmin + UniformRandom()*Zrange; 
			
			pC = create_cell( *pCD ); 
			pC->assign_position( position );
		}
	}
	std::cout << std::endl; 


	ANOTHER OPTION: Just placing the cells, and that's it
	
	// We have four different cells populations
	// All start with A active and C inactive, B is random
	// We print the value of C
	
	for (int i=0; i < 11; i+= 11)
		for (int j=0; j < 11; j+= 11){
			
			// bottom left : default
			// the formula for C is A&B. Meaning that C will only activate for half the cells
			pC = create_cell(get_cell_definition("default")); 
			pC->assign_position(-i-10, -j-10, 0.0 );

						
			// bottom middle : other
			// the formula for C is A|B. C will activate in all cells
			pC = create_cell(get_cell_definition("other")); 
			pC->assign_position(i+10, -j-10, 0.0 );

			// top left : another
			// Here we mutate the C node at zero, so it will stay there
			pC = create_cell(get_cell_definition("another")); 
			pC->assign_position(-i-10, j+10, 0.0 );
			
			// top middle : yet_another
			// Here we change the default value for the rates, acelerating the activation of C
			pC = create_cell(get_cell_definition("yet_another")); 
			pC->assign_position(i+10, j+10, 0.0 );
			
			// top right : yet_yet_another
			// Here we acelerate the activation of C by changing the scaling value
			pC = create_cell(get_cell_definition("yet_yet_another")); 
			pC->assign_position(i+110, j+10, 0.0 );
			
			// bottom right : last_one
			// Here we start with $time_scale = 0, then at the middle of the simulation we set it to 0.1
			pC = create_cell(get_cell_definition("last_one")); 
			pC->assign_position(i+110, -j-10, 0.0 );
		}





====================================================================================

 //////////////////////  COLORING FUNCTIONS

 
		if ( pCell -> phenotype.death.dead == false) 
	{
		output[1] = "rgb(201,51,242)";
		output[0] = "rgb(242,161,29)"; // 
		output[2] = "rgb(255,0,0)";

		return output;
	}

	Previous coloring function, based on the node to visualize
	
	if ( !pCell->phenotype.intracellular->get_boolean_variable_value( parameters.strings("node_to_visualize") ) )
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
		
	}
	else{
		output[0] = "rgb(0, 255,0)";
		output[2] = "rgb(0, 125,0)";
	}



=====================================================================================================

------------------------------ SPECIFIC PHENOTYPE FUNCTIONS

--- Example from the User Guide

void custom_o2_phenotype_rule( Cell* pCell, Phenotype& phenotype, double dt )
{
    // don’t bother if you’re dead
    if( pCell->phenotype.death.dead == true )
    { return; }

    // first, call the standard function
    update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);

    // next, let’s evaluate the oxygen
    static int o2_index = microenvironment.find_density_index("oxygen");
    double o2 = pCell->nearest_density_vector()[o2_index];

    if( o2 > pCell->parameters.o2_hypoxic_response )
    { return; }

    // interpolation variable
    double theta = ( pCell->parameters.o2_hypoxic_response - o2 )/
    (pCell->parameters.o2_hypoxic_response - pCell->parameters.o2_hypoxic_saturation);
    
    if( theta > 1.0 )
    { theta = 1.0; }
    // increase the speed of motiltiy up to a max of 1.5 micron/min
    phenotype.motility.is_motile = true;
    phenotype.motility.migration_speed = 1.5;
    phenotype.motility.migration_speed *= theta;
    phenotype.mechanics.cell_cell_adhesion_strength = (1.0-theta);
    phenotype.mechanics.cell_cell_adhesion_strength *=
    pCell->parameters.pReference_live_phenotype->mechanics.cell_cell_adhesion_strength;
    return;
}

--- This one is from the heterogeneity_sample

// custom cell phenotype function to scale immunostimulatory factor with hypoxia 
void tumor_cell_phenotype_with_oncoprotein( Cell* pCell, Phenotype& phenotype, double dt )
{
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	
	// if cell is dead, don't bother with future phenotype changes. 
	if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL; 		
		return; 
	}

	// multiply proliferation rate by the oncoprotein 
	
	static int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	static int oncoprotein_i = pCell->custom_data.find_variable_index( "oncoprotein" ); 

	phenotype.cycle.data.transition_rate( cycle_start_index ,cycle_end_index ) *= pCell->custom_data[oncoprotein_i] ; 
	
	return; 
}


----- Modify motility based on O2, from a PhysiCell UG template

	/*

		Modify motility based on O2

		

	if( o2 > pCell->parameters.o2_hypoxic_response )
    { return; }

    // interpolation variable
    double theta = ( pCell->parameters.o2_hypoxic_response - o2 )/
    (pCell->parameters.o2_hypoxic_response - pCell->parameters.o2_hypoxic_saturation);
    
    if( theta > 1.0 )
    { theta = 1.0; }
	
    // increase the speed of motiltiy up to a max of 1.5 micron/min
    phenotype.motility.is_motile = true;
    phenotype.motility.migration_speed = 1.5;
    phenotype.motility.migration_speed *= theta;
    phenotype.mechanics.cell_cell_adhesion_strength = (1.0-theta);
    phenotype.mechanics.cell_cell_adhesion_strength *=
    pCell->parameters.pReference_live_phenotype->mechanics.cell_cell_adhesion_strength;
	return;

	*/


	
/*

		<!-- 
			
			"other" cell type. Inherits some phenotype features from the default.

		 -->

		<cell_definition name="other" parent_type="default" ID="1">
			<phenotype>
				<cycle code="5" name="live">  
					<!-- using higher than normal significant digits to match divisions in default code -->
					<phase_transition_rates units="1/min"> 
						<rate start_index="0" end_index="0" fixed_duration="false">0.003333</rate>
					</phase_transition_rates>
				</cycle>
				
				<death>
					<model code="100" name="apoptosis"> 
						<death_rate units="1/min">0</death_rate>
						<phase_durations units="min">
							<duration index="0" fixed_duration="true">0</duration>
						</phase_durations>
						<parameters>
							<unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
							<lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
							<cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
							<nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
							<calcification_rate units="1/min">0</calcification_rate>
							<relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
						</parameters>
					</model> 

					<model code="101" name="necrosis">
						<death_rate units="1/min">0.0</death_rate>
						<phase_durations units="min">
							<duration index="0" fixed_duration="true">0</duration>
							<duration index="1" fixed_duration="true">0</duration>
						</phase_durations>
						
						<parameters>
							<unlysed_fluid_change_rate units="1/min">0.05</unlysed_fluid_change_rate>
							<lysed_fluid_change_rate units="1/min">0</lysed_fluid_change_rate>
							<cytoplasmic_biomass_change_rate units="1/min">1.66667e-02</cytoplasmic_biomass_change_rate>
							<nuclear_biomass_change_rate units="1/min">5.83333e-03</nuclear_biomass_change_rate>
							<calcification_rate units="1/min">0</calcification_rate>
							<relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
						</parameters>
					</model> 
				</death>						
				
				<volume>
					<total units="micron^3">2494</total>
					<fluid_fraction units="dimensionless">0.75</fluid_fraction>
					<nuclear units="micron^3">540</nuclear>
					
					<fluid_change_rate units="1/min">0.05</fluid_change_rate>
					<cytoplasmic_biomass_change_rate units="1/min">0.0045</cytoplasmic_biomass_change_rate>
					<nuclear_biomass_change_rate units="1/min">0.0055</nuclear_biomass_change_rate>
					
					<calcified_fraction units="dimensionless">0</calcified_fraction>
					<calcification_rate units="1/min">0</calcification_rate>
					
					<relative_rupture_volume units="dimensionless">2.0</relative_rupture_volume>
				</volume>

				<secretion>
					<substrate name="oxygen">
						<secretion_rate units="1/min">10e-08</secretion_rate>
						<secretion_target units="substrate density">1</secretion_target>
						<uptake_rate units="1/min">10e-07</uptake_rate>
						<net_export_rate units="total substrate/min">0</net_export_rate> 
					</substrate>

					<substrate name="prot_A">
						<secretion_rate units="1/min">0</secretion_rate>
						<secretion_target units="substrate density">1</secretion_target>
						<uptake_rate units="1/min">0</uptake_rate>
						<net_export_rate units="total substrate/min">0</net_export_rate> 
					</substrate> 					
				</secretion>

				<motility> 
					<speed units="micron/min">1</speed>
					<persistence_time units="min">1</persistence_time> 
					<migration_bias units="dimensionless">0</migration_bias>
					
					<options>
						<enabled>true</enabled>
						<use_2D>true</use_2D>
						<chemotaxis>
							<enabled>true</enabled>
							<substrate>oxygen</substrate>
							<direction>1</direction>
						</chemotaxis>
					</options>
				</motility>			

				<!-- 

				Leave the inherited one from the DEFAULT cell type, in order to avoid differences
				seen due to the internal nodes 


				<intracellular type="maboss">
					<bnd_filename>./config/model_0.bnd</bnd_filename>
				</intracellular>
				-->
				
			
			</phenotype>
			<custom_data/>




	-----------------------------------------------------------------------
	These are the other phenotypes that were described in the XML file
	-----------------------------------------------------------------------

		<!-- I really don't need these cell type definitions, but could be useful to try other stuff -->
		</cell_definition>
		<cell_definition name="another" parent_type="default" ID="2">
			<phenotype>
				<intracellular type="maboss">
					<mutations>
						<mutation node="C">0.0</mutation>
					</mutations>
				</intracellular>			
			</phenotype>
			<custom_data/>
		</cell_definition>
		<cell_definition name="yet_another" parent_type="default" ID="3">
			<phenotype>
				<intracellular type="maboss">
					<parameters>
						<parameter name="$time_scale">0.2</parameter>
					</parameters>
				</intracellular>
			</phenotype>
			<custom_data/>
		</cell_definition>
		<cell_definition name="yet_yet_another" parent_type="default" ID="4">
			<phenotype>
				<intracellular type="maboss">
					<scaling>0.25</scaling>
				</intracellular>
			</phenotype>
			<custom_data/>
		</cell_definition>
		<cell_definition name="last_one" parent_type="default" ID="5">
			<phenotype>
				<intracellular type="maboss">
					<parameters>
						<parameter name="$time_scale">0.0</parameter>
					</parameters>
				</intracellular>
			</phenotype>
			<custom_data/>
		</cell_definition>
	


*/



/////////// UNIFORM DISTRIBUTION WITHIN CELLS in setup_tissue(

	int n = 0;
	while (y < tumor_radius)
	{
		// static int Rc_sY_index = pCell->custom_data.find_variable_index("Rc_sY"); // [E] Receptor alone (mM)

		x = 0.0;
		if (n % 2 == 1)
		{
			x = 0.5 * cell_spacing;
		}
		x_outer = sqrt(tumor_radius * tumor_radius - y * y);

		while (x < x_outer)
		{
			pCell = create_cell(pCH); // healthy cell
			pCell = create_cell(pCC); // healthy cell
			pCell->assign_position(x, y, 0.0);
			pCell->custom_data[4] = NormalRandom(p_mean, p_sd);
			if (pCell->custom_data[4] < p_min)
			{
				pCell->custom_data[4] = p_min;
			}
			if (pCell->custom_data[4] > p_max)
			{
				pCell->custom_data[4] = p_max;
			}
			std::cout << "pCell->custom_data[3] = " << pCell->custom_data[4] << std::endl;

			if (fabs(y) > 0.01)
			{
				pCell = create_cell(); // tumor cell
				pCell->assign_position(x, -y, 0.0);
				pCell->custom_data[4] = NormalRandom(p_mean, p_sd);
				if (pCell->custom_data[4] < p_min)
				{
					pCell->custom_data[4] = p_min;
				}
				if (pCell->custom_data[4] > p_max)
				{
					pCell->custom_data[4] = p_max;
				}
				std::cout << "pCell->custom_data[3] = " << pCell->custom_data[4] << std::endl;
			}

			if (fabs(x) > 0.01)
			{
				pCell = create_cell(); // tumor cell
				pCell->assign_position(-x, y, 0.0);
				pCell->custom_data[4] = NormalRandom(p_mean, p_sd);
				if (pCell->custom_data[4] < p_min)
				{
					pCell->custom_data[4] = p_min;
				}
				if (pCell->custom_data[4] > p_max)
				{
					pCell->custom_data[4] = p_max;
				}
				std::cout << "pCell->custom_data[3] = " << pCell->custom_data[4] << std::endl;

				if (fabs(y) > 0.01)
				{
					pCell = create_cell(); // tumor cell
					pCell->assign_position(-x, -y, 0.0);
					pCell->custom_data[4] = NormalRandom(p_mean, p_sd);
					if (pCell->custom_data[4] < p_min)
					{
						pCell->custom_data[4] = p_min;
					}
					if (pCell->custom_data[4] > p_max)
					{
						pCell->custom_data[4] = p_max;
					}
					std::cout << "pCell->custom_data[3] = " << pCell->custom_data[4] << std::endl;
				}
			}
			x += cell_spacing;
		}

		y += cell_spacing * sqrt(3.0) / 2.0;
		n++;
	}

	double sum = 0.0;
	double min = 9e9;
	double max = -9e9;
	for (int i = 0; i < all_cells->size(); i++)
	{
		double r = (*all_cells)[i]->custom_data[4];
		sum += r;
		if (r < min)
		{
			min = r;
		}
		if (r > max)
		{
			max = r;
		}
	}
	double mean = sum / (all_cells->size() + 1e-15);
	// compute standard deviation

	sum = 0.0;
	for (int i = 0; i < all_cells->size(); i++)
	{
		sum += ((*all_cells)[i]->custom_data[4] - mean) * ((*all_cells)[i]->custom_data[4] - mean);
	}
	double standard_deviation = sqrt(sum / (all_cells->size() - 1.0 + 1e-15));

	std::cout << std::endl
			  << "prot_A summary: " << std::endl
			  << "===================" << std::endl;
	std::cout << "mean: " << mean << std::endl;
	std::cout << "standard deviation: " << standard_deviation << std::endl;
	std::cout << "[min max]: [" << min << " " << max << "]" << std::endl
			  << std::endl;

	// load cells from your CSV file (if enabled)
	// load_cells_from_pugixml();








///////////// UNIT SPHERE PLACEMENT OF CELLS

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
} // namespace bdm


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

	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0); // for paraview visualization

	static int sY_index = microenvironment.find_density_index("substrate_Y");

	Cell *pCell;

	// Generic variables:
	cell_defaults.custom_data.add_variable("I_sY", "amol / um^3", parameters.doubles("Initial_I_sY"));
	cell_defaults.custom_data.add_variable("E_sY_near", "amol / um^3", default_microenvironment_options.initial_condition_vector[sY_index]); // Mean nearest substrate conc.
	cell_defaults.custom_data.add_variable("E_sY", "amol/um^3", default_microenvironment_options.initial_condition_vector[sY_index]);		 // Mean external substrate conc.

	cell_defaults.custom_data.add_variable("Rc_sY", "amol/um^3", parameters.doubles("Initial_Rc_sY")); // Initial internal substrate_Y
	cell_defaults.custom_data.add_variable("Rcb_sY", "amol / um^3", parameters.doubles("Initial_Rcb_sY"));
	cell_defaults.custom_data.add_variable("total_Rc_sY", "amol / um^3", parameters.doubles("Initial_Rcb_sY") + parameters.doubles("Initial_Rc_sY")); // Bound receptor density

	cell_defaults.custom_data.add_variable("total_E_sY", "amol", 0.0); // total external sY
	cell_defaults.custom_data.add_variable("total_I_sY", "amol", 0.0); // total internal sY
	cell_defaults.custom_data.add_variable("total_sY", "amol", 0.0);   // total internal sY
	cell_defaults.custom_data.add_variable("sY_flux", "amol/min", 0.0);
	cell_defaults.custom_data.add_variable("adjusted_sY_flux", "amol/min", 0.0);

	cell_defaults.custom_data.add_variable("Initial_E_sY", "amol/um^3", default_microenvironment_options.initial_condition_vector[sY_index]); // Initial external substrate_Y
	cell_defaults.custom_data.add_variable("Initial_I_sY", "amol/um^3", parameters.doubles("Initial_I_sY"));								  // Initial internal substrate_Y

	cell_defaults.custom_data.add_variable("DC_sY", "um^3/min", microenvironment.diffusion_coefficients[sY_index]); // Diffusion coefficient
	cell_defaults.custom_data.add_variable("D_sY", "mmol/um^3", 0.0);
	// cell_defaults.custom_data.add_variable("k_sY", "um^2/min", parameters.doubles("k_sY") );

	cell_defaults.custom_data.add_variable("number_of_voxels", "dimensionless", microenvironment.number_of_voxels());

	/*
	   This parses the cell definitions in the XML config file.
	*/

	initialize_cell_definitions_from_pugixml();

	/*
	   Put any modifications to individual cell definitions here.

	   This is a good place to set custom functions.
	*/

	Cell_Definition* pCH = find_cell_definition("healthy");
	Cell_Definition* pCC = find_cell_definition("cancer_cell");

	// Introduce here the custom cell rules 

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
	static Cell_Definition* pCH = find_cell_definition("healthy");
	static Cell_Definition* pCC = find_cell_definition("cancer_cell");

	double Initial_Rc_sY = parameters.doubles("Initial_Rc_sY");
	double Initial_Rcb_sY = parameters.doubles("Initial_Rcb_sY");

	double Initial_I_sY = parameters.doubles("Initial_I_sY");
	static int sY_index = microenvironment.find_density_index("sY");

	// double Rc_sY_index = pCell->custom_data.find_variable_index("Rcb_sY"); // [E] Receptor alone (mM)

	double cell_radius = cell_defaults.phenotype.geometry.radius;
	double cell_spacing = 0.95 * 2.0 * cell_radius;
	double tumor_radius = parameters.doubles("tumor_radius"); // 250.0;

	// Parameter<double> temp;

	double x = 0.0;
	double x_outer = tumor_radius;
	double y = 0.0;

	double R_mean = parameters.doubles("R_mean");
	double R_sd = parameters.doubles("R_sd");
	double R_min = parameters.doubles("R_min");
	double R_max = parameters.doubles("R_max");


	double Xmin = microenvironment.mesh.bounding_box[0];
	double Ymin = microenvironment.mesh.bounding_box[1];
	double Zmin = microenvironment.mesh.bounding_box[2];

	double Xmax = microenvironment.mesh.bounding_box[0];
	double Ymax = microenvironment.mesh.bounding_box[1];
	double Zmax = microenvironment.mesh.bounding_box[2];

	if (default_microenvironment_options.simulate_2D == true)
	{ 
		Zmin == 0.0;
		Zmax == 0.0;
	}

	double Xrange = Xmax - Xmin;
	double Yrange = Ymax - Ymin;
	double Zrange = Zmax - Zmin;

	double center_x = (Xmin + Xmax) * 0.5;
	double center_y = (Ymin + Ymax) * 0.5;
	double center_z = (Zmin + Zmax) * 0.5;

	#pragma omp parallel for
	for(int k=0; k<parameters.ints("number_healthy"); k++) // Place some
	{
		std::vector<double> position = {0.0, 0.0, 0.0};
		double r = NormalRandom(0,1) * tumor_radius;
		double theta = UniformRandom() * M_PI * 2;

		position[0] = center_x * cos(theta);
		position[1] = center_y * sin(theta);
		position[2] = center_z;

		pCell = create_cell();
		pCell->assign_position(position);
	}

	#pragma omp parallel for
	for(int k=0; k<parameters.ints("number_cancer"); k++)
	{
		std::vector<double> position = {0.0, 0.0, 0.0};

		position[0] = Xmin + UniformRandom() * Xrange;
		position[1] = Ymin + UniformRandom() * Yrange;
		position[2] = Zmin + UniformRandom() * Zrange;

		
		pCell = create_cell();
		pCell->assign_position(position);

	}


	return;



}

/*
	Facilitated diffusion model

	ODE model for simple diffusion based on a set of differential equations.

*/



void transport_variables_setup(Cell *pCell, Phenotype &phenotype, double dt)
{
	return;
}


void facilitated_carrier_model(Cell *pCell, Phenotype &phenotype, double dt)
{

	live.phases[0].entry_function = my_mutation_function;
	live.phase_link(0,0).arrest_function = volume_arrest_function;
	// live.display(std::cout);

	// Cell_Definition* pCD = find_cell_definition( pCell->type_name );

	static int sY_index = microenvironment.find_density_index("sY"); // microenv substrate index

	static int E_sY_near_index = pCell->custom_data.find_variable_index("E_sY_near");
	static int I_sY_index = pCell->custom_data.find_variable_index("I_sY");
	static int E_sY_index = pCell->custom_data.find_variable_index("E_sY");

	double Rc_sY_index = pCell->custom_data.find_variable_index("Rc_sY");	// [E] Receptor alone (mM)
	double Rcb_sY_index = pCell->custom_data.find_variable_index("Rcb_sY"); // [ES] Bound receptor (mM)
	double total_Rc_sY_index = pCell->custom_data.find_variable_index("total_Rc_sY");

	static int DC_sY_index = pCell->custom_data.find_variable_index("DC_sY"); // Difussion coefficient
	static int D_sY_index = pCell->custom_data.find_variable_index("D_sY");	  // Concentration gradient
	static int total_E_sY_index = pCell->custom_data.find_variable_index("total_E_sY");
	static int total_I_sY_index = pCell->custom_data.find_variable_index("total_I_sY");
	static int total_sY_index = pCell->custom_data.find_variable_index("total_sY");

	static int sY_flux_index = pCell->custom_data.find_variable_index("sY_flux");
	static int adjusted_sY_flux_index = pCell->custom_data.find_variable_index("adjusted_sY_flux");

	// Constants
	static double DC_sY = microenvironment.diffusion_coefficients[sY_index];
	static double V_cell = pCell->phenotype.volume.total;
	static double R_cell = pCell->phenotype.geometry.radius;		  // Around 8, but I don't get the difference between both
	static double R_cell_2 = std::pow(V_cell * 3.0 / 4.0, 1.0 / 3.0); // 12.3, use this one for now, check the radius in PCUG
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

	// std::cout << "Rc_sY: " << Rc_sY << std::endl;

	// Substrate concentrations
	double E_sY_near = pCell->nearest_density_vector()[sY_index];							   // A density (mM)
	double I_sY = pCell->phenotype.molecular.internalized_total_substrates[sY_index] / V_cell; // Convert to density (mM)

	static double Initial_E_sY = default_microenvironment_options.initial_condition_vector[sY_index];
	static double Initial_I_sY = parameters.doubles("Initial_I_sY");

	// Obtaining the total net internal and external amounts through iterating over the cells

	double total_E_sY = 0.0;
	double total_I_sY = 0.0;

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

	// replace hard-coded initial values at each time-step with the new values

	// Substrate concentrations
	pCell->custom_data.variables[E_sY_near_index].value = E_sY_near; // Mean of near_density_vector
	pCell->custom_data.variables[E_sY_index].value = E_sY_mean;		 // Mean of all voxels
	pCell->custom_data.variables[I_sY_index].value = I_sY;			 // Density

	pCell->custom_data.variables[DC_sY_index].value = DC_sY; // Substrate diffusion coefficient
	pCell->custom_data.variables[D_sY_index].value = D_sY;	 // Substrate diffusion coefficient

	// Receptor concentrations
	pCell->custom_data.variables[total_Rc_sY_index].value = Rcb_sY + Rc_sY;

	// Separated components

	double binding_component;
	double recycling_component;
	double endocytosis_component;
	double internal_binding_component;

	binding_component = diffusion_dt * sY_binding_rate * Rc_sY * E_sY_near;

	if (binding_component > Rc_sY * V_cell || binding_component > E_sY_near * V_cell )
	{ binding_component = std::min(Rc_sY * V_cell * diffusion_dt, E_sY_near * V_cell * diffusion_dt); } // Get the limiting step 
	
	double excess_sY = binding_component - std::min(Rc_sY * V_cell * diffusion_dt, E_sY_near * V_cell * diffusion_dt);
	{ E_sY_near += excess_sY * 1/microenvironment.mesh.dV; }
	
	Rcb_sY += binding_component;
	Rc_sY -= binding_component;

	recycling_component = diffusion_dt * sY_recycling_rate * Rcb_sY;
	if (recycling_component > Rcb_sY * V_cell){ recycling_component = Rcb_sY * V_cell * diffusion_dt; } 
	Rc_sY += recycling_component;
	Rcb_sY -= recycling_component;

	endocytosis_component = diffusion_dt * sY_endocytosis_rate * Rcb_sY;
	if (endocytosis_component > Rcb_sY * V_cell){ endocytosis_component = Rcb_sY * V_cell * diffusion_dt; } 
	I_sY += endocytosis_component;
	Rcb_sY -= endocytosis_component;
	Rc_sY += endocytosis_component;

	internal_binding_component = diffusion_dt * sY_internal_binding_rate * Rc_sY * I_sY;
	if (endocytosis_component > Rcb_sY * V_cell || endocytosis_component > I_sY * V_cell )
	{ endocytosis_component = std::min(Rcb_sY * V_cell * diffusion_dt, I_sY * V_cell * diffusion_dt); } // Get the limiting step 
	
	I_sY -= internal_binding_component;
	Rcb_sY += internal_binding_component;
	Rc_sY -= internal_binding_component;

	double vf_max = sY_endocytosis_rate * (Rc_sY + Rcb_sY); // both are constant given that total receptor is constant
	double vb_max = sY_recycling_rate * (Rc_sY + Rcb_sY);
	double Km_1 = (sY_recycling_rate + sY_endocytosis_rate) / sY_binding_rate;
	double Km_2 = (sY_recycling_rate + sY_endocytosis_rate) / sY_internal_binding_rate;

	double v_formation = -((vf_max * E_sY_near) / Km_1 - (vb_max * I_sY) / Km_2) / (1 + (I_sY / Km_2) + (E_sY_near / Km_1)); // Option A
	// double v_formation = - binding_component + recycling_component - endocytosis_component + internal_binding_component;  // Option B
	// double v_formation = - endocytosis_component;																		 // Option C

	v_formation = v_formation * V_cell; // Convert mM/min to amol/min

	pCell->custom_data.variables[sY_flux_index].value = v_formation;
	pCell->custom_data.variables[Rc_sY_index].value = Rc_sY;
	pCell->custom_data.variables[Rcb_sY_index].value = Rcb_sY;
	pCell->phenotype.secretion.net_export_rates[sY_index] = v_formation / diffusion_dt; // Avoid second multiplication by diffusion_dt

	// Net total amounts
	pCell->custom_data.variables[total_E_sY_index].value = total_E_sY;
	pCell->custom_data.variables[total_I_sY_index].value = total_I_sY;
	pCell->custom_data.variables[total_sY_index].value = total_E_sY + total_I_sY;

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

	// pCell->phenotype.death.rates[ necrosis_idx ] *= I_sY

	// std::cout << "I_sY: " << I_sY << std::endl;



	// Connecting the transport model to the phenotype

	if (I_sY > parameters.doubles("cytotoxic_threshold") && parameters.doubles("cytotoxic_threshold") != 0.0)
	{
		std::cout << "Necrosis!" << std::endl;
		// pCell->start_death(necrosis_idx);
		pCell->start_death( apoptosis_idx ); // Apparently no difference between necrosis and apoptosis?
	}


	if (I_sY > parameters.doubles("growth_threshold"))
	{
		// pCell->flag_for_division();
		// std::cout << "Growth! I_sY: " << I_sY << std::endl;
		// pCell->phenotype.cycle.advance_cycle(pCell, phenotype, diffusion_dt);
		// std::cout << "RATE before: " << pCell->phenotype.cycle.data.transition_rate(0,0) << std::endl
		// << "I_sY: " << I_sY << std::endl;
		pCell->phenotype.cycle.data.transition_rate(0,0) = I_sY * 0.0001;
		// std::cout << "RATE after: " << pCell->phenotype.cycle.data.transition_rate(0,0) << std::endl;

	}

	
	if (I_sY < parameters.doubles("growth_threshold"))
	{
		// pCell->flag_for_division();
		// pCell->phenotype.cycle.advance_cycle(pCell, phenotype, diffusion_dt);
		pCell->phenotype.cycle.data.transition_rate(0,0) = 1e-08;

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

// Both functions related to the BoSS part of the simulation, leave for now
void set_input_nodes(Cell *pCell) {}

void from_nodes_to_cell(Cell *pCell, Phenotype &phenotype, double dt) {}




// Specific coloring function - specific for heterogeneity

std::vector<std::string> my_coloring_function(Cell *pCell)
{

	// immune are blue
	// std::vector< std::string > output( 4, "black" );
	
	std::vector<std::string> output = false_cell_coloring_live_dead(pCell);

	// if (pCell->type == 1)
	// {
	// 	output[0] = "blue";
	// 	output[1] = "blue";
	// 	output[2] = "blue";
	// 	output[3] = "blue";

	// 	return output;
	// }

	// if (pCell->type == 2)
	// {
	// 	output[0] = "red";
	// 	output[1] = "red";
	// 	output[2] = "red";
	// 	output[3] = "red";

	// 	return output;
	// }

	// live cells are green, but shaded by oncoprotein value
	// if (pCell->phenotype.death.dead == false)
	// {
	// 	// int oncoprotein = (int) round( (1.0/(p_max- p_min + 10.0)) * (pCell->custom_data[oncoprotein_i]-p_min) * 255.0 );
	// 	// int internal_sY = (int) round( ( std::pow( (pCell->custom_data[3] * 1000) / 255, 2.0 )) );
	// 	int internal_sY = (int)round(pCell->custom_data[1] * 1000);
	// 	char szTempString[128];
	// 	sprintf(szTempString, "rgb(%u,%u,%u)", 138, 184, internal_sY);
	// 	output[0].assign(szTempString);
	// 	output[1].assign(szTempString);
	// 	output[2].assign("black");
	// 	output[3].assign(szTempString);
	// }

        
	// // if not, dead colors

	// if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic) // Apoptotic - Red
	// {
	// 	output[0] = "rgb(255,0,0)";
	// 	output[2] = "rgb(125,0,0)";
	// }

	// // Necrotic - Brown
	// if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling ||
	// 	pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed ||
	// 	pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic)
	// {
	// 	output[0] = "rgb(0,0,0)";
	// 	output[2] = "rgb(139,69,19)";
	// }

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

void my_mutation_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Fetch from custom data instead

	double p_mean = 0.05;
	double p_sd = 0.1;
	double p_min = 0.0;
	double p_max = 4.0;

	pCell->custom_data[1] = 0.0;

	pCell->custom_data[4] = NormalRandom(p_mean, p_sd);
	if (pCell->custom_data[4] < p_min)
	{ pCell->custom_data[4] = p_min; }
	if (pCell->custom_data[4] > p_max)
	{ pCell->custom_data[4] = p_max; }

	return;
}

bool volume_arrest_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.volume.total < 2494 )
	{ return true; }
	return false;
}

