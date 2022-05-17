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


	static int sY_index = microenvironment.find_density_index("sY");
	static int cytB_index = microenvironment.find_density_index("cytB");

	Cell *pCell;

	// Generic variables:
	cell_defaults.custom_data.add_variable(parameters.strings("node_to_visualize"), "dimensionless", 0.0); // for paraview visualization

	// Densities - sY
	cell_defaults.custom_data.add_variable("I_sY", "amol / um^3", parameters.doubles("Initial_I_sY"));
	cell_defaults.custom_data.add_variable("E_sY_near", "amol / um^3", default_microenvironment_options.initial_condition_vector[sY_index]); // Mean nearest substrate conc.
	cell_defaults.custom_data.add_variable("E_sY", "amol/um^3", default_microenvironment_options.initial_condition_vector[sY_index]);		 // Mean external substrate conc.

	// Densities - cytB
	cell_defaults.custom_data.add_variable("I_cytB", "amol / um^3", parameters.doubles("Initial_I_cytB"));
	cell_defaults.custom_data.add_variable("E_cytB_near", "amol / um^3", default_microenvironment_options.initial_condition_vector[cytB_index]); // Mean nearest substrate conc.
	cell_defaults.custom_data.add_variable("E_cytB", "amol/um^3", default_microenvironment_options.initial_condition_vector[cytB_index]);

	// Receptors - sY (FDC)
	cell_defaults.custom_data.add_variable("Rc_sY", "amol/um^3", parameters.doubles("Initial_Rc_sY"));
	cell_defaults.custom_data.add_variable("Rcb_sY", "amol / um^3", parameters.doubles("Initial_Rcb_sY"));
	cell_defaults.custom_data.add_variable("total_Rc_sY", "amol / um^3", parameters.doubles("Initial_Rcb_sY") + parameters.doubles("Initial_Rc_sY")); // Bound receptor density

	// Densities - ATP
	cell_defaults.custom_data.add_variable("I_ATP", "amol / um^3", parameters.doubles("Initial_I_ATP"));
	cell_defaults.custom_data.add_variable("I_ADP", "amol / um^3", parameters.doubles("Initial_I_ADP")); 
	
	// Receptors - cytB, ATP 
	cell_defaults.custom_data.add_variable("Rcb_cytB", "amol/um^3", parameters.doubles("Initial_Rcb_cytB"));
	cell_defaults.custom_data.add_variable("Rcb_ATP", "amol/um^3", parameters.doubles("Initial_Rcb_ATP"));
	cell_defaults.custom_data.add_variable("total_Rcb_cytB", "amol / um^3", parameters.doubles("Initial_Rcb_cytB") + parameters.doubles("Initial_Rcb_ATP")); // Bound receptor density
	cell_defaults.custom_data.add_variable("Rc_cytB", "amol/um^3", parameters.doubles("Initial_Rc_cytB"));
	cell_defaults.custom_data.add_variable("total_Rc_cytB", "amol/um^3", parameters.doubles("Initial_Rc_cytB")); // Total receptor density


	// Total amounts - sY
	cell_defaults.custom_data.add_variable("total_E_sY", "amol", 0.0); // total external sY
	cell_defaults.custom_data.add_variable("total_I_sY", "amol", 0.0); // total internal sY
	cell_defaults.custom_data.add_variable("total_sY", "amol", 0.0);   // total sY

	// Total amounts - cytB
	cell_defaults.custom_data.add_variable("total_E_cytB", "amol", 0.0); // total external sY
	cell_defaults.custom_data.add_variable("total_I_cytB", "amol", 0.0); // total internal sY
	cell_defaults.custom_data.add_variable("total_cytB", "amol", 0.0);   // total sY

	// Flux - sY
	cell_defaults.custom_data.add_variable("sY_flux", "amol/min", 0.0);
	cell_defaults.custom_data.add_variable("adjusted_sY_flux", "amol/min", 0.0);

	// Flux - cytB
	cell_defaults.custom_data.add_variable("cytB_flux", "amol/min", 0.0);
	cell_defaults.custom_data.add_variable("adjusted_cytB_flux", "amol/min", 0.0);

	cell_defaults.custom_data.add_variable("pump_rate", "amol/min", 0.0);
	cell_defaults.custom_data.add_variable("pumping_probabilty", "dimensionless", 0.0); 



	// Kinetic parameters - ATP-pump
	cell_defaults.custom_data.add_variable("cytB_binding_component", "mmol/um^3", 0.0);
	cell_defaults.custom_data.add_variable("cytB_movement_component", "mmol/um^3", 0.0);
	cell_defaults.custom_data.add_variable("cytB_internal_binding_component", "mmol/um^3", 0.0);
	cell_defaults.custom_data.add_variable("cytB_recycling_component", "mmol/um^3", 0.0);

	// Kinetic parameters - GLUT carrier
	cell_defaults.custom_data.add_variable("binding_component", "mmol/um^3", 0.0);
	cell_defaults.custom_data.add_variable("endocytosis_component", "mmol/um^3", 0.0);
	cell_defaults.custom_data.add_variable("internal_binding_component", "mmol/um^3", 0.0);
	cell_defaults.custom_data.add_variable("recycling_component", "mmol/um^3", 0.0);
	
	// Kinetic parameters - ATP-ADP system
	cell_defaults.custom_data.add_variable("ATP_binding_component", "mmol/um^3", 0.0);
	cell_defaults.custom_data.add_variable("ATP_unbinding_component", "mmol/um^3", 0.0);
	cell_defaults.custom_data.add_variable("ATP_hydrolisis_component", "mmol/um^3", 0.0);
	cell_defaults.custom_data.add_variable("ADP_recycling_component", "mmol/um^3", 0.0);

	// Diffusion coefficients, kinetic parameters - sY
	cell_defaults.custom_data.add_variable("DC_sY", "um^3/min", microenvironment.diffusion_coefficients[sY_index]); // Diffusion coefficient
	cell_defaults.custom_data.add_variable("D_sY", "mmol/um^3", 0.0);

	// If we want to add noise to the kinetic parameters, we have to do it here
	// double noise_amount = 0.0;
	// cell_defaults.custom_data.add_variable("k_1", "um^2/min", NormalRandom(parameters.doubles("k_1"), noise_amount) );
	// cell_defaults.custom_data.add_variable("k_n1", "um^2/min", NormalRandom(parameters.doubles("k_n1"), noise_amount) );
	// cell_defaults.custom_data.add_variable("k_2", "um^2/min", NormalRandom(parameters.doubles("k_2"), noise_amount) );
	// cell_defaults.custom_data.add_variable("k_n2", "um^2/min", NormalRandom(parameters.doubles("k_n2"), noise_amount) );

	// std::cout << "k_1 = " << parameters.doubles("k_1") << std::endl;


	// Diffusion coefficients, kinetic parameters - cytB
	cell_defaults.custom_data.add_variable("DC_cytB", "um^3/min", microenvironment.diffusion_coefficients[cytB_index]); // Diffusion coefficient
	cell_defaults.custom_data.add_variable("D_cytB", "mmol/um^3", 0.0);

	cell_defaults.custom_data.add_variable("k_cytB", "um^2/min", parameters.doubles("k_cytB") );

	
	// Other constants (Num. of vocels, Initial conditions, etc.)


	cell_defaults.custom_data.add_variable("Initial_E_sY", "amol/um^3", default_microenvironment_options.initial_condition_vector[sY_index]); // Initial external substrate_Y
	cell_defaults.custom_data.add_variable("Initial_I_sY", "amol/um^3", parameters.doubles("Initial_I_sY"));	
	
	cell_defaults.custom_data.add_variable("Initial_E_cytB", "amol/um^3", default_microenvironment_options.initial_condition_vector[cytB_index]); 
	cell_defaults.custom_data.add_variable("Initial_I_cytB", "amol/um^3", parameters.doubles("Initial_I_cytB"));

	cell_defaults.custom_data.add_variable("Initial_I_ATP", "amol/um^3", parameters.doubles("Initial_I_ATP"));
	cell_defaults.custom_data.add_variable("Initial_I_ADP", "amol/um^3", parameters.doubles("Initial_I_ADP"));


	cell_defaults.custom_data.add_variable("total_healthy_cells", "dimensionless", 0.0);				
	cell_defaults.custom_data.add_variable("total_cancer_cells", "dimensionless", 0.0);
	cell_defaults.custom_data.add_variable("number_of_cells", "dimensionless", (*all_cells).size());

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

	// Cell_Definiton healthy; // Healthy cell
	// Cell_Definition cancer_cell; // Cancer cell

	static Cell_Definition* pCH = find_cell_definition("healthy");
	static Cell_Definition* pCC = find_cell_definition("cancer_cell");
	

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


// Setup tissue from spheroid TNF - A spheroid with 100 cells

/*
void setup_tissue(void)
{

	double R_mean = parameters.doubles("R_mean");
	double R_sd = parameters.doubles("R_sd");
	double R_min = parameters.doubles("R_min");
	double R_max = parameters.doubles("R_max");

	double Initial_I_sY = parameters.doubles("Initial_I_sY");
	static int sY_index = microenvironment.find_density_index("sY");

	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);

	for (int i = 0; i < cells.size(); i++)
	{
		std::cout << "Placing cancer cells.. " << i << std::endl;

		Cell *pCell;
		float x = cells[i].x;
		float y = cells[i].y;
		float z = cells[i].z;
		// double elapsed_time = cells[i].elapsed_time;
		
		pCell = create_cell(get_cell_definition("cancer_cell"));
		// pCell->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;

		pCell->phenotype.molecular.internalized_total_substrates[sY_index] = parameters.doubles("Initial_I_sY");

		pCell->custom_data[7] = NormalRandom(R_mean, R_sd) + 0.5; // Distro of GLUT transporter density
		if (pCell->custom_data[7] < R_min) { pCell->custom_data[7] = R_min; }
		if (pCell->custom_data[7] > R_max) { pCell->custom_data[7] = R_max; }

		std::cout << "With [Rc]: " << pCell->custom_data[7] << std::endl;

		// Probability of developing MDR resistance
		double pumping_probability = std::abs(NormalRandom(0, 0.5));
		if (pumping_probability > 1.0){ pumping_probability = 1.0; }

		pCell->assign_position(x, y, z);

		// update_monitor_variables(pCell);
	}

	return;
}

*/


void setup_tissue(void)
{
	Cell *pCell = NULL;
	static Cell_Definition* pCH = find_cell_definition("healthy");
	static Cell_Definition* pCC = find_cell_definition("cancer_cell");

	double Initial_Rc_sY = parameters.doubles("Initial_Rc_sY");
	double Initial_Rcb_sY = parameters.doubles("Initial_Rcb_sY");

	double Initial_I_sY = parameters.doubles("Initial_I_sY");
	static int sY_index = microenvironment.find_density_index("sY");

	double cell_radius = cell_defaults.phenotype.geometry.radius;
	double cell_spacing = 0.95 * 2.2 * cell_radius;
	double tumor_radius = parameters.doubles("tumor_radius"); // 250.0;

	double R_mean = parameters.doubles("R_mean");
	double R_sd = parameters.doubles("R_sd");
	double R_min = parameters.doubles("R_min");
	double R_max = parameters.doubles("R_max");

	
	double Xmin = microenvironment.mesh.bounding_box[0];
	double Ymin = microenvironment.mesh.bounding_box[1];
	double Zmin = microenvironment.mesh.bounding_box[2];

	double Xmax = - microenvironment.mesh.bounding_box[0]; // negative sign bc if not range is 0
	double Ymax = - microenvironment.mesh.bounding_box[1];
	double Zmax = - microenvironment.mesh.bounding_box[2];

	if (default_microenvironment_options.simulate_2D == true)
	{ 
		Zmin == 0.0;
		Zmax == 0.0;
	}

	double Xrange = Xmax - Xmin; // - 200.0 before
	double Yrange = Ymax - Ymin;
	double Zrange = Zmax - Zmin;

	std::cout << "Xmin: " << Xmin << std::endl;
	std::cout << "Ymin: " << Ymin << std::endl;
	std::cout << "Zmin: " << Zmin << std::endl;
	std::cout << "Xmax: " << Xmax << std::endl;
	std::cout << "Ymax: " << Ymax << std::endl;
	std::cout << "Zmax: " << Zmax << std::endl;
	std::cout << "Xrange: " << Xrange << std::endl;
	std::cout << "Yrange: " << Yrange << std::endl;
	std::cout << "Zrange: " << Zrange << std::endl;

	double center_x = (Xmin + Xmax) * 0.5;
	double center_y = (Ymin + Ymax) * 0.5;
	double center_z = (Zmin + Zmax) * 0.5;


	// OPTION A: Two clusters, one of healthy cells, the other cancer cells

	// Healthy cell cluster, centered at (-250, 0)

	double x = 0.0;
	double x_outer = tumor_radius;
	double y = 0.0;

	// int n = 0;
	// while (y < tumor_radius)
	// {
	// 	x = 0.0;
	// 	if (n % 2 == 1)
	// 	{
	// 		x = 0.5 * cell_spacing;
	// 	}
	// 	x_outer = sqrt(tumor_radius * tumor_radius - y * y);

	// 	while (x < x_outer)
	// 	{
	// 		pCell = create_cell( *pCH );
	// 		pCell->assign_position(x, y, 0.0);
	// 		pCell->phenotype.molecular.internalized_total_substrates[sY_index] = parameters.doubles("Initial_I_sY");// NormalRandom(parameters.doubles("Initial_I_sY"), 0.01); 
	// 		if (pCell->phenotype.molecular.internalized_total_substrates[sY_index] < 0.0) { pCell->phenotype.molecular.internalized_total_substrates[sY_index] = 0.0; }

	// 		pCell->custom_data[7] = NormalRandom(R_mean, R_sd); // Random distribution of transporter
	// 		if (pCell->custom_data[7] < R_min) { pCell->custom_data[7] = R_min; }
	// 		if (pCell->custom_data[7] > R_max) { pCell->custom_data[7] = R_max; }

	// 		if (fabs(y) > 0.01)
	// 		{
	// 			pCell = create_cell( *pCH ); // tumor cell
	// 			pCell->assign_position(x, -y, 0.0);
	// 			pCell->phenotype.molecular.internalized_total_substrates[sY_index] = parameters.doubles("Initial_I_sY");// NormalRandom(parameters.doubles("Initial_I_sY"), 0.01); 
	// 			if (pCell->phenotype.molecular.internalized_total_substrates[sY_index] < 0.0) { pCell->phenotype.molecular.internalized_total_substrates[sY_index] = 0.0; }

	// 			pCell->custom_data[7] = NormalRandom(R_mean, R_sd); // Random distribution of transporter
	// 			if (pCell->custom_data[7] < R_min) { pCell->custom_data[7] = R_min; }
	// 			if (pCell->custom_data[7] > R_max) { pCell->custom_data[7] = R_max; }
	// 		}

	// 		if (fabs(x) > 0.01)
	// 		{
	// 			pCell = create_cell( *pCH ); // tumor cell
	// 			pCell->assign_position(-x, y, 0.0);
	// 			pCell->phenotype.molecular.internalized_total_substrates[sY_index] = parameters.doubles("Initial_I_sY");// NormalRandom(parameters.doubles("Initial_I_sY"), 0.01); 
	// 			if (pCell->phenotype.molecular.internalized_total_substrates[sY_index] < 0.0) { pCell->phenotype.molecular.internalized_total_substrates[sY_index] = 0.0; }

	// 			pCell->custom_data[7] = NormalRandom(R_mean, R_sd); // Random distribution of transporter
	// 			if (pCell->custom_data[7] < R_min) { pCell->custom_data[7] = R_min; }
	// 			if (pCell->custom_data[7] > R_max) { pCell->custom_data[7] = R_max; }

	// 			if (fabs(y) > 0.01)
	// 			{
	// 				pCell = create_cell( *pCH ); // tumor cell
	// 				pCell->assign_position(-x, -y, 0.0);
					
	// 				pCell->phenotype.molecular.internalized_total_substrates[sY_index] = parameters.doubles("Initial_I_sY");// NormalRandom(parameters.doubles("Initial_I_sY"), 0.01); 
	// 				if (pCell->phenotype.molecular.internalized_total_substrates[sY_index] < 0.0) { pCell->phenotype.molecular.internalized_total_substrates[sY_index] = 0.0; }

	// 				pCell->custom_data[7] = NormalRandom(R_mean, R_sd); // Random distribution of transporter
	// 				if (pCell->custom_data[7] < R_min) { pCell->custom_data[7] = R_min; }
	// 				if (pCell->custom_data[7] > R_max) { pCell->custom_data[7] = R_max; }
	// 			}
	// 		}

	// 		x += cell_spacing;
	// 	}

	// 	y += cell_spacing * sqrt(3.0) / 2.0;
	// 	n++;
	// }


	

	// for(int k=0; k<parameters.ints("number_cancer"); k++)
	// {

	// 	std::cout << "Placing cancer cells.. " << k << std::endl;

	// 	std::vector<double> position = {0.0, 0.0, 0.0};
	// 	double r = NormalRandom(0, 1) * tumor_radius;
	// 	double theta = UniformRandom() * M_PI * 2;

	// 	position[0] = center_x + r*std::cos(theta);
	// 	position[1] = center_y + r*std::sin(theta);
	// 	position[2] = center_z + r*UniformRandom();

	// 	pCell = create_cell( *pCC );
	// 	pCell->assign_position( position );
      
	// 	pCell->phenotype.molecular.internalized_total_substrates[sY_index] = parameters.doubles("Initial_I_sY");

	// 	pCell->custom_data[7] = NormalRandom(R_mean, R_sd) + 0.5; // Random distro of GLUT transporter
	// 	if (pCell->custom_data[7] < R_min) { pCell->custom_data[7] = R_min; }
	// 	if (pCell->custom_data[7] > R_max) { pCell->custom_data[7] = R_max; }

	// 	// Probability of developing MDR resistance
	// 	double pumping_probability = std::abs(NormalRandom(0, 0.5));
	// 	if (pumping_probability > 1.0){ pumping_probability *= 0.01; }

	// 	pCell->custom_data[25] = pumping_probability;
	// 	std::cout << "Pumping probability: " << pumping_probability << std::endl;

	// }

	// double tumor_radius = 250.0; // 250.0; 
	
	for(int k=0; k<parameters.ints("number_cancer"); k++)
	{
		double theta = UniformRandom() * 6.283185307179586476925286766559; 
		double phi = acos( 2.0*UniformRandom() - 1.0 );  
		
		double radius = NormalRandom( tumor_radius, 25.0 ); 
		
		pCell = create_cell( *pCC );
		pCell->assign_position( radius*cos(theta)*sin(phi), radius*sin(theta)*sin(phi), 0.0 ); 
	}

	return;
}



	



/*
	Facilitated diffusion model + Simple diffusion of cytochalasin B + ATP Pumping to recover tumoral phenotype

*/

void HC_cell_rule( Cell *pCell, Phenotype &phenotype, double dt )
{
	static int sY_index = microenvironment.find_density_index( "sY" );
	double V_cell = pCell->phenotype.volume.total;
	double I_sY = pCell->phenotype.molecular.internalized_total_substrates[sY_index] / V_cell;

	double I_sY_degradation = parameters.doubles("k_3") * I_sY * diffusion_dt;

	return;

}

void CC_cell_rule( Cell *pCell, Phenotype &phenotype, double dt )
{
	static int sY_index = microenvironment.find_density_index( "sY" );
	double V_cell = pCell->phenotype.volume.total;
	double I_sY = pCell->phenotype.molecular.internalized_total_substrates[sY_index] / V_cell;
	
	double I_sY_degradation = 2 * parameters.doubles("k_3") * I_sY * diffusion_dt;

	return;

}

void facilitated_carrier_model(Cell *pCell, Phenotype &phenotype, double dt)
{
	Cell_Definition* pCD = find_cell_definition( pCell->type_name );

	// live.phases[0].entry_function = my_mutation_function;
	// live.phase_link(0,0).arrest_function = volume_arrest_function;
	// live.display(std::cout);

	static int sY_index = microenvironment.find_density_index("sY"); // microenv substrate index
	static int cytB_index = microenvironment.find_density_index("cytB"); 

	// Densities - sY
	static int I_sY_index = pCell->custom_data.find_variable_index("I_sY");
	static int E_sY_near_index = pCell->custom_data.find_variable_index("E_sY_near");
	static int E_sY_index = pCell->custom_data.find_variable_index("E_sY");

	// Densities - cytB
	static int I_cytB_index = pCell->custom_data.find_variable_index("I_cytB");
	static int E_cytB_near_index = pCell->custom_data.find_variable_index("E_cytB_near");
	static int E_cytB_index = pCell->custom_data.find_variable_index("E_cytB");

	// Densities - ATP
	static int I_ATP_index = pCell->custom_data.find_variable_index("I_ATP");
	static int I_ADP_index = pCell->custom_data.find_variable_index("I_ADP");

	// Receptors - sY
	double Rc_sY_index = pCell->custom_data.find_variable_index("Rc_sY");	// [E] Receptor alone (mM)
	double Rcb_sY_index = pCell->custom_data.find_variable_index("Rcb_sY"); // [ES] Bound receptor (mM)
	double total_Rc_sY_index = pCell->custom_data.find_variable_index("total_Rc_sY");

	// Receptors - ATP-pump of cytB
	double Rcb_cytB_index = pCell->custom_data.find_variable_index("Rcb_cytB");
	double Rcb_ATP_index = pCell->custom_data.find_variable_index("Rcb_ATP");
	double Rc_cytB_index = pCell->custom_data.find_variable_index("Rc_cytB");
	double total_Rcb_cytB_index = pCell->custom_data.find_variable_index("total_Rcb_cytB");
	double total_Rc_cytB_index = pCell->custom_data.find_variable_index("total_Rc_cytB");

	static int DC_sY_index = pCell->custom_data.find_variable_index("DC_sY"); // Difussion coefficient
	static int D_sY_index = pCell->custom_data.find_variable_index("D_sY");	  // Concentration gradient

	static int DC_cytB_index = pCell->custom_data.find_variable_index("DC_cytB"); // Difussion coefficient
	static int D_cytB_index = pCell->custom_data.find_variable_index("D_cytB");	  // Concentration gradient

	// static int sY_k1_index = pCell->custom_data.find_variable_index("k_1");	  
	// static int sY_kn1_index = pCell->custom_data.find_variable_index("k_n1");	  
	// static int sY_k2_index = pCell->custom_data.find_variable_index("k_2");	 
	// static int sY_kn2_index = pCell->custom_data.find_variable_index("k_n2");

	static int total_E_sY_index = pCell->custom_data.find_variable_index("total_E_sY");
	static int total_I_sY_index = pCell->custom_data.find_variable_index("total_I_sY");
	static int total_sY_index = pCell->custom_data.find_variable_index("total_sY");

	static int total_E_cytB_index = pCell->custom_data.find_variable_index("total_E_cytB");
	static int total_I_cytB_index = pCell->custom_data.find_variable_index("total_I_cytB");
	static int total_cytB_index = pCell->custom_data.find_variable_index("total_cytB");

	static int sY_flux_index = pCell->custom_data.find_variable_index("sY_flux");
	static int adjusted_sY_flux_index = pCell->custom_data.find_variable_index("adjusted_sY_flux");

	static int cytB_flux_index = pCell->custom_data.find_variable_index("cytB_flux");
	static int adjusted_cytB_flux_index = pCell->custom_data.find_variable_index("adjusted_cytB_flux");

	static int pump_rate_index = pCell->custom_data.find_variable_index("pump_rate");
	static int pumping_probability_index = pCell->custom_data.find_variable_index("pumping_probability");

	static int total_healthy_cells_index = pCell->custom_data.find_variable_index("total_healthy_cells");
	static int total_cancer_cells_index = pCell->custom_data.find_variable_index("total_cancer_cells");

	static int time_index = pCell->custom_data.find_variable_index("time");
	pCell->custom_data.variables[time_index].value = PhysiCell_globals.current_time;
	static int total_cells_index = pCell->custom_data.find_variable_index("number_of_cells");
	pCell->custom_data.variables[total_cells_index].value = (*all_cells).size();


	// Constants
	static double DC_sY = microenvironment.diffusion_coefficients[sY_index];
	static double DC_cytB = microenvironment.diffusion_coefficients[cytB_index];

	// Kinetic parameters - cytB-ATP-pump
	static int cytB_binding_component_index = pCell->custom_data.find_variable_index("cytB_binding_component");
	static int cytB_movement_component_index = pCell->custom_data.find_variable_index("cytB_movement_component");
	static int cytB_internal_binding_component_index = pCell->custom_data.find_variable_index("cytB_internal_binding_component");
	static int cytB_recycling_component_index = pCell->custom_data.find_variable_index("cytB_recycling_component");

	// Kinetic parameters - GLUT Carrier
	static int binding_component_index = pCell->custom_data.find_variable_index("binding_component");
	static int endocytosis_component_index = pCell->custom_data.find_variable_index("endocytosis_component");
	static int internal_binding_component_index = pCell->custom_data.find_variable_index("internal_binding_component");
	static int recycling_component_index = pCell->custom_data.find_variable_index("recycling_component");

	// Kinetic parameters - ATP-ADP system
	static int ATP_binding_component_index = pCell->custom_data.find_variable_index("ATP_binding_component");
	static int ATP_unbinding_component_index = pCell->custom_data.find_variable_index("ATP_unbinding_component");
	static int ATP_hydrolisis_component_index = pCell->custom_data.find_variable_index("ATP_hydrolisis_component");
	static int ADP_recycling_component_index = pCell->custom_data.find_variable_index("ADP_recycling_component");


	double V_cell = pCell->phenotype.volume.total;
	double R_cell = pCell->phenotype.geometry.radius;		
	double A_cell = 4 * M_PI * std::pow(R_cell, 2);

	double V_voxel = (microenvironment.mesh.dV);
	double total_num_voxels = microenvironment.number_of_voxels();


	/*

		Kinetic parameters of sY entry if we want to add noise to the sY flux


	if( pCell->custom_data.variables[sY_k1_index].value < 0.0) { pCell->custom_data.variables[sY_k1_index].value = 0.0; }
	if( pCell->custom_data.variables[sY_kn1_index].value < 0.0) { pCell->custom_data.variables[sY_kn1_index].value = 0.0; }
	if( pCell->custom_data.variables[sY_k2_index].value < 0.0) { pCell->custom_data.variables[sY_k2_index].value = 0.0; }
	if( pCell->custom_data.variables[sY_kn2_index].value < 0.0) { pCell->custom_data.variables[sY_kn2_index].value = 0.0; }

	double sY_binding_rate = pCell->custom_data.variables[sY_k1_index].value; // parameters.doubles("k_1");
	double sY_recycling_rate = pCell->custom_data.variables[sY_kn2_index].value; //parameters.doubles("k_n1");
	double sY_endocytosis_rate = pCell->custom_data.variables[sY_k2_index].value; //parameters.doubles("k_2");
	double sY_internal_binding_rate = pCell->custom_data.variables[sY_kn2_index].value; //parameters.doubles("k_n2");

	double sY_internal_degradation_rate = parameters.doubles("k_i");

	pCell->custom_data.variables[sY_k1_index].value = sY_binding_rate;
	pCell->custom_data.variables[sY_kn1_index].value = sY_recycling_rate;
	pCell->custom_data.variables[sY_k2_index].value = sY_endocytosis_rate;
	pCell->custom_data.variables[sY_kn2_index].value = sY_internal_binding_rate;

	*/

	// Kinetic parameters of sY entry without noise

	double sY_binding_rate = parameters.doubles("k_1"); 
	double sY_recycling_rate = parameters.doubles("k_n1"); 
	double sY_endocytosis_rate = parameters.doubles("k_2");
	double sY_internal_binding_rate = parameters.doubles("k_n2");

	// cytB binding to ATP pump kinetics
	double cytB_binding_rate = parameters.doubles("cytB_k_1B");
	double cytB_recycling_rate = parameters.doubles("cytB_k_n1B");
	double cytB_movement_rate = parameters.doubles("cytB_k_2B");
	double cytB_internal_binding_rate = parameters.doubles("cytB_k_n2B"); // not employed, the transport is one-directional

	// ATP-ADP kinetic parameters
	double ATP_binding_rate = parameters.doubles("k_1A");
	double ATP_recycling_rate = parameters.doubles("k_n1A");
	double ATP_hydrolisis_rate = parameters.doubles("k_2A");
	double ADP_recycling_rate = parameters.doubles("k_n2A");

	// Receptor concentrations
	double Rc_sY = pCell->custom_data.variables[Rc_sY_index].value;	  // [E] Receptor alone (mM)
	double Rcb_sY = pCell->custom_data.variables[Rcb_sY_index].value; // [ES] Bound receptor (mM)

	// Receptor densities - ATP pump
	double Rcb_cytB = pCell->custom_data.variables[Rcb_cytB_index].value;
	double Rcb_ATP = pCell->custom_data.variables[Rcb_ATP_index].value;
	double total_Rcb_cytB = pCell->custom_data.variables[total_Rcb_cytB_index].value;
	double Rc_cytB = pCell->custom_data.variables[Rc_cytB_index].value;

	// Substrate concentrations
	double E_sY_near = pCell->nearest_density_vector()[sY_index];							   // A density (mM)
	double I_sY = pCell->phenotype.molecular.internalized_total_substrates[sY_index] / V_cell; // Convert to density (mM)
	static double Initial_E_sY = default_microenvironment_options.initial_condition_vector[sY_index];
	static double Initial_I_sY = parameters.doubles("Initial_I_sY");

	double E_cytB_near = pCell->nearest_density_vector()[cytB_index];							   // A density (mM)
	double I_cytB = pCell->phenotype.molecular.internalized_total_substrates[cytB_index] / V_cell; // Convert to density (mM)
	static double Initial_E_cytB = default_microenvironment_options.initial_condition_vector[cytB_index];
	static double Initial_I_cytB = parameters.doubles("Initial_I_cytB");

	double I_ATP = parameters.doubles("Initial_I_ATP");
	double I_ADP = parameters.doubles("Initial_I_ADP");


	// Obtaining the total net internal and external amounts through iterating over the cells

	double total_E_sY = 0.0;
	double total_I_sY = 0.0;
	double total_E_cytB = 0.0;
	double total_I_cytB = 0.0;

	std::vector<double> I_sY_vector(microenvironment.number_of_voxels());
	std::vector<double> I_cytB_vector(microenvironment.number_of_voxels());

	int total_healthy_cells = 0.0;
	int total_cancer_cells = 0.0;

	// Counting net internal amount
	#pragma omp parallel for
	for (int i = 0; i < (*all_cells).size(); i++)
	{
		Cell *pC = (*all_cells)[i];
		I_sY_vector[i] = pC->phenotype.molecular.internalized_total_substrates[sY_index];
		total_I_sY += pC->phenotype.molecular.internalized_total_substrates[sY_index];

		I_cytB_vector[i] = pC->phenotype.molecular.internalized_total_substrates[cytB_index];
		total_I_cytB += pC->phenotype.molecular.internalized_total_substrates[cytB_index];

		// Get total healthy and cancer cells
		if(pC->type == 1){ total_healthy_cells +=1; }
		if(pC->type == 2){ total_cancer_cells += 1; }

	}

	// std::cout << total_cancer_cells << " cancer cells" << std::endl;
	pCell->custom_data.variables[total_cancer_cells_index].value = total_cancer_cells;
	pCell->custom_data.variables[total_healthy_cells_index].value = total_healthy_cells;

	total_cancer_cells = 0.0;
	total_healthy_cells = 0.0;

	double I_sY_sum = std::accumulate(I_sY_vector.begin(), I_sY_vector.end(), 0.0);
	double I_cytB_sum = std::accumulate(I_cytB_vector.begin(), I_cytB_vector.end(), 0.0);

	// I_sY_sum and total_I_sY provide the same exact result. They reflect the actual total amount of
	// internal net sY

	std::vector<double> E_sY_vector(microenvironment.number_of_voxels());
	std::vector<double> E_cytB_vector(microenvironment.number_of_voxels());

	#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		E_sY_vector[n] = (microenvironment.density_vector(n)[sY_index]);
		total_E_sY += (microenvironment.density_vector(n)[sY_index]) * V_voxel; // A net amount of sY

		E_cytB_vector[n] = (microenvironment.density_vector(n)[cytB_index]);
		total_E_cytB += (microenvironment.density_vector(n)[cytB_index]) * V_voxel; 
	}

	double E_sY_sum = std::accumulate(E_sY_vector.begin(), E_sY_vector.end(), 0.0);
	double E_sY_mean = E_sY_sum / total_num_voxels; // NOT CORRECT, the accumulate just adds up the number of voxels

	double E_cytB_sum = std::accumulate(E_cytB_vector.begin(), E_cytB_vector.end(), 0.0);
	double E_cytB_mean = E_cytB_sum / total_num_voxels; // NOT CORRECT, the accumulate just adds up the number of voxels

	double D_sY = I_sY - E_sY_near;
	double D_cytB = I_cytB - E_cytB_near;

	
	
	// ATP pump that secretes cytB from cancer cells, with some stochasticity

	double cytB_binding_component;
	double cytB_recycling_component;
	double cytB_movement_component;
	double cytB_internal_binding_component;
	double cytB_internal_degradation_component;

	double ATP_binding_component; 
	double ATP_unbinding_component;
	double ATP_hydrolisis_component;
	double ADP_recycling_component;

	double pumping_probability = pCell->custom_data.variables[pumping_probability_index].value; // Fetch the pumping probability provided in the setup_tissue()

	pumping_probability = std::exp(-std::exp(I_cytB - (I_cytB/50.0) ));
	if(pumping_probability > 1.0){pumping_probability = 1.0;}

	// std::cout << "pumping probability: " << pumping_probability << std::endl;

	pCell->custom_data.variables[pumping_probability_index].value = pumping_probability;

	double pump_rate = 0.0;

	if(pCell->type == 2 && I_cytB > 0.0 ) // && pumping_probability > 0.1
	{
		// std::cout << "pumping_probability: " << pumping_probability << std::endl;

		// step 1: ATP binding to the receptor
		ATP_binding_component = diffusion_dt * ATP_binding_rate * I_ATP * Rc_cytB;
		if (ATP_binding_component > Rc_cytB * diffusion_dt || ATP_binding_component > I_ATP * diffusion_dt )
		{ ATP_binding_component = std::min(Rc_cytB * diffusion_dt, I_ATP  * diffusion_dt); } 

		I_ATP -= ATP_binding_component;
		Rcb_ATP += ATP_binding_component;
		Rc_cytB -= ATP_binding_component;

		if(I_ATP < 0.0){ I_ATP = 0.0; }
		if(Rc_cytB < 0.0){ Rc_cytB = 0.0; }

		// step 2: Some of the ATP unbinds
		ATP_unbinding_component = diffusion_dt * ATP_recycling_rate * Rcb_ATP;
		if(ATP_unbinding_component > Rcb_ATP * diffusion_dt){ ATP_unbinding_component = Rcb_ATP   * diffusion_dt; }

		I_ATP += ATP_unbinding_component;
		Rcb_ATP -= ATP_unbinding_component;
		Rc_cytB += ATP_unbinding_component;
		if(Rcb_ATP < 0.0){ Rcb_ATP = 0.0; }


		// step 3: cytB binds the ATP bomb (from inside)
		cytB_binding_component = diffusion_dt * cytB_binding_rate * I_cytB * Rc_cytB;
		if (cytB_binding_component > Rc_cytB * diffusion_dt || cytB_binding_component > I_cytB * diffusion_dt )
		{ cytB_binding_component = std::min(Rc_cytB   * diffusion_dt, I_cytB   * diffusion_dt); } // Get the limiting step 

		I_cytB -= cytB_binding_component;
		Rcb_cytB += cytB_binding_component;
		Rc_cytB -= cytB_binding_component;

		if(I_cytB < 0.0){ I_cytB = 0.0;}
		if(Rc_cytB < 0.0){ Rc_cytB = 0.0; }

		// step 4: cytB unbinds the ATP bomb (from inside)
		cytB_recycling_component = diffusion_dt * cytB_recycling_rate * Rcb_cytB;
		if (cytB_recycling_component > Rcb_cytB * diffusion_dt){ cytB_recycling_component = Rcb_cytB   * diffusion_dt; } 

		I_cytB += cytB_recycling_component;
		Rcb_cytB -= cytB_recycling_component;
		Rc_cytB += cytB_recycling_component;
		if(Rcb_cytB < 0.0){ Rcb_cytB = 0.0; }
		if(Rc_cytB < 0.0){ Rc_cytB = 0.0; }

		// step 5: ATP hydrolisis & cytB expulsion
		ATP_hydrolisis_component = diffusion_dt * ATP_hydrolisis_rate * Rcb_ATP;
		if( ATP_hydrolisis_component > Rcb_ATP * diffusion_dt ){ ATP_hydrolisis_component = Rcb_ATP * diffusion_dt; }

		cytB_movement_component = diffusion_dt * ATP_hydrolisis_rate * Rcb_cytB;
		if(cytB_movement_component > Rcb_cytB * diffusion_dt){ cytB_movement_component = Rcb_cytB * diffusion_dt; }
		
		I_ATP -= ATP_hydrolisis_component;
		I_ADP += ATP_hydrolisis_component;
		I_cytB -= cytB_movement_component;
		E_cytB_near += cytB_movement_component;

		
		if(I_ATP < 0.0){ I_ATP = 0.0; }
		if(I_cytB < 0.0){ I_cytB = 0.0; }


		Rcb_cytB -= cytB_movement_component;
		Rcb_ATP -= ATP_hydrolisis_component;
		
		if(Rcb_cytB < 0.0){ Rcb_cytB = 0.0; }
		if(Rcb_ATP < 0.0){ Rcb_ATP = 0.0; }

		Rc_cytB += cytB_movement_component;
		Rc_cytB += ATP_hydrolisis_component;
		

		// step 6: ADP recycling
		ADP_recycling_component = diffusion_dt * ADP_recycling_rate * I_ADP;
		if (ADP_recycling_component > I_ADP  ){ ADP_recycling_component = I_ADP  ; }

		I_ADP -= ADP_recycling_component;
		if(I_ADP < 0.0){ I_ADP = 0.0; }
		I_ATP += ADP_recycling_component;

		double cytB_v_max = cytB_movement_rate * (Rc_cytB + Rcb_cytB + Rcb_ATP); 
		double cytB_Km = (cytB_recycling_rate + ATP_hydrolisis_component) / cytB_binding_rate; // k_n1 + k_2 / k_1
		double cytB_monod = cytB_v_max * (I_cytB / (cytB_Km + I_cytB));
		cytB_monod *= diffusion_dt;

		pump_rate = (cytB_v_max * I_cytB) / ( cytB_Km + I_cytB); // Option A
		pump_rate *= V_cell;

		
	}




	// Simple diffusion of cytB
	double k_cytB = 0.0; // Initialize at 0.0, change depending on user parameters at specific time after

	if (PhysiCell_globals.current_time > parameters.ints("time_add_cytB") )
	{
		// k_cytB = NormalRandom(parameters.doubles("k_cytB"), 1);
		k_cytB = parameters.doubles("k_cytB");
	}

	double flux_cytB  = (A_cell) * k_cytB * (D_cytB); // (amol/min)
	double adjusted_flux_cytB = flux_cytB;

	if ( flux_cytB < 0.0 ) { // Uptake
		adjusted_flux_cytB = - std::min( std::abs(flux_cytB), std::abs(E_cytB_near * V_cell) );
	}
	else if ( flux_cytB > 0.0 ) { // Secretion
		adjusted_flux_cytB = std::min(flux_cytB, std::abs(I_cytB * V_cell ) );
	}

	pCell->custom_data.variables[cytB_flux_index].value = flux_cytB * diffusion_dt; // this is what happens in the solver, actually
	pCell->custom_data.variables[adjusted_cytB_flux_index].value = adjusted_flux_cytB * diffusion_dt;

	pCell->phenotype.secretion.net_export_rates[cytB_index] = adjusted_flux_cytB + pump_rate ; 


	
	// Separated components - FDC

	double binding_component;
	double recycling_component;
	double endocytosis_component;
	double internal_binding_component;
	double internal_degradation_component;

	// step 0: Modification on the amount of GLUT receptors based on the amount of cytochalasin B

	if (pCell -> type == 2)
	{ 
		Rc_sY += ( (adjusted_flux_cytB + pump_rate / V_cell)) * diffusion_dt; 
		Rcb_sY -= ( (adjusted_flux_cytB + pump_rate / V_cell)) * diffusion_dt;

		if(Rc_sY > parameters.doubles("Initial_Rc_sY") )
		{ 
			Rc_sY = parameters.doubles("Initial_Rc_sY");
			Rcb_sY = 0.0;
		}
		if(Rcb_sY > parameters.doubles("Initial_Rc_sY") )
		{ 
			Rcb_sY = parameters.doubles("Initial_Rc_sY");
			Rc_sY = 0.0;
		
		}
		
		if(Rc_sY < 0.0){ Rc_sY = 0.0; }
		if(Rcb_sY < 0.0){ Rcb_sY = 0.0; }
	
	}

	// if (pCell -> type == 1){ Rc_sY += ( (adjusted_flux_cytB /V_cell)) * diffusion_dt; }

	 // /V_cell makes the process really slow, but maybe more realistic 
	// Rcb_sY -= ( (adjusted_flux_cytB + pump_rate/V_cell) * diffusion_dt);
	if(Rc_sY < 0.0) { Rc_sY = 0.0; }
	if(Rcb_sY < 0.0) { Rcb_sY = 0.0; }

	binding_component = diffusion_dt * sY_binding_rate * Rc_sY * E_sY_near;

	if (binding_component > Rc_sY   || binding_component > E_sY_near   )
	{ binding_component = std::min(Rc_sY * V_cell  , E_sY_near * V_voxel  ); } // Get the limiting step 
	
	Rcb_sY += binding_component;
	Rc_sY -= binding_component;	


	recycling_component = diffusion_dt * sY_recycling_rate * Rcb_sY;
	if (recycling_component > Rcb_sY  ){ recycling_component = Rcb_sY * V_cell  ; } 
	Rc_sY += recycling_component;
	Rcb_sY -= recycling_component;
	if(Rcb_sY < 0.0) { Rcb_sY = 0.0; }


	endocytosis_component = diffusion_dt * sY_endocytosis_rate * Rcb_sY;
	if (endocytosis_component > Rcb_sY  ){ endocytosis_component = Rcb_sY * V_cell  ; }
	I_sY += endocytosis_component;
	Rcb_sY -= endocytosis_component;
	if(Rcb_sY < 0.0) { Rcb_sY = 0.0; }
	Rc_sY += endocytosis_component;


	internal_binding_component = diffusion_dt * sY_internal_binding_rate * Rc_sY * I_sY;
	if (internal_binding_component > Rcb_sY   || internal_binding_component > I_sY   )
	{ internal_binding_component = std::min(Rcb_sY * V_cell  , I_sY * V_cell  ); } // Get the limiting step 
	
	I_sY -= internal_binding_component;
	if(I_sY < 0.0) { I_sY = 0.0; }
	Rcb_sY += internal_binding_component;
	Rc_sY -= internal_binding_component;


	double vf_max = sY_endocytosis_rate * (Rc_sY + Rcb_sY); 
	double vb_max = sY_recycling_rate * (Rc_sY + Rcb_sY);
	double Km_1 = (sY_recycling_rate + sY_endocytosis_rate) / sY_binding_rate;
	double Km_2 = (sY_recycling_rate + sY_endocytosis_rate) / sY_internal_binding_rate;

	double sY_monod = vf_max * (I_sY / (Km_1 + I_sY)); // Using the "forward" (entry) kinetic values
	sY_monod *= diffusion_dt;


	double v_formation = - ((vf_max * E_sY_near) / Km_1 - (vb_max * I_sY) / Km_2) / (1 + (I_sY / Km_2) + (E_sY_near / Km_1)); // Option A
	v_formation = v_formation * V_cell; // Convert mM/min to amol/min

	
	// Internal degradation of sY - Set to kill HC in 1000 min if Initial_I_sY is 0.3 mM

	// double I_sY_degradation;
	// if(pCell->type == 1){ I_sY_degradation = parameters.doubles("k_3") * I_sY * diffusion_dt; }
	// if(pCell->type == 2){ I_sY_degradation = parameters.doubles("k_3") * I_sY * diffusion_dt; } // CC consume glucose 5 times faster - Similar to a Warburg effect
	
	// I_sY -= I_sY_degradation + (sY_monod * V_cell); // Base degration rate + Consumption associated to growth
	// if(I_sY < 1e-10){ I_sY = 0.0; }


	// Substrate concentrations

	// For sY
	pCell->custom_data.variables[E_sY_near_index].value = E_sY_near; // Mean of near_density_vector
	pCell->custom_data.variables[E_sY_index].value = E_sY_mean;		 // Mean of all voxels
	pCell->custom_data.variables[I_sY_index].value = I_sY;			 // Density

	pCell->custom_data.variables[DC_sY_index].value = DC_sY; // Substrate diffusion coefficient
	pCell->custom_data.variables[D_sY_index].value = D_sY;	 // Substrate diffusion coefficient

	// For cytB
	if(E_cytB_near < 0.0){ E_cytB_near = 0.0; }
	pCell->custom_data.variables[E_cytB_near_index].value = E_cytB_near; // Mean of near_density_vector
	pCell->custom_data.variables[E_cytB_index].value = E_cytB_mean;		 // Mean of all voxels
	pCell->custom_data.variables[I_cytB_index].value = I_cytB;			 // Density

	pCell->custom_data.variables[DC_cytB_index].value = DC_cytB; // Substrate diffusion coefficient
	pCell->custom_data.variables[D_cytB_index].value = D_cytB;	 // Substrate diffusion coefficient

	
	// Fluxes

	pCell->custom_data.variables[sY_flux_index].value = v_formation;
	pCell->custom_data.variables[pump_rate_index].value = pump_rate;

	pCell->phenotype.secretion.net_export_rates[sY_index] = v_formation / diffusion_dt; // Avoid second multiplication by diffusion_dt

	// Receptor concentrations
	pCell->custom_data.variables[Rc_sY_index].value = Rc_sY;
	pCell->custom_data.variables[Rcb_sY_index].value = Rcb_sY;
	pCell->custom_data.variables[total_Rc_sY_index].value = Rcb_sY + Rc_sY;

	pCell->custom_data.variables[Rcb_ATP_index].value = Rcb_ATP;
	pCell->custom_data.variables[Rcb_cytB_index].value = Rcb_cytB;
	pCell->custom_data.variables[total_Rcb_cytB_index].value = Rcb_ATP + Rcb_cytB;
	pCell->custom_data.variables[Rc_cytB_index].value = Rc_cytB;
	pCell->custom_data.variables[total_Rc_cytB_index].value = Rc_cytB + Rcb_cytB + Rcb_ATP;



	// Net total amounts
	pCell->custom_data.variables[total_E_sY_index].value = total_E_sY;
	pCell->custom_data.variables[total_I_sY_index].value = total_I_sY;
	pCell->custom_data.variables[total_sY_index].value = total_E_sY + total_I_sY;

	pCell->custom_data.variables[total_E_cytB_index].value = total_E_cytB;
	pCell->custom_data.variables[total_I_cytB_index].value = total_I_cytB;
	pCell->custom_data.variables[total_cytB_index].value = total_E_cytB + total_I_cytB;

	// Kinetics
	pCell->custom_data.variables[binding_component_index].value = binding_component;
	pCell->custom_data.variables[endocytosis_component_index].value = endocytosis_component;
	pCell->custom_data.variables[recycling_component_index].value = recycling_component;
	pCell->custom_data.variables[internal_binding_component_index].value = internal_binding_component;

	pCell->custom_data.variables[cytB_binding_component_index].value = cytB_binding_component;
	pCell->custom_data.variables[cytB_recycling_component_index].value = cytB_recycling_component;
	pCell->custom_data.variables[cytB_movement_component_index].value = cytB_movement_component;
	pCell->custom_data.variables[cytB_internal_binding_component_index].value = cytB_internal_binding_component;

	pCell->custom_data.variables[ATP_binding_component_index].value = ATP_binding_component;
	pCell->custom_data.variables[ATP_unbinding_component_index].value = ATP_unbinding_component;
	pCell->custom_data.variables[ATP_hydrolisis_component_index].value = ATP_hydrolisis_component;
	pCell->custom_data.variables[ADP_recycling_component_index].value = ADP_recycling_component;

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

	int glucose_flux_threshold = 0.0;

	// Connecting the transport model to the phenotype

	// Monod function that modifies the rate of duplication

	// if (I_sY < parameters.doubles("cytotoxic_threshold") && parameters.doubles("cytotoxic_threshold") != 0.0)
	// {
	// 	std::cout << "Not enough glucose... " << std::endl;
	// 	pCell->start_death( necrosis_idx ); // "Custom" apoptosis 
	// }

	if (pCell->type == 2) { pCell->phenotype.cycle.data.transition_rate(0,0) = sY_monod * V_cell; }
	// if (pCell->type == 1) { pCell->phenotype.cycle.data.transition_rate(0,0) = 0.0; } // HC don't divide
	// if ( I_sY < 1e-06 ){ pCell->start_death( necrosis_idx ); std::cout << "Necrosis" << std::endl; }
	
	// if (I_sY < parameters.doubles("growth_threshold"))
	// {
	// 	// pCell->flag_for_division();
	// 	// pCell->phenotype.cycle.advance_cycle(pCell, phenotype, diffusion_dt);
	// 	pCell->phenotype.cycle.data.transition_rate(0,0) = 1e-08;
	// }

	// if (I_sY < 0.01)
	// {
	// 	glucose_flux_threshold += 1;
	// 	if (glucose_flux_threshold > 3)
	// 	{
	// 		std::cout << "Not enough glucose flux..." << std::endl;
	// 		pCell->start_death( necrosis_idx );
	// 	}

	// }


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

	Cell_Definition* pCD = find_cell_definition( pCell->type_name );
	// int internal_sY = (int)round(pCell->custom_data[1] * 200);

	static int sY_index = microenvironment.find_density_index("sY");
	static int cytB_index = microenvironment.find_density_index("cytB");
	double I_sY = pCell->phenotype.molecular.internalized_total_substrates[sY_index] * 300;
	double I_cytB = pCell->phenotype.molecular.internalized_total_substrates[cytB_index] * 300;

	int internal_Rc = (int)round(pCell->custom_data[7] * 200); // GLUT receptors
	int internal_cytB = (int)round(pCell->custom_data[4] * 200); 


	char szTempString[128];

	std::vector<std::string> output = false_cell_coloring_live_dead(pCell);

	if (pCell->type == 1) // Healthy - Blue
	{

		output[1] = "blue"; // outer cytoplasm
		
		sprintf(szTempString, "rgb(%u,%u,%u)", 173, 173, 173 + I_sY);
		output[0].assign(szTempString); // cytoplasm

		sprintf(szTempString, "rgb(%u,%u,%u)", 0, internal_Rc, 0);
		output[2] = "rgb(230,230,230)"; // nucleus
		output[3].assign(szTempString); // outer nucleus

	}

	if(pCell->type == 2) // Cancer - Red
	{
		output[1]= "red";

		sprintf(szTempString, "rgb(%u,%u,%u)", 173 + I_sY, 173, 173);
		output[0].assign(szTempString);

		sprintf(szTempString, "rgb(%u,%u,%u)", 0, internal_Rc, 0);
		output[3].assign(szTempString); // outer nucleus

		sprintf(szTempString, "rgb(%u,%u,%u)", I_cytB ,I_cytB, I_cytB);
		output[2] = "rgb(230,230,230)"; // nucleus


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



	/////////// ANOTHER OPTION

	// OPTION B: Random disposition of HC and CC

	// for( int n = 0 ; n < parameters.ints("number_healthy") ; n++ )
	// {
	// 	std::cout << "Creating healthy cell " << n << std::endl;
	// 	std::vector<double> position = {0.0, 0.0, 0.0};
	// 	position[0] = Xmin + UniformRandom()*Xrange;
	// 	position[1] = Ymin + UniformRandom()*Yrange;
	// 	position[2] = Zmin + UniformRandom()*Zrange;
	// 	pCell = create_cell( *pCH );
	// 	// std::cout << "this is the position of the cell: " << position[0] << " " << position[1] << " " << position[2] << std::endl;
	// 	pCell->assign_position( position );
      
	// 	pCell->phenotype.molecular.internalized_total_substrates[sY_index] = 0.0;

	// 	pCell->custom_data[7] = NormalRandom(R_mean, R_sd);
	// 	if (pCell->custom_data[7] < R_min) { pCell->custom_data[7] = R_min; }
	// 	if (pCell->custom_data[7] > R_max) { pCell->custom_data[7] = R_max; }

	// 	// std::cout << "And its I_sY is: " << pCell->custom_data[7] << std::endl;
	
	// }

	// for( int n = 0 ; n < parameters.ints("number_cancer") ; n++ )
	// {
	// 	// std::cout << "Creating cancer cell " << n << std::endl;
	// 	std::vector<double> position = {0.0, 0.0, 0.0};
	// 	position[0] = Xmin + UniformRandom()*Xrange;
	// 	position[1] = Ymin + UniformRandom()*Yrange;
	// 	position[2] = Zmin + UniformRandom()*Zrange;
	// 	pCell = create_cell( *pCC );
	// 	// std::cout << "this is the position of the cell: " << position[0] << " " << position[1] << " " << position[2] << std::endl;
	// 	pCell->assign_position( position );

	// 	pCell->custom_data[7] = NormalRandom(R_mean, R_sd) * 3.0;
	// 	if (pCell->custom_data[7] < R_min) { pCell->custom_data[7] = R_min; }
	// 	if (pCell->custom_data[7] > R_max) { pCell->custom_data[7] = R_max; }

	// 	// std::cout << "And its I_sY is: " << pCell->custom_data[7] << std::endl;

	// }

	// #pragma omp parallel for
	// for( int k=1; k < cell_definitions_by_index.size() ; k++ ) // k starts at 1 to avoid "default"-type cells
	// {
	// 	Cell_Definition* pCD = cell_definitions_by_index[k];
	// 	std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
	// 	for( int n = 0 ; n < parameters.ints("cells_number") ; n++ )
	// 	{
	// 		std::vector<double> position = {0,0,0};
	// 		position[0] = Xmin + UniformRandom()*Xrange;
	// 		position[1] = Ymin + UniformRandom()*Yrange;
	// 		position[2] = Zmin + UniformRandom()*Zrange;
	// 		pCell = create_cell( *pCD );
	// 		pCell->assign_position( position );
	// 	}
	// }


	// for(int k=0; k<parameters.ints("number_healthy"); k++)
	// {
	// 	std::cout << "Placing healthy cell.. " << k << std::endl;

	// 	std::vector<double> position = {0.0, 0.0, 0.0};
	// 	double r = NormalRandom(0, 1) * tumor_radius;
	// 	double theta = UniformRandom() * M_PI * 2;

	// 	position[0] = center_x + r*std::cos(theta);
	// 	position[1] = center_y + r*std::sin(theta);
	// 	position[2] = center_z;

	// 	pCell = create_cell( *pCH );
	// 	pCell->assign_position( position );
      
	// 	pCell->phenotype.molecular.internalized_total_substrates[sY_index] = 0.0;

	// 	pCell->custom_data[7] = NormalRandom(R_mean, R_sd);
	// 	if (pCell->custom_data[7] < R_min) { pCell->custom_data[7] = R_min; }
	// 	if (pCell->custom_data[7] > R_max) { pCell->custom_data[7] = R_max; }
	// }

	

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


// Funtion to read init files created with PhysiBoSSv2
//unneaded
std::vector<init_record> read_init_file(std::string filename, char delimiter, bool header)
{
	// File pointer
	std::fstream fin;
	std::vector<init_record> result;

	// Open an existing file
	fin.open(filename, std::ios::in);

	// Read the Data from the file
	// as String Vector
	std::vector<std::string> row;
	std::string line, word;

	if (header)
		getline(fin, line);

	do
	{
		row.clear();

		// read an entire row and
		// store it in a string variable 'line'
		getline(fin, line);

		// used for breaking words
		std::stringstream s(line);

		// read every column data of a row and
		// store it in a string variable, 'word'
		while (getline(s, word, delimiter))
		{

			// add all the column data
			// of a row to a vector
			row.push_back(word);
		}

		init_record record;
		record.x = std::stof(row[2]);
		record.y = std::stof(row[3]);
		record.z = std::stof(row[4]);
		record.radius = std::stof(row[5]);
		// record.phase = std::stoi(row[13]);
		// record.elapsed_time = std::stod(row[14]);

		result.push_back(record);
	} while (!fin.eof());

	return result;
}


void my_mutation_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Fetch from custom data instead

	double R_mean = parameters.doubles("R_mean");
	double R_sd = parameters.doubles("R_sd");
	double R_min = parameters.doubles("R_min");
	double R_max = parameters.doubles("R_max");

	pCell->custom_data[7] = NormalRandom(parameters.doubles("Initial_Rc_sY"), R_sd);
	if (pCell->custom_data[7] < R_min) { pCell->custom_data[7] = R_min; }
	if (pCell->custom_data[7] > R_max) { pCell->custom_data[7] = R_max; }

	return;
}

bool volume_arrest_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	if( phenotype.volume.total < 2494 )
	{ return true; }
	return false;
}

void inject_cytB(int density_index, double concentration) 
{
	std::cout << "Adding cytochalasin b... " << std::endl;

	#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		auto current_voxel = microenvironment.voxels(n);
		// std::vector<double> cent = {current_voxel.center[0], current_voxel.center[1], current_voxel.center[2]};
		microenvironment.density_vector(n)[density_index] = concentration;
	}
}


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

	// if not, dead colors

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