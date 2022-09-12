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

#include "./custom.h"

// declare cell definitions here
void create_cell_types(void)
{
	// use the same random seed so that future experiments have the
	// same initial histogram of oncoprotein, even if threading means
	// that future division and other events are still not identical
	// for all runs

	SeedRandom(parameters.ints("random_seed")); // or specify a seed here

	// housekeeping
	std::cout << cell_defaults.name << std::endl;

	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	// no migration_bias needed
	cell_defaults.functions.update_migration_bias = NULL;
	// No custom rule needed
	cell_defaults.functions.custom_cell_rule = NULL;

	cell_defaults.custom_data.add_variable("time", "min", PhysiCell_globals.current_time );

	/*
		update_pc_parameters_O2_based flag controls how is the update_phenotypes
	 	if update_pc_parameters_O2_based is set to true	the model will call
		tumor_cell_phenotype_with_signaling which just calls two functions 
		sequentially: 1) update_cell_and_death_parameters_O2_based and 
		2) tnf_bm_interface_main; the former updates growth and death rates 
		based on oxygen while the second is the	function that update the boolean model. 
		If the flag is false then only the tnf_bm_interface_main is invoked
	*/
	if (parameters.bools("update_pc_parameters_O2_based"))
	{
		cell_defaults.functions.update_phenotype = tumor_cell_phenotype_with_signaling;
	}
	else
	{
		cell_defaults.functions.update_phenotype = ags_bm_interface_main; // ags_bm_interface_main
	}

	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
	cell_defaults.functions.calculate_distance_to_membrane = NULL;
	cell_defaults.functions.set_orientation = NULL;

	/*
	   This parses the cell definitions in the XML config file. 
	*/
	initialize_cell_definitions_from_pugixml();

	// initialize tnf
	// the following should be changed
	//commented today
	ags_receptor_model_setup();
	ags_boolean_model_interface_setup();
	submodel_registry.display(std::cout);

	// Needs to initialize one of the receptor state to the total receptor value

	// cell_defaults.custom_data["unbound_external_TNFR"] = cell_defaults.custom_data["TNFR_receptors_per_cell"];
	// cell_defaults.custom_data["bound_external_TNFR"] = 0;
	// cell_defaults.custom_data["bound_internal_TNFR"] = 0;

	build_cell_definitions_maps();
	display_cell_definitions(std::cout);

	return;
}

void setup_microenvironment(void)
{
	// make sure to override and go back to 2D
	// if (default_microenvironment_options.simulate_2D == true)
	// {
	// 	std::cout << "Warning: overriding XML config option and setting to 3D!" << std::endl;
	// 	default_microenvironment_options.simulate_2D = false;
	// }

	// initialize BioFVM
	initialize_microenvironment();

	return;
}

void setup_tissue( void )
{
	// place a cluster of tumor cells at the center 
	
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.8 * 2.0 * cell_radius; 
	
	double tumor_radius = 
		parameters.doubles("tumor_radius"); // 250.0; 
	
	Cell* pCell = NULL; 
	
	// std::vector<std::vector<double>> positions = create_cell_sphere_positions(cell_radius, tumor_radius); 
	std::vector<std::vector<double>> positions = create_cell_circle_positions(cell_radius, tumor_radius); 
	std::cout << "creating " << positions.size() << " closely-packed tumor cells ... " << std::endl; 

	// #pragma omp parallel for 
	for( int i=0; i < positions.size(); i++ )
	{
		
		pCell = create_cell(); 
		pCell->assign_position( positions[i] );
		// std::cout << "Creating cell in position " << positions[i] << std::endl;

		pCell->phenotype.intracellular->start();

		const int max_iter = 1;

		#pragma omp parallel for 
		for (int j=0; j < max_iter; j++){
			pCell->phenotype.intracellular->update(); }
		
		// pCell->phenotype.intracellular->print_current_nodes();
		
	
		update_monitor_variables(pCell);
	
	}

	return; 
}



// void setup_tissue(void)
// {

// 	std::vector<init_record> cells = read_init_file(parameters.strings("init_cells_filename"), ';', true);

// 	for (int i = 0; i < cells.size(); i++)
// 	{
// 		Cell *pCell;
// 		float x = cells[i].x;
// 		float y = cells[i].y;
// 		float z = cells[i].z;
// 		double elapsed_time = cells[i].elapsed_time;
		
// 		pCell = create_cell(get_cell_definition("default"));
// 		pCell->phenotype.cycle.data.elapsed_time_in_phase = elapsed_time;
// 		pCell->assign_position(x, y, z);

// 		update_monitor_variables(pCell);
// 	}

// 	return;
// }

// custom cell phenotype function to run PhysiBoSS when is needed
void tumor_cell_phenotype_with_signaling(Cell *pCell, Phenotype &phenotype, double dt)
{
	if (phenotype.death.dead == true)
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

	update_cell_and_death_parameters_O2_based(pCell, phenotype, dt);
	ags_bm_interface_main(pCell, phenotype, dt);

}

//different in physiboss_cell_lines
// cell coloring function for ploting the svg files
std::vector<std::string> my_coloring_function(Cell *pCell)
{
	// start with live coloring
	std::vector<std::string> output = false_cell_coloring_live_dead(pCell);

	// dead cells
	if (pCell->phenotype.death.dead == false)
	{
		static double V_cell = pCell->phenotype.volume.total; 
		static int dA_index = microenvironment.find_density_index("dA");
		static int dB_index = microenvironment.find_density_index("dB");

		double I_dA = pCell->phenotype.molecular.internalized_total_substrates[dA_index] / V_cell;
		double I_dB = pCell->phenotype.molecular.internalized_total_substrates[dB_index] / V_cell;

		float activation_threshold = pCell->custom_data.find_variable_index("activation threshold");

		// int Rcb_dA = (int)round((pCell->custom_data[Rcb_dA_ix] / activation_threshold) * 255.0);

		if (I_dA > 0 || I_dB > 0)
		{
			char szTempString[128];
			sprintf(szTempString, "rgb(%u,%u,%u)", 255, 255 - I_dB, 255 - I_dA);
			output[0].assign("black");
			output[1].assign(szTempString);
			output[2].assign("black");
			output[3].assign(szTempString);
		}
	}

	return output;
}

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
		record.radius = std::stof(row[5]);
		record.phase = std::stoi(row[13]);
		record.z = std::stof(row[4]);
		record.elapsed_time = std::stod(row[14]);

		result.push_back(record);
	} while (!fin.eof());

	return result;
}

// aded from example 
void set_input_nodes(Cell* pCell) {}

void from_nodes_to_cell(Cell* pCell, Phenotype& phenotype, double dt) {}


void color_node(Cell* pCell){
	std::string node_name = parameters.strings("node_to_visualize");
	pCell->custom_data[node_name] = pCell->phenotype.intracellular->get_boolean_variable_value(node_name);
}
void inject_density_sphere(int density_index, double concentration, double membrane_lenght)
{
// Inject given concentration on the extremities only
#pragma omp parallel for
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
	{
		auto current_voxel = microenvironment.voxels(n);
		std::vector<double> cent = {current_voxel.center[0], current_voxel.center[1], current_voxel.center[2]};

		if ((membrane_lenght - norm(cent)) <= 0)
			microenvironment.density_vector(n)[density_index] = concentration;
	}
}

void remove_density(int density_index)
{
	for (int n = 0; n < microenvironment.number_of_voxels(); n++)
		microenvironment.density_vector(n)[density_index] = 0;
	// std::cout << "Removal done" << std::endl;
}

// std::string cells_message_builder(std::vector<Cell*> all_cells, double timepoint){
// 	// implement count different situation of cells
// 	int alive_no,apoptotic_no,necrotic_no=0;
// 	// int p_id = all_cells[i]->ID
// 	static int TUMOR_TYPE=0;
// 	std::string message;
// 	for(int i=0;i<all_cells.size();i++)
// 		{
// 			if( all_cells[i]->phenotype.cycle.pCycle_Model )
// 			{
// 				int code= all_cells[i]->phenotype.cycle.current_phase().code;
// 				if (code ==PhysiCell_constants::Ki67_positive_premitotic || code==PhysiCell_constants::Ki67_positive_postmitotic || code==PhysiCell_constants::Ki67_positive || code==PhysiCell_constants::Ki67_negative || code==PhysiCell_constants::live)
// 					// _nameCore="LIVE";
// 					alive_no +=1;
// 				else if (code==PhysiCell_constants::apoptotic)
// 					// _nameCore="APOP";
// 					apoptotic_no +=1;
// 				else if (code==PhysiCell_constants::necrotic_swelling || code==PhysiCell_constants::necrotic_lysed || code==PhysiCell_constants::necrotic)
// 					// _nameCore="NEC";
// 					necrotic_no += 1;
// 				// else if (code==PhysiCell_constants::debris)
// 				// 	_nameCore="DEBR";
// 				// else
// 				// 	_nameCore="MISC";
// 			}
// 			else if(all_cells[i]->type==TUMOR_TYPE)
// 				// _nameCore="LIVE";
// 				alive_no+=1;
// 		}
// 	// message = std::to_string(p_id) + ';' + std::to_string(timepoint) + ';' + std::to_string(alive_no) + ';' + std::to_string(apoptotic_no) + ';' + std::to_string(necrotic_no) + ';';
// 	message = std::to_string(timepoint) + ';' + std::to_string(alive_no) + ';' + std::to_string(apoptotic_no) + ';' + std::to_string(necrotic_no) + ';';
// 	return message;
// }

double total_live_cell_count()
{
	double out = 0.0;

	for (int i = 0; i < (*all_cells).size(); i++)
	{
		if ((*all_cells)[i]->phenotype.death.dead == false && (*all_cells)[i]->type == 0)
		{
			out += 1.0;
		}
	}

	return out;
}

double total_dead_cell_count()
{
	double out = 0.0;

	for (int i = 0; i < (*all_cells).size(); i++)
	{
		if ((*all_cells)[i]->phenotype.death.dead == true && (*all_cells)[i]->phenotype.death.current_death_model_index == 0)
		{
			out += 1.0;
		}
	}

	return out;
}

double total_necrosis_cell_count()
{
	double out = 0.0;

	for (int i = 0; i < (*all_cells).size(); i++)
	{
		if ((*all_cells)[i]->phenotype.death.dead == true && (*all_cells)[i]->phenotype.death.current_death_model_index == 1)
		{
			out += 1.0;
		}
	}

	return out;
}



std::vector<std::vector<double>> create_cell_circle_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*sqrt(3);

	std::vector<double> tempPoint(3,0.0);
	
	for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
	{
		for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
		{
			tempPoint[1]=y + (xc%2) * cell_radius;
			tempPoint[0]=x;
			tempPoint[2]=0;
			if(sqrt(norm_squared(tempPoint))< sphere_radius)
			{ cells.push_back(tempPoint); }
		}
	}
	return cells;
}

std::vector<std::vector<double>> create_cell_sphere_positions(double cell_radius, double sphere_radius)
{
	std::vector<std::vector<double>> cells;
	int xc=0,yc=0,zc=0;
	double x_spacing= cell_radius*sqrt(3);
	double y_spacing= cell_radius*2;
	double z_spacing= cell_radius*sqrt(3);

	std::vector<double> tempPoint(3,0.0);
	// std::vector<double> cylinder_center(3,0.0);

	for(double z=-sphere_radius;z<sphere_radius;z+=z_spacing, zc++)
	{
		for(double x=-sphere_radius;x<sphere_radius;x+=x_spacing, xc++)
		{
			for(double y=-sphere_radius;y<sphere_radius;y+=y_spacing, yc++)
			{
				tempPoint[0]=x + (zc%2) * 0.5 * cell_radius;
				tempPoint[1]=y + (xc%2) * cell_radius;
				tempPoint[2]=z;

				if(sqrt(norm_squared(tempPoint))< sphere_radius)
				{ cells.push_back(tempPoint); }
			}
		}
	}
	return cells;

}