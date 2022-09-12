
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