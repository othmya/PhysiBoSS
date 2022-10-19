/*
 * ags_boolean_model_interface.cpp
 *
 */


#include "./boolean_model_interface.h"
#include <math.h>

using namespace PhysiCell; 

Submodel_Information bm_interface_info;

double Hill_response_function( double s, double half_max , double hill_power )
{ 
    // newer. only one expensive a^b operation. 45% less computationl expense. 

	// give an early exit possibility to cut cost on "empty" rules
	if(s < 1e-16 ) // maybe also try a dynamic threshold: 0.0001 * half_max 
	{ return 0.0; } 

	// operations to reduce a^b operations and minimize hidden memory allocation / deallocation / copy operations. 
	// Hill = (s/half_max)^hill_power / ( 1 + (s/half_max)^hill_power  )
	double temp = s; // s 
	temp /= half_max; // s/half_max 
	double temp1 = pow(temp,hill_power); // (s/half_max)^h 
	temp = temp1;  // (s/half_max)^h 
	temp +=1 ;  // (1+(s/half_max)^h ); 
	temp1 /= temp; // (s/half_max)^h / ( 1 + s/half_max)^h) 
    if(temp1 < 1e-16)
        temp1 = 0.0;
	return temp1; 
}


void boolean_model_interface_setup()
{
    bm_interface_info.name = "AGS Boolean model interface"; 
	bm_interface_info.version = "0.0.1";
	
    bm_interface_info.main_function = ags_bm_interface_main; 

	// These are just auxiliary variables to keep track of some BN nodes

    bm_interface_info.cell_variables.push_back( "mek_node" );
    bm_interface_info.cell_variables.push_back( "pi3k_node" );
    bm_interface_info.cell_variables.push_back( "tak1_node" );
    bm_interface_info.cell_variables.push_back( "akt_node" );

    bm_interface_info.cell_variables.push_back( "anti_mek_node" );
    bm_interface_info.cell_variables.push_back( "anti_pi3k_node" );
    bm_interface_info.cell_variables.push_back( "anti_tak1_node" );
    bm_interface_info.cell_variables.push_back( "anti_akt_node" );


    bm_interface_info.cell_variables.push_back( "prosurvival_b1_node" );
    bm_interface_info.cell_variables.push_back( "prosurvival_b2_node" );
    bm_interface_info.cell_variables.push_back( "prosurvival_b3_node" );

    bm_interface_info.cell_variables.push_back( "antisurvival_b1_node" );
    bm_interface_info.cell_variables.push_back( "antisurvival_b2_node" );
    bm_interface_info.cell_variables.push_back( "antisurvival_b3_node" );


    // Could add here output of transfer functions

	bm_interface_info.register_model();
}

void update_boolean_model_inputs( Cell* pCell, Phenotype& phenotype, double dt )
{

    if( pCell->phenotype.death.dead == true )
	{ return; }

    anti_node_mapping_function(pCell, "drug_X", parameters.strings("drug_X_target"), parameters.doubles("drug_X_half_max"), parameters.doubles("drug_X_Hill_coeff"));
    anti_node_mapping_function(pCell, "drug_Y", parameters.strings("drug_Y_target"), parameters.doubles("drug_Y_half_max"), parameters.doubles("drug_Y_Hill_coeff"));

    return;
}


void update_cell_from_boolean_model(Cell* pCell, Phenotype& phenotype, double dt)
{	
    static int drug_X_ix = microenvironment.find_density_index( "drug_X" );
	static int drug_Y_ix = microenvironment.find_density_index( "drug_Y" );

    // static int drug_X_export_rate = pCell->custom_data.find_variable_index( "drug_X_net_production_rate" );

    static int apoptosis_model_index = phenotype.death.find_death_model_index( "Apoptosis" );
    static int necrosis_model_index = phenotype.death.find_death_model_index( "Necrosis" );
    
    // Getting the state of the boolean model readouts (Readout can be in the XML)

    bool casp37_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "Caspase37_b1" );
    bool casp37_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "Caspase37_b2" );
    bool FOXO = pCell->phenotype.intracellular->get_boolean_variable_value( "FOXO" );

    bool antisurvival_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b1" );
    bool antisurvival_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b2" );
    bool antisurvival_b3 = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b3" );

    double anti_w1 = 0.1; // Introduce them from within the XML file
    double anti_w2 = 0.2;
    double anti_w3 = 0.7;
    double S_anti = (anti_w1*antisurvival_b1) + (anti_w2 * antisurvival_b2) + (anti_w3 * antisurvival_b3);
    double S_anti_real = (anti_w1*casp37_b1) + (anti_w2 * casp37_b2) + (anti_w3 * FOXO);

    bool prosurvival_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b1" );
    bool prosurvival_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b2" );
    bool prosurvival_b3 = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b3" );

    bool cMYC = pCell->phenotype.intracellular->get_boolean_variable_value( "cMYC" );
    bool CCND_b1 = pCell->phenotype.intracellular->get_boolean_variable_value( "CCND1_b1" );
    bool CCND_b2 = pCell->phenotype.intracellular->get_boolean_variable_value( "CCND1_b2" );

    double pro_w1 = 0.2;
    double pro_w2 = 0.3;
    double pro_w3 = 0.5;
    double S_pro = (pro_w1*prosurvival_b1) + (pro_w2 * prosurvival_b2) + (pro_w3 * prosurvival_b3);
    double S_pro_real = (pro_w1*cMYC) + (pro_w2 * CCND_b1) + (pro_w3 * CCND_b2);

    // Connect output from model to actual cell variables

    double apoptosis_rate_basal = parameters.doubles("apoptosis_rate_basal");
    double maximum_apoptosis_rate = parameters.doubles("max_apoptosis_rate");
    double hill_coeff_apoptosis = parameters.doubles("hill_coeff_apoptosis");
    double K_half_apoptosis = parameters.doubles("K_half_apoptosis");

    double basal_growth_rate = parameters.doubles("basal_growth_rate");
    double hill_coeff_growth = parameters.doubles("hill_coeff_growth");
    double K_half_growth = parameters.doubles("K_half_growth");
    double growth_value = growth_mapping_logistic(basal_growth_rate, hill_coeff_growth, K_half_growth, S_pro);
    double growth_value_Hill = Hill_response_function(S_pro, K_half_growth, hill_coeff_growth); // Max value is 1 

    // Apoptosis mapping
    pCell-> phenotype.death.rates[apoptosis_model_index] =
         apoptosis_mapping_logistic(apoptosis_rate_basal, maximum_apoptosis_rate, hill_coeff_apoptosis, K_half_apoptosis, S_anti);
    
    // Growth mapping
    phenotype.cycle.data.transition_rate(0, 0) = growth_value_Hill * basal_growth_rate;

    
    return;
}


void update_monitor_variables(Cell* pCell )
{
	static int mek_node_ix = pCell->custom_data.find_variable_index("mek_node");
	static int akt_node_ix = pCell->custom_data.find_variable_index("akt_node");
	static int pi3k_node_ix = pCell->custom_data.find_variable_index("pi3k_node");
	static int tak1_node_ix = pCell->custom_data.find_variable_index("tak1_node");

    static int anti_mek_node_ix = pCell->custom_data.find_variable_index("anti_mek_node");
    static int anti_akt_node_ix = pCell->custom_data.find_variable_index("anti_akt_node");
    static int anti_pi3k_node_ix = pCell->custom_data.find_variable_index("anti_pi3k_node");
    static int anti_tak1_node_ix = pCell->custom_data.find_variable_index("anti_tak1_node");

    static int antisurvival_b1_ix = pCell->custom_data.find_variable_index("antisurvival_b1_node");
    static int antisurvival_b2_ix = pCell->custom_data.find_variable_index("antisurvival_b2_node");
    static int antisurvival_b3_ix = pCell->custom_data.find_variable_index("antisurvival_b3_node");

    static int prosurvival_b1_ix = pCell->custom_data.find_variable_index("prosurvival_b1_node");
    static int prosurvival_b2_ix = pCell->custom_data.find_variable_index("prosurvival_b2_node");
    static int prosurvival_b3_ix = pCell->custom_data.find_variable_index("prosurvival_b3_node");

	pCell->custom_data[mek_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "MEK" );
    pCell->custom_data[akt_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "AKT" );
    pCell->custom_data[pi3k_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "PI3K" );
    pCell->custom_data[tak1_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "TAK1" );

    pCell->custom_data[anti_mek_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_MEK" );
    pCell->custom_data[anti_akt_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_AKT" );
    pCell->custom_data[anti_pi3k_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_PI3K" );
    pCell->custom_data[anti_tak1_node_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "anti_TAK1" );

    pCell->custom_data[antisurvival_b1_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b1" );
    pCell->custom_data[antisurvival_b2_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b2" );
    pCell->custom_data[antisurvival_b3_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Antisurvival_b3" );

    pCell->custom_data[prosurvival_b1_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b1" );
    pCell->custom_data[prosurvival_b2_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b2" );
    pCell->custom_data[prosurvival_b3_ix] = pCell->phenotype.intracellular->get_boolean_variable_value( "Prosurvival_b3" );

    return;
}


void ags_bm_interface_main(Cell* pCell, Phenotype& phenotype, double dt)
{

    if( phenotype.death.dead == true )
	{
		pCell->functions.update_phenotype = NULL;
		return;
	}

    if ( pCell->phenotype.intracellular->need_update() )
    {

        // First we update the Boolean Model inputs
        update_boolean_model_inputs(pCell, phenotype, dt );
    
        // Run maboss to update the boolean state of the cell
        pCell->phenotype.intracellular->update();

        // update the cell fate based on the boolean outputs
        update_cell_from_boolean_model(pCell, phenotype, dt);

        // Get track of some boolean node values for debugging
        update_monitor_variables(pCell);
    }

    return;
}


bool boolean_node_deactivation_prob(double drug_density, double scaling, double GI50 ){

    // Logistic function

    double node_probability = 1 / (1 + std::exp(- log10(drug_density / GI50 ) * scaling));
    
    std::default_random_engine generator;
    std::bernoulli_distribution distro(node_probability);

    if (distro(generator)){ node_probability = 1;} else { node_probability = 0;}
    // if (node_probability <= 0.5){ node_probability = 0;} else { node_probability = 1;}

    return node_probability;

} 

double growth_mapping_logistic(double doubling_time, double hill_coeff, double K_half, double S_value){

    // double growth_logistic_function = doubling_time / (1 + std::exp(- log10(readout_value) * scaling));

    double growth_logistic_function;
    growth_logistic_function = (doubling_time * std::pow(S_value, hill_coeff ) ) / (K_half + std::pow(S_value, hill_coeff) ) ;

    return growth_logistic_function;

}

double apoptosis_mapping_logistic(double basal_apoptosis_rate, double maximum_apoptosis_rate, double hill_coeff, double K_half, double S_value){


    double apoptosis_mapping_function;
    apoptosis_mapping_function = (maximum_apoptosis_rate * std::pow(S_value, hill_coeff) ) / (K_half + std::pow(S_value, hill_coeff) ) ;

    return apoptosis_mapping_function +  basal_apoptosis_rate;
    
}

double pressure_effect_growth_rate(double pressure, double hill_coeff, double pressure_half){

    // Suggestion: Employ the Hill_effect function from PhysiCell instead of this one, just to align with the other mapping functions

    // double pressure_exponential_function = std::pow(6e-03, pressure);
    double pressure_exponential_function =  std::pow(pressure, hill_coeff) / (pressure_half + std::pow(pressure, hill_coeff));
    // if (pressure_exponential_function > 1) pressure_exponential_function = 1.0;
    return pressure_exponential_function;
}


void anti_node_mapping_function( Cell* pCell, std::string drug_name, std::string target_node, double drug_half_max, double drug_Hill_coeff){

    double cell_volume = pCell->phenotype.volume.total;
	static int drug_idx = microenvironment.find_density_index( drug_name );

    double drug_int = pCell->phenotype.molecular.internalized_total_substrates[drug_idx];
    drug_int /= cell_volume; // Convert to density (mM)

    double target_inactivate_p = Hill_response_function(drug_int, drug_half_max, drug_Hill_coeff );

    // std::cout << "For drug " << drug_name << " uniform_random value " << uniform_random() << " target inactivation prob " << target_inactivate_p << std::endl;

    if (target_node != "none"){
        if ( uniform_random() <  target_inactivate_p ){ // Added normalization by maximum rand value
            pCell->phenotype.intracellular->set_boolean_variable_value(target_node, 1);
        } else { 
            pCell->phenotype.intracellular->set_boolean_variable_value(target_node, 0);
        }
    }

    return;
}