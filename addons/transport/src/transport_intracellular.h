#ifndef _transport_Intracellular_h_
#define _transport_Intracellular_h_

#include <string>
#include <map>
#include <iomanip>   // for setw

#include "../../../core/PhysiCell.h"
#include "../../../core/PhysiCell_phenotype.h"
#include "../../../core/PhysiCell_cell.h"
#include "../../../modules/PhysiCell_pugixml.h"
// #include "maboss_network.h"

// #ifdef ADDON_ROADRUNNER
// These are for C
// #define STATIC_RRC
// #include "rrc_api.h"
// #include "rrc_types.h"
// #include "../roadrunner/include/rr/C/rrc_api.h"
// #include "../roadrunner/include/rr/C/rrc_types.h"
// #include "rrc_utilities.h"
// extern "C" rrc::RRHandle createRRInstance();
// #endif

using namespace std;


struct kinetic_parm
{
	string name;
	string untis;
	float value;
};

struct initial_conds
{
	int density_index;
	string name;
	string untis;
	float value;
};

struct transport_model 
{
	int density_index;

	initial_conds I_density_initial;
	initial_conds E_density_initial;

	kinetic_parm k1;
	kinetic_parm kn1;
	kinetic_parm k2;
	kinetic_parm kn2;

	kinetic_parm Km_f;
	kinetic_parm Km_r;
	kinetic_parm Vmax_f;
	kinetic_parm Vmax_r;

};

class transportIntracellular : public PhysiCell::Intracellular 
{
 private:
 public:
	
	static long counter;

	// No external files needed
	
	std::map<std::string, double> initial_values;
	std::map<std::string, double> parameters;
	std::map<std::string, std::string> substrate_species;

    transportIntracellular();

	transportIntracellular(pugi::xml_node& node);
	
	transportIntracellular(transportIntracellular* copy);
	
    // rwh: review this
	Intracellular* clone()
    {
		// return static_cast<Intracellular*>(new RoadRunnerIntracellular(this));
		transportIntracellular* clone = new transportIntracellular(this);
		// clone->sbml_filename = this->sbml_filename;
		clone->substrate_species = this->substrate_species;
        // clone->phenotype_species = this->phenotype_species;
		// clone->custom_data_species = this->custom_data_species;
		return static_cast<Intracellular*>(clone);
	}

	Intracellular* getIntracellularModel() 
    {
        std::cout << "------ transportIntracellular: getIntracellularModel called\n";
		return static_cast<Intracellular*>(this);
	}
	
	void initialize_intracellular_from_pugixml(pugi::xml_node& node);
	
    // Need 'int' return type to avoid bizarre compile errors? But 'void' to match MaBoSS.
	void start();

	bool need_update();

    // Need 'int' return type to avoid bizarre compile errors.
	void update();
    
    int update_phenotype_parameters(PhysiCell::Phenotype& phenotype);
    int validate_PhysiCell_tokens(PhysiCell::Phenotype& phenotype);
    // int validate_SBML_species();
    // int create_custom_data_for_SBML(PhysiCell::Phenotype& phenotype);
	
	double get_parameter_value(std::string name);
	void set_parameter_value(std::string name, double value);
	
	std::string get_state();

    // // for now, define dummy methods for these in the abstract parent class
    // bool has_variable(std::string name) { return false; }
	// bool get_boolean_variable_value(std::string name) { return false; }
	// void set_boolean_variable_value(std::string name, bool value)  {}
    // void print_current_nodes() {}

	// static void save_libRR(std::string path, std::string index);
};

#endif