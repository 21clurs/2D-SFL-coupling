// much of this adapted with little modification from SimOptions.h in Fang Da's implementation of "Surface-Only Liquids"
#include <fstream>
#include <sstream>

#include "simoptions.h"

std::map<std::string, SimOptions::SimOption> SimOptions::sim_options;

void SimOptions::addStringOption(const std::string & key, const std::string & default_value){
    assert(sim_options.find(key) == sim_options.end()); // verify this option doesn't already exit
    
    SimOption o;
    o.key = key;
    o.type = STRING;
    o.str_value = default_value;
    
    sim_options[key] = o;
}

void SimOptions::addIntegerOption(const std::string & key, int default_value)
{
    assert(sim_options.find(key) == sim_options.end()); // verify this option doesn't already exit
    
    SimOption o;
    o.key = key;
    o.type = INTEGER;
    o.int_value = default_value;
    
    sim_options[key] = o;
}

void SimOptions::addDoubleOption(const std::string & key, double default_value)
{
    assert(sim_options.find(key) == sim_options.end()); // verify this option doesn't already exit
    
    SimOption o;
    o.key = key;
    o.type = DOUBLE;
    o.double_value = default_value;
    
    sim_options[key] = o;
}

void SimOptions::addBooleanOption(const std::string & key, bool default_value)
{
    assert(sim_options.find(key) == sim_options.end()); // verify this option doesn't already exit
    
    SimOption o;
    o.key = key;
    o.type = BOOLEAN;
    o.bool_value = default_value;
    
    sim_options[key] = o;
}

const std::string & SimOptions::strValue(const std::string & key)
{
    assert(sim_options.find(key) != sim_options.end()); // verify this option exists
    assert(sim_options[key].type == STRING);          // verify this option has the correct type
    return sim_options[key].str_value;
}

int SimOptions::intValue(const std::string & key)
{
    assert(sim_options.find(key) != sim_options.end()); // verify this option exists
    assert(sim_options[key].type == INTEGER);         // verify this option has the correct type
    return sim_options[key].int_value;
}

double SimOptions::doubleValue(const std::string & key)
{
    assert(sim_options.find(key) != sim_options.end()); // verify this option exists
    assert(sim_options[key].type == DOUBLE);          // verify this option has the correct type
    return sim_options[key].double_value;
}

bool SimOptions::boolValue(const std::string & key)
{
    assert(sim_options.find(key) != sim_options.end()); // verify this option exists
    assert(sim_options[key].type == BOOLEAN);         // verify this option has the correct type
    return sim_options[key].bool_value;
}


// adapted from Options::ParseOptionFile() from Da 2016 code
bool SimOptions::loadSimOptions(std::string infileName){
    // load sim options file
    std::ifstream infile(infileName);
    if (!infile.is_open()) { 
        std::cerr << "Unable to open options file!" << std::endl; 
        assert(!"Unable to open options file!");
    }

    std::string line;
    while(!infile.eof()){
        std::getline(infile, line);
        std::stringstream ss(line);

        std::string key;
        ss >> key;
        if (key == "#" || key == "" || ss.eof())    // skip comment lines and empty lines
            continue;
        
        std::map<std::string, SimOption>::iterator i = sim_options.find(key);
        if (i == sim_options.end())
        {
            std::cout << "Unrecognized option: " << key << " in option file " << infileName << "." << std::endl;
            assert(!"Unrecognized option in file");
        }

        switch (i->second.type)
        {
            case STRING:
                ss >> i->second.str_value;
                break;
            case INTEGER:
                ss >> i->second.int_value;
                break;
            case DOUBLE:
                ss >> i->second.double_value;
                break;
            case BOOLEAN:
                ss >> i->second.bool_value;
                break;
            default:
                assert(!"Unexpected option type");
                break;
        }
    }
    infile.close();
    return true;
}