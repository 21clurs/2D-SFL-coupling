// simoptions.h
// modeled after SimOptions.h in Fang Da's implementation of "Surface-Only Liquids"

#ifndef SIM_OPTIONS_H
#define SIM_OPTIONS_H

#include <iostream>
#include <string>
#include <vector>
#include <map>

class SimOptions{
    public:
        enum Type
        {
            STRING,
            INTEGER,
            DOUBLE,
            BOOLEAN
        };

        static void addStringOption(const std::string & key, const std::string & default_value);
        static void addIntegerOption(const std::string & key, int defaut_value);
        static void addDoubleOption(const std::string & key, double default_value);
        static void addBooleanOption(const std::string & key, bool default_value);

        static const std::string & strValue(const std::string & key);
        static int                 intValue(const std::string & key);
        static double              doubleValue(const std::string & key);
        static bool                boolValue(const std::string & key);

        static bool loadSimOptions(std::string infileName);
    
    protected:
        class SimOption
        {
        public:
            std::string key;        // option key
            Type type;              // option type

            std::string str_value; 
            int         int_value;
            double      double_value;
            bool        bool_value;
        };
        
        // all the sim options!!
        static std::map<std::string, SimOption> sim_options;
};

#endif