#pragma once

#include <iostream>
#include <iomanip>  
#include <fstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

class CovMatDatabase{
    public:
        CovMatDatabase();
        CovMatDatabase(const std::string datafile);
        ~CovMatDatabase();

        void ReadDatabaseFromFile(const std::string covmat_database);
        void SaveDatabaseToFile(const std::string output_name);

        void AddCovMat(const std::string mat_name, const std::vector<float> matrix);
        void RemoveCovMat(const std::string mat_name);

        std::vector<float> GetCovMat(const std::string mat_name) const;        
        std::unordered_map<std::string, std::vector<float> > GetDataMap() const;  
        std::vector<std::string> GetListOfMatNames() const;      

        void PrintCovMat(const std::string mat_name) const;
        void PrintListOfMatNames() const;
        void PrintDatabase() const;

        std::string GetLastInputFileName() const;
        std::string GetLastOutputFileName() const;

        bool CheckMatrix(std::string mat_name, bool silent=false) const;
        bool CheckDatabase(bool silent=false) const;
        
    private:
        std::unordered_map<std::string, std::vector<float> > data_map;
        std::string last_input_file_name;
        std::string last_output_file_name;

};