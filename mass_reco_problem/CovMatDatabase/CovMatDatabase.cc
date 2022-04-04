#include "CovMatDatabase.h"

CovMatDatabase::CovMatDatabase(){
    last_input_file_name = "";
    last_output_file_name = "";
}

CovMatDatabase::CovMatDatabase(const std::string datafile){
    ReadDatabaseFromFile(datafile);
    last_input_file_name = datafile;
    last_output_file_name = "";
}

CovMatDatabase::~CovMatDatabase(){

}

void CovMatDatabase::AddCovMat(const std::string mat_name, const std::vector<float> matrix){
    data_map[mat_name] = matrix;
}

void CovMatDatabase::RemoveCovMat(const std::string mat_name){
    data_map.erase(mat_name);
    //data_map.erase(data_map.find(mat_name));
}

std::vector<float> CovMatDatabase::GetCovMat(const std::string mat_name) const{
    std::unordered_map<std::string, std::vector<float> >::const_iterator it = data_map.find(mat_name);
    if (it == data_map.end()){
        std::cout << "There is no matrix with name \"" << mat_name << "\" in the database.\n";
        std::cout << "Empty vector is returned.\n";
        return std::vector<float>();
    } else {
        return it->second;
    }
}

std::unordered_map<std::string, std::vector<float> > CovMatDatabase::GetDataMap() const {
    return data_map;
}

std::vector<std::string> CovMatDatabase::GetListOfMatNames() const {
    std::vector<std::string> list;
    std::unordered_map<std::string, std::vector<float> >::const_iterator iter = data_map.begin();
    while(iter!=data_map.end()){
        list.push_back(iter->first);
        ++iter;
    }
    return list;
}

void CovMatDatabase::ReadDatabaseFromFile(const std::string covmat_database){
    std::ifstream input_file(covmat_database);
    if(!input_file){
		std::cout << "Cannot open file \"" << covmat_database << "\": no such file. Exit." << std::endl;
		return;
	}
    if (!data_map.empty())
        data_map.clear();

    std::string mat_name = "";
    size_t num_of_indep_el = 0;

    std::string s;
    while(getline(input_file, s)){
        if (s.size() == 0)
            continue;
        if (s[0] == '#')
            continue;
        if (s[0] == '%')
            continue;
        if (s[0] == '/' && s[1] == '/')
            continue;
        std::stringstream ss(s);
        ss >> mat_name >> num_of_indep_el;
        std::vector<float> cov_mat(num_of_indep_el);
        for (size_t i = 0; i < num_of_indep_el; i++){
            ss >> cov_mat[i];
        }
        data_map[mat_name] = cov_mat;
    }
    input_file.close();
    last_input_file_name = covmat_database;
}

void CovMatDatabase::PrintCovMat(const std::string mat_name) const {
    std::cout << "\n  Covariance matrix \"" << mat_name << "\": " << std::endl;
    size_t max_row_counter = 1;
    size_t present_row_number = 0;
    if (data_map.count(mat_name) == 0){
        std::cout << "There is no matrix with name \"" << mat_name << "\" in the database.\n";
        return;
    }
    for (int i = 0; i < data_map.at(mat_name).size(); i++){
        std::cout << data_map.at(mat_name).at(i) << std::setw(15);
        present_row_number++;
        if (present_row_number == max_row_counter){
            std::cout << "\n";
            max_row_counter++;
            present_row_number = 0;
        }            
    }
    std::cout << "\n";
}

void CovMatDatabase::PrintListOfMatNames() const {
    std::cout << "\n  List of matrixes:  " << std::endl;
    std::unordered_map<std::string, std::vector<float> >::const_iterator iter = data_map.begin();
    while(iter!=data_map.end()){
        std::cout << iter->first << std::endl;
        ++iter;
    }
    return;
}

void CovMatDatabase::PrintDatabase() const {
    std::cout << "  Database:  " << std::endl;
    std::unordered_map<std::string, std::vector<float> >::const_iterator iter = data_map.begin();
    while(iter!=data_map.end()){
        std::cout << iter->first << ":  \n";
        PrintCovMat(iter->first);
        ++iter;
    }
    return;
}

void CovMatDatabase::SaveDatabaseToFile(const std::string output_name){
    std::ofstream out;
    out.open(output_name);
    if (out.is_open())
    {   
        out << "# Name     Number of independent elements of matrix       Matrix elements\n";
        std::unordered_map<std::string, std::vector<float> >::const_iterator iter = data_map.begin();
        while(iter!=data_map.end()){
            out << iter->first << "     ";
            out << iter->second.size() << "     ";
            if (iter->second.size() == 0){
                std::cout << "Warning! You are trying to save 0-dimensional covariance matrix \"" << iter->first << "\" !" << std::endl;
                ++iter;
                continue;
            }
            if (iter->second.size() == 1){
                out << iter->second.at(0) << "\n";
                ++iter;
                continue;
            }
            for (size_t i = 0; i < iter->second.size() - 1; i++){
                out << iter->second.at(i) << "  ";
            }
            out << iter->second.at(iter->second.size() - 1) << "\n";
            ++iter;
        }   
    }
    std::cout << "Covariance matrixes database saved to " << output_name << "\n";
    out.close(); 
    last_output_file_name = output_name;
    return;
}

std::string CovMatDatabase::GetLastInputFileName() const {
    return last_input_file_name;
}

std::string CovMatDatabase::GetLastOutputFileName() const{
    return last_output_file_name;
}

bool CovMatDatabase::CheckMatrix(std::string mat_name, bool silent) const{
    std::unordered_map<std::string, std::vector<float> >::const_iterator matrix_it = data_map.find(mat_name);
    if (matrix_it == data_map.end()){
        if (!silent)
            std::cout << "No matrix with name \"" << mat_name << "\" in the database. Return false.\n";
        return false;
    }
    size_t l = (matrix_it->second).size();
    float n = (-1 + sqrt(1+8*l))/2;
    if (n != floor(n)){
        if (!silent)
            std::cout << "Warning! Covariance matrix \"" << mat_name << "\" has invalid format. Its length L = " << l << ", hence it cannot be interpreted as a square matrix. A square simmetrical NxN matrix must be converted into vector with the length L=(N+1)*N/2.\n";
        return false;
    }

    size_t n_int = (size_t)n;
    size_t row_num = 0;    
    while(row_num < n_int){
        size_t diag_el_num = (row_num+1)*row_num/2 + row_num;
        if((matrix_it->second).at(diag_el_num)<0){
            if (!silent)
                std::cout << "Warning! Covariance matrix \"" << mat_name << "\" has negative diagonal element #" << diag_el_num << " with the value " << (matrix_it->second).at(diag_el_num) << ". Please, check the matrix.\n";
            return false;
        }
        row_num++;
    }
    if (!silent)
        std::cout << "Covariance matrix \"" << mat_name << "\" is OK!\n";
    return true;
}

bool CovMatDatabase::CheckDatabase(bool silent) const{
    if (!silent)
        std::cout << "\n    Checking the database matrixes format...\n";
    bool flag = true;
    std::unordered_map<std::string, std::vector<float> >::const_iterator iter = data_map.begin();
    while(iter!=data_map.end()){
        if (!CheckMatrix(iter->first,silent))
            flag = false;
        ++iter;
    }
    if (flag){
        if (!silent)
            std::cout << "\nAll matrixes in the database are OK! (they're square and have non-negative diagonal elements only)\n";
    } else {
        if (!silent)
            std::cout << "\nSome matrixes have invalid format and/or negative diagonal elements. Check the information above.\n";
    }
    return flag;
}