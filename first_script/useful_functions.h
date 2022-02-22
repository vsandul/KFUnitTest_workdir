#pragma once

#include <vector>
#include <cmath>
#include <iostream>

double GetVectorLength(std::vector<double>& vec){
    double res = 0;
    for (const auto& i:vec){
        res += i*i;
    }

    return sqrt(res);
}

std::vector<double> GetUnitVector(double x, double y, double z){
    std::vector<double> unit_vector = {x ,y, z};
    double norm = GetVectorLength(unit_vector);

    for (auto& i:unit_vector){
        i /= norm;
    }

    return unit_vector;
}

double GetAngleBetweenVectors(std::vector<double>& v1, std::vector<double>& v2){
    if (v1.size() != v2.size()){
        std::cout << "ERROR! The vectors have different sizes. Returning 0." << std::endl;
        return 0;
    } else {
        if (v1.size() == 0){
            std::cout << "ERROR! The vectors are zero-dimensional. Cannot calculate angle and returning 0." << std::endl;
        return 0;
        } else { 
            double result = 0;
            double v1_len = 0;
            double v2_len = 0;
            for (int i = 0; i < v1.size(); i++){
                result += v1.at(i)*v2.at(i);
                v1_len += v1.at(i)*v1.at(i);
                v2_len += v2.at(i)*v2.at(i);
            }

            if (v1_len*v1_len == 0){
                std::cout << "ERROR! One of the vectors has zero-length. Cannot calculate angle and returning 0." << std::endl;
                return 0;
            } else {
                return acos(result/sqrt(v1_len)/sqrt(v2_len));
            }
        }
    }


}