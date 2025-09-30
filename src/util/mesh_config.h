#ifndef MESH_CONFIG_H
#define MESH_CONFIG_H

#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include "util/constants.h"

/**
 * Simple JSON parser for mesh configuration.
 * Reads mesh_config.json and provides mesh parameters.
 */
class MeshConfig {
private:
    int mulFac;
    double Re;  // Reynolds number
    int maxIterations;  // Maximum solver iterations
    double bg_x0, bg_y0, bg_theta_deg, bg_length, bg_width;
    int bg_Nx_base, bg_Ny_base;
    double comp_x0, comp_y0, comp_theta_deg, comp_length, comp_width;
    int comp_Nx_base, comp_Ny_base;

    // Simple helper to extract value from JSON line like "  "key": value,"
    template<typename T>
    T extractValue(const std::string& line) {
        size_t colonPos = line.find(':');
        if (colonPos == std::string::npos) return T();
        
        std::string valueStr = line.substr(colonPos + 1);
        // Remove spaces, commas, quotes using algorithm
        valueStr.erase(std::remove_if(valueStr.begin(), valueStr.end(), 
                       [](char c) { return c == ' ' || c == ',' || c == '"'; }), 
                       valueStr.end());
        
        std::istringstream iss(valueStr);
        T value;
        iss >> value;
        return value;
    }

    void parseConfig(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open mesh_config.json");
        }

        std::string line;
        std::string currentMesh = "";
        
        while (std::getline(file, line)) {
            // Detect which mesh section we're in
            if (line.find("\"bg\"") != std::string::npos) {
                currentMesh = "bg";
            } else if (line.find("\"comp\"") != std::string::npos) {
                currentMesh = "comp";
            }
            
            // Parse mulFac
            if (line.find("\"mulFac\"") != std::string::npos) {
                mulFac = extractValue<int>(line);
            }
            
            // Parse Re (Reynolds number)
            if (line.find("\"Re\"") != std::string::npos) {
                Re = extractValue<double>(line);
            }
            
            // Parse maxIterations
            if (line.find("\"maxIterations\"") != std::string::npos) {
                maxIterations = extractValue<int>(line);
            }
            
            // Parse background mesh parameters
            if (currentMesh == "bg") {
                if (line.find("\"x0\"") != std::string::npos) bg_x0 = extractValue<double>(line);
                else if (line.find("\"y0\"") != std::string::npos) bg_y0 = extractValue<double>(line);
                else if (line.find("\"theta_deg\"") != std::string::npos) bg_theta_deg = extractValue<double>(line);
                else if (line.find("\"length\"") != std::string::npos) bg_length = extractValue<double>(line);
                else if (line.find("\"width\"") != std::string::npos) bg_width = extractValue<double>(line);
                else if (line.find("\"Nx_base\"") != std::string::npos) bg_Nx_base = extractValue<int>(line);
                else if (line.find("\"Ny_base\"") != std::string::npos) bg_Ny_base = extractValue<int>(line);
            }
            
            // Parse component mesh parameters
            if (currentMesh == "comp") {
                if (line.find("\"x0\"") != std::string::npos) comp_x0 = extractValue<double>(line);
                else if (line.find("\"y0\"") != std::string::npos) comp_y0 = extractValue<double>(line);
                else if (line.find("\"theta_deg\"") != std::string::npos) comp_theta_deg = extractValue<double>(line);
                else if (line.find("\"length\"") != std::string::npos) comp_length = extractValue<double>(line);
                else if (line.find("\"width\"") != std::string::npos) comp_width = extractValue<double>(line);
                else if (line.find("\"Nx_base\"") != std::string::npos) comp_Nx_base = extractValue<int>(line);
                else if (line.find("\"Ny_base\"") != std::string::npos) comp_Ny_base = extractValue<int>(line);
            }
        }
        file.close();
    }

public:
    MeshConfig(const std::string& filename = "mesh_config.json") {
        // Set default values
        maxIterations = 60000;  // Default max iterations
        Re = 100.0;  // Default Reynolds number
        
        parseConfig(filename);
    }

    // Getters for background mesh
    double getBgX0() const { return bg_x0; }
    double getBgY0() const { return bg_y0; }
    double getBgTheta() const { return bg_theta_deg * PI / 180.0; }
    double getBgLength() const { return bg_length; }
    double getBgWidth() const { return bg_width; }
    int getBgNx() const { return bg_Nx_base * mulFac; }
    int getBgNy() const { return bg_Ny_base * mulFac; }

    // Getters for component mesh
    double getCompX0() const { return comp_x0; }
    double getCompY0() const { return comp_y0; }
    double getCompTheta() const { return comp_theta_deg * PI / 180.0; }
    double getCompLength() const { return comp_length; }
    double getCompWidth() const { return comp_width; }
    int getCompNx() const { return comp_Nx_base * mulFac; }
    int getCompNy() const { return comp_Ny_base * mulFac; }

    int getMulFac() const { return mulFac; }
    double getReynolds() const { return Re; }
    int getMaxIterations() const { return maxIterations; }
};

#endif // MESH_CONFIG_H
