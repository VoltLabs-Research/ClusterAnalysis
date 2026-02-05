#pragma once

#include <opendxa/core/opendxa.h>
#include <opendxa/core/lammps_parser.h>
#include <opendxa/core/particle_property.h>
#include <opendxa/analysis/cluster_analysis.h>
#include <nlohmann/json.hpp>
#include <memory>
#include <string>

namespace OpenDXA{
using json = nlohmann::json;

class ClusterAnalysisWrapper{
public:
    ClusterAnalysisWrapper();

    void setCutoff(double cutoff);
    void setOptions(
        bool sortBySize,
        bool unwrapParticleCoordinates,
        bool computeCenterOfMass,
        bool computeRadiusOfGyration
    );

    json compute(const LammpsParser::Frame& frame, const std::string &outputFilename);

private:
    std::shared_ptr<Particles::ParticleProperty> createPositionProperty(const LammpsParser::Frame& frame);
    double _cutoff;
    bool _sortBySize;
    bool _unwrapParticleCoordinates;
    bool _computeCentersOfMass;
    bool _computeRadiusOfGyration;
};

}