//
//  contactNumber.hpp
//  gmx_contactNumber
//
//  Created by Yiming Tang on 17/03/2018.
//  Copyright © 2018 Yiming Tang. All rights reserved.
//

#ifndef contactNumber_hpp
#define contactNumber_hpp

#include <stdio.h>
#include <gromacs/trajectoryanalysis.h>
#include <string>


using namespace std;
using namespace gmx;

class contactNumber : public TrajectoryAnalysisModule
{
public:
    contactNumber();
    
    virtual void initOptions(IOptionsContainer             *options,
                             TrajectoryAnalysisSettings    *settings);
    
    virtual void initAnalysis(const TrajectoryAnalysisSettings  &settings,
                              const TopologyInformation         &top);
    
    virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                              TrajectoryAnalysisModuleData *pdata);
    
    virtual void finishAnalysis(int nframes);
    
    virtual void writeOutput();
    
private:
    
    class ModuleData;
    
    bool isCarbon(string testString);
    bool isHydrogen(string testString);
    
    AnalysisData                     dataContact_;
    AnalysisDataAverageModulePointer avem_;
    
    std::string      fnContact_;
    
    double           cutoff_carbon_;
    double           cutoff_noncarbon_;
    
    bool             consider_hydrogen_;
    
    Selection        ref_;
    SelectionList    sel_;
    
    t_topology      *top_;
    t_atoms          atoms_;
    
};


#endif /* contactNumber_hpp */
