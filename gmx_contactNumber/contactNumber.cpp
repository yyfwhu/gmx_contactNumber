//
//  contactNumber.cpp
//  gmx_contactNumber
//
//  Created by Yiming Tang on 17/03/2018.
//  Copyright © 2018 Yiming Tang. All rights reserved.
//

#include "contactNumber.hpp"
#include <iostream>

using namespace std;

bool contactNumber::isCarbon(string testString)
{
    if (testString[0] == 'C')
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool contactNumber::isHydrogen(string testString)
{
    if (testString[0] == 'H')
    {
        return true;
    }
    else
    {
        return false;
    }
}

contactNumber::contactNumber() : cutoff_carbon_(0.54), cutoff_noncarbon_(0.46)
{
    registerAnalysisDataset(&dataContact_, "contactNumber");
}


void contactNumber::initOptions(gmx::IOptionsContainer  *options, gmx::TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] =
    {
        "Contact Number Calculater by Yiming Tang @ Fudan. \n",
        "This tool calculates contact numbers as a function of time between",
        "two groups identified as reference group (ref) and selection group (sel).",
        "The default tolerable range of exisiting contacts is 0.54 nm for carbon-carbon",
        "interactions and 0.46 nm for carbon-noncarbon or noncarbon-noncarbon",
        "interactions.\n",
        "IMPORTANT: For non-protein, element field is usually missing so that you should,",
        "set -noproteinonly flag! "
    };
    
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    
    settings->setHelpText(desc);
    
    options->addOption(BooleanOption("proteinonly")
                       .store(&only_protein_).required().defaultValue(true)
                       .description("If your selection contains only protein. "));
    
    options->addOption(FileNameOption("o")
                       .filetype(eftPlot).outputFile()
                       .store(&fnContact_).defaultBasename("contact_number")
                       .description("contact number as a function of time"));
    
    options->addOption(SelectionOption("reference")
                       .store(&ref_).required()
                       .description("Reference group to calculate contact from"));
    
    options->addOption(SelectionOption("select")
                       .storeVector(&sel_).required().multiValue()
                       .description("Groups to calculate contacts to"));
    
    options->addOption(BooleanOption("hydrogen")
                       .store(&consider_hydrogen_).required().defaultValue(false)
                       .description("Whether hydrogen atoms are taken into consideration"));
    
    options->addOption(DoubleOption("cc").store(&cutoff_carbon_)
                       .description("Cutoff for carbon-carbon interactions"));
    
    options->addOption(DoubleOption("cnc").store(&cutoff_noncarbon_)
                       .description("Cutoff for carbon-noncarbon and noncarbon-noncarbon interactions"));
    
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    
}

void contactNumber::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
                                 const gmx::TopologyInformation &top)
{
    dataContact_.setColumnCount(0, (int)sel_.size());
    
    avem_.reset(new AnalysisDataAverageModule());
    dataContact_.addModule(avem_);
    
    if (!fnContact_.empty())
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnContact_);
        plotm->setTitle("Contact Number");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Contact Number");
        dataContact_.addModule(plotm);
    }
    
    top_   = top.topology();
    atoms_ = top.topology()->atoms;
    
}

void contactNumber::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                               TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle        dh     = pdata->dataHandle(dataContact_);
    const Selection          &ref    = pdata->parallelSelection(ref_);
    t_atoms                   atoms  = atoms_;
    
    double cutoff_carbon    = cutoff_carbon_;
    double cutoff_noncarbon = cutoff_noncarbon_;
    
    dh.startFrame(frnr, fr.time);
    
    double tempDistance;
    double tempContactNumber;
    string elemRef;
    string elemSel;
    
    
    for (size_t g = 0; g < sel_.size(); ++g)
    {
        tempContactNumber = 0;
        
        const Selection &sel   = pdata->parallelSelection(sel_[g]);
        
        ConstArrayRef<rvec>::iterator iter_coordinate_ref;
        ConstArrayRef<int>::iterator  iter_atomIndice_ref;
        ConstArrayRef<rvec>::iterator iter_coordinate_sel;
        ConstArrayRef<int>::iterator  iter_atomIndice_sel;
        
        for(iter_coordinate_ref  = ref.coordinates().begin(), iter_atomIndice_ref  = ref.atomIndices().begin();
            iter_coordinate_ref != ref.coordinates().end() && iter_atomIndice_ref != ref.atomIndices().end();
            ++iter_coordinate_ref, ++iter_atomIndice_ref)
        {
            for(iter_coordinate_sel  = sel.coordinates().begin(), iter_atomIndice_sel  = sel.atomIndices().begin();
                iter_coordinate_sel != sel.coordinates().end() && iter_atomIndice_sel != sel.atomIndices().end();
                ++iter_coordinate_sel, ++iter_atomIndice_sel)
            {
                tempDistance = sqrt(  pow(iter_coordinate_ref[0][0] - iter_coordinate_sel[0][0], 2.0)
                                    + pow(iter_coordinate_ref[0][1] - iter_coordinate_sel[0][1], 2.0)
                                    + pow(iter_coordinate_ref[0][1] - iter_coordinate_sel[0][2], 2.0));
                
                // printf("%f\n", tempDistance);

                if(only_protein_)
                {
                    elemRef = atoms.atom[*iter_atomIndice_ref].elem;
                    elemSel = atoms.atom[*iter_atomIndice_sel].elem;
                }
                else
                {
                    elemRef = *atoms.atomname[*iter_atomIndice_ref];
                    elemSel = *atoms.atomname[*iter_atomIndice_sel];
                }
                
                
                if( consider_hydrogen_ || ( !isHydrogen(elemRef) && !isHydrogen(elemSel)))
                {
                    /*cout<<elemRef;
                    printf("\t");
                    cout<<elemSel;
                    printf("\n");*/
                    
                    if ( isCarbon(elemRef) && isCarbon(elemSel))
                    {
                        if(tempDistance <= cutoff_carbon) { tempContactNumber++;  }
                    }
                    else
                    {
                        if(tempDistance <= cutoff_noncarbon) { tempContactNumber++;}
                    }
                }
            }
        }
        dh.setPoint((int)g, tempContactNumber);
    }
    
    dh.finishFrame();
}

void contactNumber::finishAnalysis(int /*nframes*/)
{
}

void contactNumber::writeOutput()
{
    for (size_t g = 0; g < sel_.size(); ++g)
    {
        fprintf(stderr, "Average contact number from %s to '%s': %.3f nm\n",
                ref_.name(), sel_[g].name(), avem_->average(0, (int)g));
    }
}

    


                
     
