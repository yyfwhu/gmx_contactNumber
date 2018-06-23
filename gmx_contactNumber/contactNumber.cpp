//
//  contactNumber.cpp
//  gmx_contactNumber
//
//  Created by Yiming Tang on 17/03/2018.
//  Copyright © 2018 Yiming Tang. All rights reserved.
//

#include "contactNumber.hpp"
#include <gromacs/fileio/matio.h>
#include <gromacs/fileio/gmxfio.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <string.h>

using namespace std;
using namespace gmx;

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
        "Contact Number Calculator 2.0 by Yiming Tang @ Fudan.\n",
        "I hate manuals so I'm not going to talk about how this program is constructed.",
        "Nevertheless, it is useful to tell you a little about how to make full use of",
        "this program to calculate contact number between groups.\n",
        "Firstly, you should use command line to select three sets of groups: (1) -all: this",
        "should always contain all atoms in your trajectory file. (2) -reference and (3) -select:",
        "these two sets of groups support multivalue-selection. The contacts are calculated between",
        "each group-pair between these two sets reference1-select1, reference1-select2, ...,",
        "reference2-select1, etc.\n",
        "Secondly, you should tell the program what output you want. Three kinds of output files can",
        "be processed: (1) -map: the frame-averaged contact number/probability map in XPM file which can be further transfered",
        "into eps using gmx xpm2ps command. (3) -verbose: the frame-averaged contact number/probability between each group-pair as",
        "a function of time.\n",
        "If -probability flag is specified, a Contact Probability Map is calculated instead of contact number. ",
        "In Contact Probability Map, two residues are said to be in contact if there exist at least one contact ",
        "between them.\n",
        "NOTICE: As this program contains two group-selections that support multivalue, you should always",
        "use conmmand line option of select instead of using user-interface.\n",
        "EXAMPLE: gmx_contactNumber -f traj.xtc -s topology.tpr -all protein -reference $(seq 7 206) -select",
        "$(seq 7 206), -verbose contactnumber.\n",
        "For further information contact Yiming Tang (tym@tymworld.com) but there's significant chance that I",
        "ignore you. Bite me!"
    };

    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    
    settings->setHelpText(desc);
    
    options->addOption(BooleanOption("probability")
                       .store(&probability_).required().defaultValue(false)
                       .description("Calculate contact probability instead of contact number"));
    
    options->addOption(FileNameOption("verbose")
                       .filetype(eftPlot).outputFile()
                       .store(&fnContact_).defaultBasename("contact_number")
                       .description("contact number as a function of time"));
    
    options->addOption(FileNameOption("map").filetype(eftUnknown).legacyType(efXPM).outputFile()
                       .store(&fnMap_).description("Fraction of Contact Number Map"));
    
    options->addOption(FileNameOption("mapraw").filetype(eftGenericData).outputFile()
                       .store(&fnMapRaw_).description("Raw data for Fraction of Contact Number Map"));
    
    options->addOption(SelectionOption("all")
                       .store(&all_).required()
                       .description("Group containing all atoms in trajectory for neighborhood searching"));
    
    options->addOption(SelectionOption("reference")
                       .storeVector(&ref_).required().multiValue()
                       .description("Reference groups to calculate contact from"));
    
    options->addOption(SelectionOption("select")
                       .storeVector(&sel_).required().multiValue()
                       .description("Groups to calculate contacts to"));
    
    options->addOption(BooleanOption("hydrogen")
                       .store(&consider_hydrogen_).required().defaultValue(false)
                       .description("Whether hydrogen atoms are taken into consideration"));
    
    options->addOption(BooleanOption("inter_residue")
                       .store(&not_same_residue_).required().defaultValue(false)
                       .description("Whether only contacts between different residues are considered."));
    
    options->addOption(DoubleOption("cc").store(&cutoff_carbon_)
                       .description("Cutoff for carbon-carbon interactions"));
    
    options->addOption(DoubleOption("cnc").store(&cutoff_noncarbon_)
                       .description("Cutoff for carbon-noncarbon and noncarbon-noncarbon interactions"));
    
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    
    
}

void contactNumber::initAnalysis(const gmx::TrajectoryAnalysisSettings &settings,
                                 const gmx::TopologyInformation &top)
{
    // IMPORTANT: dataset count: sel_; column count: ref_
    //cout << sel_.size() << '\t' << ref_.size() << '\t' << (int)sel_.size() * (int)ref_.size() << endl;
    
    dataContact_.setColumnCount(0, (int)sel_.size() * (int)ref_.size());
    //dataContact_.setColumnCount(0, 200);
    

    if (!fnContact_.empty())
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnContact_);
        plotm->setTitle("Contact Number");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Contact Number");
        
        //cout << "Set:filename,title,x,y" << endl;
        
        string **setName = new string*[sel_.size()];
        
        for(int gsel = 0; gsel < sel_.size(); gsel++ )
        {
            
            setName[gsel] = new string[ref_.size()];
            
            for(int gref = 0; gref < ref_.size(); gref++)
            {
                setName[gsel][gref] = sel_[gsel].selectionText();
                setName[gsel][gref].append("-");
                setName[gsel][gref].append(ref_[gref].selectionText());

                plotm->appendLegend(setName[gsel][gref]);
            }
        }


        dataContact_.addModule(plotm);
    }
    
    //cout << "finished: set plot" << endl;
    
    avemContact_.reset(new AnalysisDataAverageModule());
    dataContact_.addModule(avemContact_);
    
    top_   = top.topology();
    atoms_ = top.topology()->atoms;
    
    nb_.setCutoff(max(cutoff_carbon_,cutoff_noncarbon_));
    /*
    // Set original IDs in all in accordance to "selection"
    for (int i_all = 0; i_all < all_.atomCount(); i_all++)
    {
        all_.setOriginalId(i_all, (int)sel_.size() + 2);
    }
    
    for (int gsel = 0; gsel < sel_.size(); gsel++ )
    {
        for(ConstArrayRef<int>::iterator iter_atomIndice_sel = sel_[gsel].atomIndices().begin();
            iter_atomIndice_sel != sel_[gsel].atomIndices().end(); iter_atomIndice_sel++)
        {
            all_.setOriginalId(*iter_atomIndice_sel, gsel);
        }
    }
    // Finish: Set original IDs in all in accordance to "selection"
    */
    //cout << "FINISH: initialize" << endl;
    
    
}

void contactNumber::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                 TrajectoryAnalysisModuleData *pdata)
{
    
    // Initialize temp contact number matrix
    int **tempContactNumber = new int*[sel_.size()];
    for(int i = 0; i<sel_.size(); i++)
    {
        tempContactNumber[i] = new int[ref_.size()];
        for(int j = 0; j<ref_.size(); j++)
        {
            tempContactNumber[i][j] = 0;
        }
    }
    // END: Initialize temp contact number matrix
    
    AnalysisDataHandle         dh       = pdata->dataHandle(dataContact_);
    dh.startFrame(frnr, fr.time);
    
    for (int gsel = 0; gsel < sel_.size(); gsel++)
    {

        AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, sel_[gsel]);
        
        for (int gref = 0; gref < ref_.size(); gref++)
        {
            ConstArrayRef<rvec>::iterator iter_coordinate_ref;
            ConstArrayRef<int>::iterator  iter_atomIndice_ref;
            
            for(iter_coordinate_ref  = ref_[gref].coordinates().begin(), iter_atomIndice_ref  = ref_[gref].atomIndices().begin();
                iter_coordinate_ref != ref_[gref].coordinates().end() && iter_atomIndice_ref != ref_[gref].atomIndices().end();
                iter_coordinate_ref++, iter_atomIndice_ref++)
            {
                AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(*iter_coordinate_ref);
                AnalysisNeighborhoodPair pair;
                while(pairSearch.findNextPair(&pair))
                {
                    
                    string refName  = *(atoms_.atomname[*iter_atomIndice_ref]);
                    string selName  = *(atoms_.atomname[sel_[gsel].atomIndices()[pair.refIndex()]]);
                    double distance = sqrt(pair.distance2());
                    //cout << refName << '\t' << selName << '\t' << distance << endl;
                    
                    // Stript out contacts from same residue
                    if(not_same_residue_ && (atoms_.atom[*iter_atomIndice_ref].resind) == atoms_.atom[sel_[gsel].atomIndices()[pair.refIndex()]].resind)
                    {
                        continue;
                    }
                    // FINISH: Stript out contacts from same residue
                    
                    if( consider_hydrogen_ || ( !isHydrogen(refName) && !isHydrogen(selName) ) )
                    {
                        if( isCarbon(refName) && isCarbon(selName))
                        {
                            if(distance <= cutoff_carbon_)
                            {
                                //cout << atoms_.atom[*iter_atomIndice_ref].resind << '\t' << *(atoms_.atomname[*iter_atomIndice_ref]) << '\t'
                                //     << atoms_.atom[sel_[gsel].atomIndices()[pair.refIndex()]].resind << '\t'
                                //     << *(atoms_.atomname[sel_[gsel].atomIndices()[pair.refIndex()]]) << '\t' << distance << endl;
                                tempContactNumber[gsel][gref] += 1;
                            }
                        }
                        else
                        {
                            if(distance <= cutoff_noncarbon_)
                            {
                                //cout << atoms_.atom[*iter_atomIndice_ref].resind << '\t' << *(atoms_.atomname[*iter_atomIndice_ref]) << '\t'
                                //<< atoms_.atom[sel_[gsel].atomIndices()[pair.refIndex()]].resind << '\t'
                                //<< *(atoms_.atomname[sel_[gsel].atomIndices()[pair.refIndex()]]) << '\t' << distance << endl;
                                tempContactNumber[gsel][gref] += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    
    // If two groups are identical, divided it by 2
    
    for (int gsel = 0; gsel < sel_.size(); gsel++)
    {
        for (int gref = 0; gref < ref_.size(); gref++)
        {
            //cout << sel_.size() << '\t' << (int)ref_.size() << '\t' << ref_.size() << '\t' << (int)ref_.size() << endl;
            
            if(probability_)
            {
                dh.setPoint(gsel * (int)ref_.size() + gref, tempContactNumber[gsel][gref] > 0 ? 1 : 0);
            }
            else
            {
                if ( strcmp(ref_[gref].selectionText(), sel_[gsel].selectionText()) == 0)
                {
                    //dh.selectDataSet(gsel);
                    //dh.setPoint(gref, tempContactNumber[gsel][gref] / 2) ;
                    dh.setPoint(gsel * (int)ref_.size() + gref, tempContactNumber[gsel][gref] / 2);
                }
                else
                {
                    //dh.selectDataSet(gsel);
                    //dh.setPoint(gref, tempContactNumber[gsel][gref]);
                    dh.setPoint(gsel * (int)ref_.size() + gref, tempContactNumber[gsel][gref]);
                }
                //cout << gsel * (int)ref_.size() + gref << '\t' << tempContactNumber[gsel][gref] / 2 << endl;
            }
            

        }
    }
    dh.finishFrame();
}


void contactNumber::finishAnalysis(int /*nframes*/)
{
    
}

void contactNumber::writeOutput()
{
    
    if(!fnMap_.empty() || !fnMapRaw_.empty())
    {
        /*
        // Get Sum of Contact Numbers
        real sum_contact = 0;
        for(int gsel = 0; gsel < sel_.size(); gsel++)
        {
            for (int gref = 0; gref < ref_.size(); gref++)
            {
                sum_contact += avemContact_->average(0, gsel * (int)ref_.size() + gref);
            }
        }*/
        // Get matrix for Averaged Contact Number/Probability
        real **matContact = new real*[sel_.size()];
        for(int i = 0; i<sel_.size(); i++)
        {
            matContact[i] = new real[ref_.size()];
        }
        
        for(int gsel = 0; gsel < sel_.size(); gsel++)
        {
            for (int gref = 0; gref < ref_.size(); gref++)
            {
                matContact[gsel][gref] = avemContact_->average(0, gsel * (int)ref_.size() + gref);
            }
        }
        
        // Now output!
        
        if(!fnMap_.empty())
        {
            if(fnMap_.compare(".xpm"))  {}
            else                        {fnMap_ += ".xpm";}
            
            FILE *fpMapXPM = fopen(fnMap_.c_str(), "w");
            
            t_rgb rlo, rhi;
            rlo.r = 0.0; rlo.g = 0.0; rlo.b = 1.0;
            rhi.r = 1.0; rhi.g = 0.0; rhi.b = 0.0;
            int nlevels  = 400;
            
            real *refVector = new real[ref_.size()];
            real *selVector = new real[sel_.size()];
            
            for(int i = 0; i<ref_.size(); i++)
            {
                refVector[i] = i;
            }
            
            for(int i = 0; i<sel_.size(); i++)
            {
                selVector[i] = i;
            }
            
            int color_max = probability_ ? 1 : 10;

            write_xpm(fpMapXPM, 0, "Contact Number"
                      , "Contact Number", "Reference Group Number", "Selection Group Number"
                      , (int)(sel_.size()), (int)(ref_.size()), selVector, refVector
                      , matContact, 0, color_max, rlo, rhi, &nlevels);
            
            fclose(fpMapXPM);
        }
        
        if(!fnMapRaw_.empty())
        {
            if(fnMapRaw_.compare(".dat"))  {}
            else                        {fnMapRaw_ += ".dat";}
            
            FILE *fpMapRaw = fopen(fnMapRaw_.c_str(), "w");
            
            fprintf(fpMapRaw, "Selection Groups: ");
            for(int i = 0; i<sel_.size(); i++)
            {
                fprintf(fpMapRaw, sel_[i].selectionText());
                fprintf(fpMapRaw, "\t");
            }
            fprintf(fpMapRaw, "\n");
            
            fprintf(fpMapRaw, "Reference Groups: ");
            for(int i = 0; i<ref_.size(); i++)
            {
                fprintf(fpMapRaw, ref_[i].selectionText());
                fprintf(fpMapRaw, "\t");
            }
            
            
            for(int gsel = 0; gsel < sel_.size(); gsel++)
            {
                for (int gref = 0; gref < ref_.size(); gref++)
                {
                    fprintf(fpMapRaw, "%.3f ", matContact[gsel][gref]);
                }
                fprintf(fpMapRaw, "\n");
            }

            
            fprintf(fpMapRaw, "\n");
            
        }
        
        
    }
    

}


