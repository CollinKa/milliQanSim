//goal: Covert sim PMT number to digitizer channel number. Do I need to convert the sim channel number to offline channel number?
//Jan 30, 2026, Collin Zheng
//todo follow pmtPosition in flatlightwithphotonslab.py to get the front back pair in the offline summing amp channel.(done)
//create the mapping function for scintillator number to offline channel number.(Slab number * 2 + type number) by comparing the SLAB number map, it matches tht trigger board number map in google sheet.(done)
//create the pulse height scaling for each channel->            double scaleFactor = areaSum / (828.03) * (44.8573/59.0878); //sample area, then scale by the pulse height for SPEs in data
//Need to fix the max_pulse_height for each channel.(no rush like calibration factor)
//add Data for the SPE pulse height study: & fix scaling facor follow my slide, I don't understand what 45/55 is                double scaleFactor = event_area / (828.03) * (44.8573/59.0878) ;
//priority check the chnnale maping by run the scrip
//----------




#include "TCanvas.h"
#include "TTree.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TChain.h"
#include "../include/mqROOTEvent.hh"
#include "../include/mqPMTRHit.hh"
//#include "/net/cms18/cms18r0/cms26r0/zheng/barsim/milliQanSim/include/mqROOTEvent.hh"
//#include "/net/cms18/cms18r0/cms26r0/zheng/barsim/milliQanSim/include/mqPMTRHit.hh"
//#include "milliQanSim/include/mqROOTEvent.hh"
//#include "milliQanSim/include/mqPMTRHit.hh"
#include "TGraph.h"
#include "TVector.h"
#include "TVectorD.h"
#include "TVectorF.h"
#include "TH1.h"
#include "TH1F.h"
#include "TString.h"
#include "TChain.h"
#include "TMultiGraph.h"
#include "TF1.h"
#include "TRandom3.h"
#include <vector>
#include <map>
#include "TSystem.h"
//R__LOAD_LIBRARY(/homes/tianjiad/milliQanSim/build/libBenchCore.so)
//R__LOAD_LIBRARY(/net/cms26/cms26r0/zheng/barsim/milliQanSim/build/libMilliQanCore.so)
R__LOAD_LIBRARY(../build/libMilliQanCore.dylib)
//R__LOAD_LIBRARY(../include/libMilliQanCore.so)

using namespace std;

//----------
#include <cassert>

//Returns 0 if the PMT is 18 or 19 (front)
//Returns 1 if the PMT is 20 or 21 (back)
int pmtPosition(int simChannel, int column, int row, int layer) {
    simChannel -= (4 * ((4 * column)+ (12 * row) + layer));
    assert((simChannel == 18 || simChannel == 19 || simChannel == 20 || simChannel == 21) && "Conversion done incorrectly.");
    return (simChannel == 18 || simChannel == 19) ? 0 : 1;
}


// Function to convert sim copy number to data copy number for Scint
int slabSimToDataScint(int simChannel) {
    // for slabs, we have a different mapping
    // given by 18+(4*y+12*z+i). i is the layer number, y is the column number, z is the row number
    // 0 is in the bottom right of the lowest layer, 8 is in the lower left, 36 is in the top right, 44 is in the top left
    // 1 is in the bottom right of the second layer, 9 is in the lower left, 37 is in the top right, 45 is in the top left
    // and so on
    // we can use this to map the sim channel to the data channel

    simChannel = simChannel - 18;
    int simlayer = simChannel % 4;
    int simrow = (simChannel - simlayer) / 12; // 0 at the bottom, 3 at the top (integer division truncates, same as floor for positive numbers)
    int simcolumn = ((simChannel - simlayer) % 12) / 4; // 0 at the right, 2 at the left

    int datalayer = simlayer; //do any manipulations between the mappings that is necessary
    int datarow = simrow;
    int datacolumn = simcolumn;

    return 12 * datalayer + 4 * datacolumn + datarow;
}

// Function to convert sim PMT branch copy number to data copy number for PMTs
int slabSimToDataPMT(int simChannel) {
    // The PMTs are numbered as (num)+4*(4*y+12*z+i), num=18,19,20, or 21
    // the slabs are numbered as 18+4*y+12*z+i
    // we want to convert the pmt number to the data slab number. So we subtract off the num first, then divide by 4
    simChannel = simChannel - 18;
    simChannel = simChannel - simChannel % 4;
    simChannel = simChannel / 4;
    return slabSimToDataScint(simChannel + 18);
}







void waveinject_slab(
    TString inputFile = "../build/beamMuonSlab_1kEvent.root",
    TString outputFile = "beamMuonSlab_1kEvent_waveinjected.root",
    TString waveformFile = "modified_waveform.root") {

    std::cout << "Injecting file " << inputFile << std::endl;
    std::cout << "Outputting file " << outputFile << std::endl;
    std::cout << "Using waveform template " << waveformFile << std::endl;

    TChain rootEvents("Events");
    rootEvents.Add(inputFile);
    mqROOTEvent* myROOTEvent = new mqROOTEvent();
    rootEvents.SetBranchAddress("ROOTEvent", &myROOTEvent);
    TFile* outfile = new TFile(outputFile, "RECREATE");

    const int nDigitizers = 6;
    const int nChannelsPerDigitizer = 16;
    const int nBins = 1024;
    double binWidth = 2.5;
    double rms_noise = 1;

    Float_t waveform[nDigitizers][nChannelsPerDigitizer][nBins] = {{{0}}};
    Double_t eventWeight = -1;

    TTree* injectedTree = new TTree("Events", "Tree with digitizer waveform data");
    injectedTree->Branch("waveform", waveform, Form("waveform[%d][%d][%d]/F", nDigitizers, nChannelsPerDigitizer, nBins));
    injectedTree->Branch("eventWeight", &eventWeight, "eventWeight/D");


    //might need to fix the the fit parameters for the slab?
    TF1 *fit = new TF1("fit", "gaus(0)", 0, 5000);
    fit->SetParameter(0, 7.23967e-02);
    fit->SetParameter(1, 1.48539e+03);
    fit->SetParameter(2, 2.90976e+02);


    //get the pulse SPE template
    TFile* f = new TFile(waveformFile);
    TH1F* pulse_shape = (TH1F*)f->Get("average_waveform");

    // Calibration array (scaled by dividing by 11)
    //might need to fix the calibration array for the slab
    //unlike bars PMT calibration, there are 4 slab PMT connections per slab scintillator, I don't know how to handle this yet.
    //Current solution is to make this all be 11.0 and matching the number of PMT from slab simulation(TODO: fix this)
    std::vector<double> cali = {1.74, 1.72, 1.91, 2.10, 2.05,1.97, 2.13, 1.99, 2.22, 1.99,2.06, 1.86, 1.95, 2.06, 2.34,2.20, 2.03, 2.13, 2.09, 2.16,2.07, 2.01, 1.94, 2.00, 1.95,2.09, 2.30, 2.24, 2.19, 2.18,2.23, 1.83, 1.86, 2.45, 2.10,2.00, 2.20, 2.11, 1.74, 1.72,1.73, 2.22, 2.22, 2.23, 2.11,2.00, 2.16, 2.70, 1.88, 1.78,2.02, 2.16, 1.49, 1.71, 1.58,1.78, 1.63, 2.18, 2.07, 2.09,1.85, 2.13, 1.88, 2.00, 1.67,1.44, 1.80, 1.93, 1.76, 2.07,1.87, 1.88, 1.92, 1.95, 2.15,2.31, 1.91, 2.48, 1.67, 2.28,1.77, 1.83, 1.91, 1.42, 1.74,2.33, 2.05, 1.80, 1.74, 2.11,1.96, 1.93, 2.28, 1.91, 2.07,1.88};
    for (auto& cal : cali) cal /= 10.00; // Divide each value by 10.0

    //max values for the slab digitizer pulse height investigation remains to be dealt with
    //change this to 96 entries of 1250(TODO: fix this)
    std::vector<double> maxValues = {
        1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250,
        1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250,
        1250, 1250, 1250, 1255, 1250, 1250, 1250, 1250, 1250, 1250,
        1250, 1250, 1250, 1250, 1250, 1250, 1250, 1255, 1255, 1250,
        1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250,
        1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250, 1250,
        1250, 1250, 1250, 1250,
        1250, 1250, 1250, 1250,
        1250, 1250, 1250, 1250,
        1250, 1250, 1250, 1250,
        1250, 1250, 1250, 1250,
        1250, 1250, 1250, 1250,
        1250, 1250, 1250, 1250,
        1250, 1250, 1250, 1250,
    };

    //Using Dariush's study from https://docs.google.com/presentation/d/1HEDyquZmwJO6FHQtmqME4q2W1po3xwkPLokhCiQpXkY/edit?slide=id.g34389da2a84_0_10#slide=id.g34389da2a84_0_10
    /*
    std::vector<double> SPE_height = {
        //channel 0-15
        27.4, 53.9, 52.1, 39.3, 27.7, 47.7, 57.9, 29.5, 64.8, 69.7, 29.1, 48.0, 31.4, 47.6, 45.1, 47.8,

        //channel 16 - 31
        54.2, 54.4, 62.2, 29.9, 39.0, 39.4, 34.7, 69.5, 39.4, 50.4, 52.7, 43.9, 61.5, 46.5, 61.7, 30.9,

        //channel 32 - 47
        49.5, 53.9, 38.6, 37.4, 47.9, 51.5, 44.3, 47.7, 54.3, 40.8, 68.2, 43.5, 56.3, 45.3, 61.6, 61.0,
        
        //channel 48 - 63
        55.7, 46.1, 35.5, 37.3, 53.5, 64.0, 47.9, 46.9, 33.1, 45.9, 52.0, 54.6, 36.1, 44.7, 34.9, 42.3,

        //64-79
        32.1, 44.0, 48.9, 44.7, 49.0, 46.4, 31.3,45.9, 27.9, 39.9, 46.2, 30.8, 29.2, 56.3, 29.6, 43.2,


        //80 -96
        28.1, 33.4, 30.4, 43.8, 54.4, 50.4, 73.9, 28.5, 47.7,38.8, 40.1, 37.8,31.9, 26.9, 56.1, 37.5
         
    };
    */
    //average summing amp SPE From Ryan's bench test
    std::vector<double> SPE_height = {
        65.5, 39, 38, 53.5, 33.5,65, 65, 37, 80.5, 72,63.5, 60.5, 55.5, 34, 69.5,57, 70.5, 67, 64.5, 35,47, 40, 45.5, 73.5, 47,59, 56, 45, 73.5, 52,68, 38.5, 66.5, 52, 41,47, 58, 56, 63.5, 59.5,68, 52, 69.5, 45.5, 66,53.5, 63, 54.5, 62.5, 59,45, 40, 65, 66.5, 56.5,54, 69.5, 58, 57.5, 60.5,52, 44.5, 36, 43.5, 45,35.5, 66, 53.5, 72, 49,36, 48.5, 33, 50, 64,41.5, 46, 53, 47, 39.5,52.5, 29, 49, 34.5, 63.5,58.5, 77, 47, 67.5, 45,47, 48.5, 38.5, 36, 46,78
    };
    TRandom3 randGen(2004);

    Long64_t nentries = rootEvents.GetEntries();
    std::cout << "Entries: " << nentries << std::endl;

    //for (Long64_t i = 0; i < nentries; i++) {
    for (Long64_t i = 0; i < 1000; i++) {
        if (i % (nentries / 100) == 0) std::cout << "Processing Event " << i << "..." << std::endl;
        rootEvents.GetEntry(i);
        memset(waveform, 0, sizeof(waveform));
        eventWeight = myROOTEvent->GetEventWeight();

        // Map to group PMT hits by PMT number
        std::map<int, std::vector<mqPMTRHit*>> pmtHitsMap;

        for (int j = 0; j < myROOTEvent->GetPMTRHits()->size(); j++) {
            mqPMTRHit* PMTRHit = myROOTEvent->GetPMTRHits()->at(j);
            int PMT_number = PMTRHit->GetPMTNumber();
            pmtHitsMap[PMT_number].push_back(PMTRHit);  //should I group the PMT that connects to the same Summing amp/digitizer together at here?yes assuming the gain are the same so far. But it may be difficult to consider case where the gains are difference.
        }
        
       
        for (const auto& pair : pmtHitsMap) {
            int PMT_number = pair.first;
            const std::vector<mqPMTRHit*>& hits = pair.second;
        //for validation: following the procedure from line 158-175 in the flatlightwithphotonslab.py
        //wish to covert sim PMT number offline channel number, then convert to digitizer channel number
        int ScintChannel = slabSimToDataPMT(PMT_number);  //get the scintillator number
        int channelType = pmtPosition(PMT_number, int(ScintChannel/4) % 3, int(ScintChannel%4), int(ScintChannel/12)); //get the front back pair in the offline summing amp channel.
        int remappedPMT = channelType + ScintChannel*2; //get the offline digitizer channel number
        int digitizer = remappedPMT / nChannelsPerDigitizer;  //digitizer number (raw)
        int channel = remappedPMT % nChannelsPerDigitizer;  //digitizer channel number (raw)
        //mapping debug
        //cout << "raw PMT_number" << PMT_number << endl;
        //cout << "ScintChannel" << ScintChannel << endl;
        //cout << "channelType" << channelType << endl;
        //cout << "remappedPMT" << remappedPMT << endl;
        //cout << "digitizer" << digitizer << endl;
        //cout << "digi channel" << channel << endl;  
   
        // Extract and sort hit times
        std::vector<double> hitTimes;
        for (const auto& hit : hits) {
        if(hit->GetFirstHitTime() > 500.0 ) continue;
    hitTimes.push_back(hit->GetFirstHitTime());
        }


	 if(hitTimes.size()==0) continue;
         std::sort(hitTimes.begin(), hitTimes.end());

         // Calculate median hit time
/*
	 double median_hit_time;
         size_t size = hitTimes.size();
         if (size % 2 == 0) {
            median_hit_time = 0.5 * (hitTimes[size / 2 - 1] + hitTimes[size / 2]);
         } else {
            median_hit_time = hitTimes[size / 2];
         }
*/
         double initial_hit_time = hitTimes[0];
	 double calibration = cali[remappedPMT]; //SPE height for the slab PMT
     //double calibration = 0.682; //Temporary fix for the calibration factor, should be changed to the actual calibration factor later.
     if (hits.size() > 5000) {
            double areaSum = 0.0;
            for (size_t k = 0; k < hits.size(); ++k) {
                if(randGen.Uniform() <= calibration) areaSum += fit->GetRandom();
            }

            TH1F* new_waveform = (TH1F*)pulse_shape->Clone();
            //new_waveform->Scale(areaSum * (1077.24 / 828.03) / new_waveform->Integral(480, 640));
            
	    // Scale only the bins within [500, 660]
            //double scaleFactor = areaSum / (828.03) * (44.8573/59.0878); //sample area, then scale by the pulse height for SPEs in data
            double scaleFactor = areaSum / (828.03) * (SPE_height[remappedPMT] / 59.0878);
            for (int bin = 500; bin <= 660; ++bin) {
                double binContent = new_waveform->GetBinContent(bin);
                new_waveform->SetBinContent(bin, binContent * scaleFactor);
            }

	    int integer_shift = static_cast<int>(initial_hit_time / binWidth);
            double fractional_shift = fmod(initial_hit_time, binWidth) / binWidth;

	    //uses median time rather than first for large hits
	    //int integer_shift = static_cast<int>(median_hit_time / binWidth);
            //double fractional_shift = fmod(median_hit_time, binWidth) / binWidth;

            for (int bin = 0; bin < nBins; ++bin) waveform[digitizer][channel][bin] = 0.0;

            for (int bin = 0; bin < nBins; ++bin) {
                int shifted_bin = bin - integer_shift;
                if (shifted_bin >= nBins || shifted_bin < 1) continue;

                double value = (1.0 - fractional_shift) * new_waveform->GetBinContent(shifted_bin) +
                               fractional_shift * new_waveform->GetBinContent(shifted_bin - 1);
                double noise = randGen.Gaus(0, rms_noise);
                waveform[digitizer][channel][bin] += (value + noise);

                if (waveform[digitizer][channel][bin] > maxValues[remappedPMT]) waveform[digitizer][channel][bin] = maxValues[remappedPMT];
                if (waveform[digitizer][channel][bin] < -50) waveform[digitizer][channel][bin] = -50;
            }

            for (int bin = 0; bin < integer_shift; ++bin) waveform[digitizer][channel][bin] = 0.0;

            delete new_waveform;
         } else {
            for (mqPMTRHit* PMTRHit : hits) {
               if(randGen.Uniform() > calibration) continue;
               double initial_hit_time = PMTRHit->GetFirstHitTime();
	       if(initial_hit_time>500) {continue;}
               TH1F* new_waveform = (TH1F*)pulse_shape->Clone();
               double event_area = fit->GetRandom();

               //new_waveform->Scale(event_area * (1077.24 / 828.03) / new_waveform->Integral(480, 640));
	       
	       // Scale only the bins within [500, 660]
               //double scaleFactor = event_area / (828.03) * (44.8573/59.0878);
               double scaleFactor = event_area  / (828.03) * (SPE_height[remappedPMT] / 59.0878);
               for (int bin = 500; bin <= 660; ++bin) {
                   double binContent = new_waveform->GetBinContent(bin);
                   new_waveform->SetBinContent(bin, binContent * scaleFactor);
               }

               int integer_shift = static_cast<int>(initial_hit_time / binWidth);
	       double fractional_shift = fmod(initial_hit_time, binWidth) / binWidth;

               //for (int bin = 0; bin < nBins; ++bin) waveform[digitizer][channel][bin] = 0.0;

               //for (int bin = 500; bin < 660; ++bin) {
               for (int bin = 0; bin < nBins; ++bin) {
                  int shifted_bin = bin - integer_shift;
                  if (shifted_bin >= nBins || shifted_bin < 1) continue;

                  double value = (1.0 - fractional_shift) * new_waveform->GetBinContent(shifted_bin) +
                                 fractional_shift * new_waveform->GetBinContent(shifted_bin - 1);
		  waveform[digitizer][channel][bin] += (value);

               }

               for (int bin = 0; bin < integer_shift; ++bin) waveform[digitizer][channel][bin] = 0.0;

               delete new_waveform;
            }
                  for (int bin = 0; bin < nBins; bin++) {
		        if(bin<500 || bin>660) waveform[digitizer][channel][bin]=0;
                        double noise = randGen.Gaus(0, rms_noise);
			waveform[digitizer][channel][bin] += noise;
                  	if (waveform[digitizer][channel][bin] > maxValues[remappedPMT]) waveform[digitizer][channel][bin] = maxValues[remappedPMT];
                  	if (waveform[digitizer][channel][bin] < -50) waveform[digitizer][channel][bin] = -50;
		  }

	   }
      }
      injectedTree->Fill();
   }

      cout << "final" << endl;
   outfile->cd();
   injectedTree->Write();
   f->Close();
   outfile->Close();
}

