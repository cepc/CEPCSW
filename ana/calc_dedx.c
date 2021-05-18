#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

void get_mean( TH1D * h_energy_deposit, double& mean, int& nhit)
{
    TH1D* h_energy_deposit_wlengthCut_90 = new TH1D("h_energy_deposit_lcut_90", "h_energy_deposit_lcut_90", 200, 0.0, 20.0 );     

    // Calculate 90%
    int total_entries = h_energy_deposit->GetEntries();
    
    int count = 0;
    int max_bin = 0;
    for( int i=0; i<200; i++)
    {
        count += h_energy_deposit->GetBinContent(i+1);
        if( count > total_entries * 0.90 )
        {
            max_bin = i;
            break;
        }
    }

    int number_of_hit = 0;
    for( int i=0; i < max_bin ; i++)
    {
        h_energy_deposit_wlengthCut_90->SetBinContent(i+1,  h_energy_deposit->GetBinContent(i+1) );
        number_of_hit += h_energy_deposit->GetBinContent(i+1);
    }

    TCanvas *MyC1e = new TCanvas("MyC1e", "",600,600);    
    MyC1e->SetTicks(1,1);   

    h_energy_deposit_wlengthCut_90->Draw();


    mean = h_energy_deposit_wlengthCut_90->GetMean();
    nhit    = number_of_hit;

    delete h_energy_deposit_wlengthCut_90;
    delete MyC1e;
}

int calc_dedx()
{

    //gStyle->SetOptStat(000);
    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(111111);        

    // Read information about input/output filename
    std::ifstream ifs( "info_input.txt" );
    if(!ifs.is_open())
    {
        std::cout << "File open error!" << std::endl;
        return 0;
    }
    
    std::string rootfile;       // input rootfile
    std::string condition;    // tag name to be used for fig.
    std::string outputfile;   // output filename 
    for( int i=0; i<1; i++)
    {
        ifs >> rootfile >> condition >> outputfile ;
        std::cout << rootfile << std::endl;
        std::cout << condition << std::endl;
        std::cout << outputfile << std::endl;        
    }
    ifs.close();    
    


    TFile *f1 = new TFile( rootfile.c_str() );

    TTree* tree = dynamic_cast<TTree*>(gDirectory->Get("events"));

    // Get info. of number of entries in this tree
    int NEntry = tree->GetEntries();
    std::cout << std::endl;    
    std::cout << "Number of entries (events) in this tree : " <<  NEntry << std::endl;
    std::cout << std::endl;
    

    // Histograms
    TH1D* h_EDep = new TH1D("h_EDep", "h_EDep", 2000, 0.0, 20.0 );
    TH1I* h_Nhit = new TH1I("h_Nhit", "h_Nhit", 300, 0.0, 300.0 );
    TH1I* h_Nhit_all = new TH1I("h_Nhit_all", "h_Nhit_all", 300, 0.0, 300.0 );    
    
    
    // Only pick up MDC hit information
    int nhit;
    double mean;
    for( int i=0; i<NEntry; i++)        
    {
        TH1D* h_energy_deposit_wlengthCut = new TH1D("h_energy_deposit_wlengthCut", "h_energy_deposit_wlengthCut", 200, 0.0, 20.0 );

        // MDC1
        tree->Draw("DCHCollection.position.EDep * 1e6 / DCHCollection.pathLength * 10 >> h_energy_deposit_wlengthCut",
                   " abs(DCHCollection.pathLength-10) < 1 && DCHCollection.time > 1 && DCHCollection.time < 20.0", "", 1, i );   
        
        // MDC2
        tree->Draw("DCHCollection2.position.EDep * 1e6 / DCHCollection2.pathLength * 10 >> + h_energy_deposit_wlengthCut",
                   " abs(DCHCollection2.pathLength-10) < 1 && DCHCollection2.time > 1 && DCHCollection2.time < 20.0", "", 1, i );   
    
        get_mean(  h_energy_deposit_wlengthCut, mean, nhit );

        std::cout << "event number: " << i << std::endl;
        std::cout << "truncated mean = " << mean << std::endl;
        std::cout << "number of hits used for averaging : " << nhit << std::endl;
        std::cout << std::endl;
        
        delete h_energy_deposit_wlengthCut;

        // If number of hits per event is too low, it is discarded.
        if( nhit > 10 )
        {
            h_EDep->Fill( mean );
            h_Nhit->Fill( nhit );
        }
        h_Nhit_all->Fill( nhit );        
    }
    
    TCanvas *MyC_EDep = new TCanvas("MyC_EDep", "",600,600);    
    MyC_EDep->SetTicks(1,1);   

    h_EDep->Draw();
    
    // Fitting
    double edep_mean = h_EDep->GetMean();
    TF1 * fgaus = new TF1("fgaus", "gaus", 0, 20 );
    fgaus->SetParameter(1, edep_mean);
    fgaus->SetParLimits(1, 0.0, 20.0);     
    fgaus->SetParameter(2, edep_mean*0.05);

    fgaus->SetLineWidth(2);
    fgaus->SetLineStyle(2);

    // Set Fit Range
    double tmp_min_pre;
    double tmp_max_pre;
    if( condition == std::string("kaon_0.50GeV") )  // Multi-peaks due to (possibly) decays 
    {
        tmp_min_pre = 2.08;   // temporally set narrow range
        tmp_max_pre = 2.4;    //
    }
    else
    {
        tmp_min_pre = edep_mean - 1.5;    
        tmp_max_pre = edep_mean + 1.5;  
        if( tmp_min_pre < 0 )tmp_min_pre=0;        
    }

    int bin_min_pre = h_EDep->GetXaxis()->FindBin( tmp_min_pre );
    int bin_max_pre = h_EDep->GetXaxis()->FindBin( tmp_max_pre );
    h_EDep->GetXaxis()->SetRange(bin_min_pre, bin_max_pre);
    
    h_EDep->Fit("fgaus");
    double res[3] = {
        fgaus->GetParameter(0),
        fgaus->GetParameter(1),
        fgaus->GetParameter(2),
    };
    

    std::cout << "par[0] = " << res[0]  << " , par[1] = " << res[1] << " , par[2] = " << res[2] << std::endl;

    
    h_EDep->GetXaxis()->SetTitleFont(22);
    h_EDep->GetYaxis()->SetTitleFont(22);
    h_EDep->GetXaxis()->SetTitleOffset(1.0);
    //h_EDep->GetYaxis()->SetTitleOffset(1.6);
    h_EDep->GetYaxis()->SetTitleOffset(1.2);
    h_EDep->GetXaxis()->SetTitleSize(0.045);
    h_EDep->GetYaxis()->SetTitleSize(0.045);
    
    h_EDep->GetXaxis()->SetLabelSize(0.04);
    h_EDep->GetYaxis()->SetLabelSize(0.04);

    /*
    h_EDep->SetName("");
    //h_EDep->SetTitle("");
    std::string hist_title = condition;    
    h_EDep->SetTitle( hist_title.c_str() );    
    */
    
    h_EDep->GetXaxis()->SetTitle("Energy deposit [keV]");
    h_EDep->GetYaxis()->SetTitle("Number of events");
    
    h_EDep->GetXaxis()->CenterTitle();
    h_EDep->GetYaxis()->CenterTitle();
    
    h_EDep->GetXaxis()->SetNdivisions(505);
    h_EDep->GetYaxis()->SetNdivisions(505);    
  
    //gPad->SetBottomMargin(0.17);
    gPad->SetLeftMargin(0.14);

    std::string figname_edep = "./fig/fitting/edep_" + condition + ".gif" ;
    MyC_EDep->Print( figname_edep.c_str() );
      
    TCanvas *MyC_Nhit = new TCanvas("MyC_Nhit", "",600,600);    
    MyC_Nhit->SetTicks(1,1);   

    h_Nhit->Draw();    
    double nhit_mean = h_Nhit->GetMean();
    
    std::string figname_nhit= "./fig/fitting/nhit_" + condition + ".gif" ;
    MyC_Nhit->Print( figname_nhit.c_str() );


    // Print output parameters to a file
    std::ofstream ofs( outputfile.c_str() , ios::app );
    ofs << condition << " " << res[0] << " " << res[1] << " " << res[2] << " " << nhit_mean << std::endl;    
    ofs.close();

    
    return 0;
}
