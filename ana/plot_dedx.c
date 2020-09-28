#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>


int plot_dedx()
{

    gStyle->SetOptStat(000);
    //gStyle->SetOptStat(111111);
        
    std::string condition;    // tag name to be used for fig.
    double pion_dedx[8], pion_dedx_sigma[8], pion_nhit[8];
    double kaon_dedx[8], kaon_dedx_sigma[8], kaon_nhit[8];    
    double gaussian_height;
        
    // 1. Pion dedx
    std::ifstream ifs1( "table_pion.txt" );
    if(!ifs1.is_open())
    {
        std::cout << "File open error!" << std::endl;
        return 0;
    }    
    for( int i=0; i<8; i++)
    {
        ifs1 >> condition >> gaussian_height >> pion_dedx[i] >> pion_dedx_sigma[i] >> pion_nhit[i] ;
    }
    ifs1.close();        

    // 2. Kaon dedx
    std::ifstream ifs2( "table_kaon.txt" );
    if(!ifs2.is_open())
    {
        std::cout << "File open error!" << std::endl;
        return 0;
    }    
    for( int i=0; i<8; i++)
    {
        ifs2 >> condition >> gaussian_height >> kaon_dedx[i] >> kaon_dedx_sigma[i] >> kaon_nhit[i] ;
    }
    ifs2.close();            


    // Momentum vector
    double momentum[8] = { 0.50, 0.75, 1.0, 3.0, 5.0, 10.0, 50.0, 100.0};


    // Plot I.   dE/dx 
    TCanvas *MyC1 = new TCanvas("MyC1","",800,600);

    TH1F *h_dedx = new TH1F("h_dedx","",1200,0,120.0);
    h_dedx->Draw();
    h_dedx->SetMinimum(0.0);
    h_dedx->SetMaximum(10.0);

    int xmin1 = h_dedx->GetXaxis()->FindBin(0.4);
    int xmax1 = h_dedx->GetXaxis()->FindBin(105);   
    h_dedx->GetXaxis()->SetRange(xmin1,xmax1);

    gPad->SetLogx(1);   
    //gPad->SetLogy(1);
           
    MyC1->SetGrid();
    
    h_dedx->SetName("");
    h_dedx->SetTitle("");
    
    h_dedx->GetXaxis()->SetTitle("Momentum [GeV/c]");
    h_dedx->GetYaxis()->SetTitle("dE/dx per unit length [keV/cm]");
    h_dedx->GetXaxis()->SetTitleFont(22);
    h_dedx->GetYaxis()->SetTitleFont(22);
    h_dedx->GetXaxis()->SetTitleOffset(1.2);
    h_dedx->GetYaxis()->SetTitleOffset(1.1);
    h_dedx->GetXaxis()->SetTitleSize(0.045);
    h_dedx->GetYaxis()->SetTitleSize(0.045);
    h_dedx->GetXaxis()->SetLabelFont(22);
    h_dedx->GetYaxis()->SetLabelFont(22);
    
    h_dedx->GetXaxis()->CenterTitle();
    h_dedx->GetYaxis()->CenterTitle();   
    
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.13);      
    

    double ex[8] = { 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0, 0.0, 0.0 };
        
    TGraph *gr1 = new TGraphErrors(8, momentum, pion_dedx, ex, pion_dedx_sigma);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerColor(2);
    gr1->SetMarkerSize(1.25);   
    gr1->SetLineWidth(1);
    gr1->SetLineStyle(2);      
    gr1->SetLineColor(2);   
    gr1->Draw("SAME,LP");


    TGraph *gr2 = new TGraphErrors(8, momentum, kaon_dedx, ex, kaon_dedx_sigma);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(1);
    gr2->SetMarkerSize(1.25);   
    gr2->SetLineWidth(1);
    gr2->SetLineStyle(2);      
    gr2->SetLineColor(1);   
    gr2->Draw("SAME,LP");

    double h_xmin_dedx = 0.4;
    double h_xmax_dedx = 110.0;
    
    double h_ymin_dedx = h_dedx->GetMinimum();
    double h_ymax_dedx = h_dedx->GetMaximum();         

    double line_dedx_x_ratio1 = 0.73;
    double line_dedx_x_ratio2 = 0.78;
    double line_dedx_x_ratio3 = 0.70;   // for text position     
       
    double line_dedx_x1 = std::pow(10.0,  line_dedx_x_ratio1 * ( std::log10( h_xmax_dedx ) - std::log10( h_xmin_dedx ) ) + std::log10( h_xmin_dedx ) );
    double line_dedx_x2 = std::pow(10.0,  line_dedx_x_ratio2 * ( std::log10( h_xmax_dedx ) - std::log10( h_xmin_dedx ) ) + std::log10( h_xmin_dedx ) );
    double line_dedx_x3 = std::pow(10.0,  line_dedx_x_ratio3 * ( std::log10( h_xmax_dedx ) - std::log10( h_xmin_dedx ) ) + std::log10( h_xmin_dedx ) );             
    
    double line_dedx_y_ratio1 = 0.85;
    double line_dedx_y_ratio2 = 0.75;
    double line_dedx_y_ratio3 = 0.80;
    double line_dedx_y_ratio4 = 0.75;   

    double line_dedx_y1 = line_dedx_y_ratio1 * h_ymax_dedx;
    double line_dedx_y2 = line_dedx_y_ratio2 * h_ymax_dedx;
    double line_dedx_y3 = line_dedx_y_ratio3 * h_ymax_dedx;
    double line_dedx_y4 = line_dedx_y_ratio4 * h_ymax_dedx;
       
    double label_x[1] = {line_dedx_x1};
    double label_y1[1] = {line_dedx_y1};
    double label_y2[1] = {line_dedx_y2};
    double label_y3[1] = {line_dedx_y3};
    double label_y4[1] = {line_dedx_y4};
       
    TGraph *label1 = new TGraph(1, label_x, label_y1);
    label1->SetMarkerStyle(20);
    label1->SetMarkerColor(1);
    label1->SetMarkerSize(1.25);      
    label1->SetLineWidth(1);
    label1->SetLineStyle(2);      
    label1->SetLineColor(1);      
    label1->Draw("SAME,LP");
    
    
    TGraph *label2 = new TGraph(1, label_x, label_y2);
    label2->SetMarkerStyle(20);
    label2->SetMarkerColor(2);
    label2->SetMarkerSize(1.25);      
    label2->SetLineWidth(1);
    label2->SetLineStyle(2);      
    label2->SetLineColor(1);      
    label2->Draw("SAME,LP");


    TLatex t_line1;
    t_line1.SetTextSize(0.04);
    t_line1.SetTextFont(22);
    t_line1.SetTextAlign(12);
    t_line1.DrawLatex( line_dedx_x2, line_dedx_y1, "K");
    

    TLatex t_line2;
    t_line2.SetTextSize(0.04);
    t_line2.SetTextFont(22);
    t_line2.SetTextAlign(12);
    t_line2.DrawLatex( line_dedx_x2, line_dedx_y2, "#pi");

    
       
    MyC1->Print("./fig/dedx_plot/dedx.gif");


    // Plot II.   dE/dx Resolution
    TCanvas *MyC2 = new TCanvas("MyC2","",800,600);

    TH1F *h_dedx_res = new TH1F("h_dedx_res","",1200,0,120.0);
    h_dedx_res->Draw();
    h_dedx_res->SetMinimum(0.0);
    h_dedx_res->SetMaximum(10.0);

    //int xmin1 = h_dedx->GetXaxis()->FindBin(0.3);
    //int xmax1 = h_dedx->GetXaxis()->FindBin(105);   
    h_dedx_res->GetXaxis()->SetRange(xmin1,xmax1);

    gPad->SetLogx(1);   
    //gPad->SetLogy(1);
           
    MyC2->SetGrid();
    
    h_dedx_res->SetName("");
    h_dedx_res->SetTitle("");
    
    h_dedx_res->GetXaxis()->SetTitle("Momentum [GeV/c]");
    h_dedx_res->GetYaxis()->SetTitle("dE/dx resolution [%]");
    h_dedx_res->GetXaxis()->SetTitleFont(22);
    h_dedx_res->GetYaxis()->SetTitleFont(22);
    h_dedx_res->GetXaxis()->SetTitleOffset(1.2);
    h_dedx_res->GetYaxis()->SetTitleOffset(1.1);
    h_dedx_res->GetXaxis()->SetTitleSize(0.045);
    h_dedx_res->GetYaxis()->SetTitleSize(0.045);
    h_dedx_res->GetXaxis()->SetLabelFont(22);
    h_dedx_res->GetYaxis()->SetLabelFont(22);
    
    h_dedx_res->GetXaxis()->CenterTitle();
    h_dedx_res->GetYaxis()->CenterTitle();   
    
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.13);      
    

    double pion_dedx_res[8];
    double kaon_dedx_res[8];    

    for( int i=0; i<8; i++)
    {
        pion_dedx_res[i] = pion_dedx_sigma[i] / pion_dedx[i] * 100.0 ;
        kaon_dedx_res[i] = kaon_dedx_sigma[i] / kaon_dedx[i] * 100.0 ;        
    }
    
        
    TGraph *gr3 = new TGraph(8, momentum, pion_dedx_res);
    gr3->SetMarkerStyle(20);
    gr3->SetMarkerColor(2);
    gr3->SetMarkerSize(1.25);   
    gr3->SetLineWidth(1);
    gr3->SetLineStyle(2);      
    gr3->SetLineColor(2);   
    gr3->Draw("SAME,LP");


    TGraph *gr4 = new TGraph(8, momentum, kaon_dedx_res); 
    gr4->SetMarkerStyle(20);
    gr4->SetMarkerColor(1);
    gr4->SetMarkerSize(1.25);   
    gr4->SetLineWidth(1);
    gr4->SetLineStyle(2);      
    gr4->SetLineColor(1);   
    gr4->Draw("SAME,LP");


    TGraph *label11 = new TGraph(1, label_x, label_y1);
    label11->SetMarkerStyle(20);
    label11->SetMarkerColor(1);
    label11->SetMarkerSize(1.25);      
    label11->SetLineWidth(1);
    label11->SetLineStyle(2);      
    label11->SetLineColor(1);      
    label11->Draw("SAME,LP");
    
    
    TGraph *label12 = new TGraph(1, label_x, label_y2);
    label12->SetMarkerStyle(20);
    label12->SetMarkerColor(2);
    label12->SetMarkerSize(1.25);      
    label12->SetLineWidth(1);
    label12->SetLineStyle(2);      
    label12->SetLineColor(1);      
    label12->Draw("SAME,LP");


    TLatex t_line11;
    t_line11.SetTextSize(0.04);
    t_line11.SetTextFont(22);
    t_line11.SetTextAlign(12);
    t_line11.DrawLatex( line_dedx_x2, line_dedx_y1, "K");
    

    TLatex t_line12;
    t_line12.SetTextSize(0.04);
    t_line12.SetTextFont(22);
    t_line12.SetTextAlign(12);
    t_line12.DrawLatex( line_dedx_x2, line_dedx_y2, "#pi");

    
    MyC2->Print("./fig/dedx_plot/dedx_res.gif");

    

    // Plot III.   <S>
    TCanvas *MyC3 = new TCanvas("MyC3","",800,600);

    TH1F *h_dedx_separation = new TH1F("h_dedx_separation","",1200,0,120.0);
    h_dedx_separation->Draw();
    h_dedx_separation->SetMinimum(0.0);
    h_dedx_separation->SetMaximum(10.0);

    //int xmin1 = h_dedx->GetXaxis()->FindBin(0.3);
    //int xmax1 = h_dedx->GetXaxis()->FindBin(105);   
    h_dedx_separation->GetXaxis()->SetRange(xmin1,xmax1);

    gPad->SetLogx(1);   
    //gPad->SetLogy(1);
           
    MyC3->SetGrid();
    
    h_dedx_separation->SetName("");
    h_dedx_separation->SetTitle("");
    
    h_dedx_separation->GetXaxis()->SetTitle("Momentum [GeV/c]");
    h_dedx_separation->GetYaxis()->SetTitle("S-value");
    h_dedx_separation->GetXaxis()->SetTitleFont(22);
    h_dedx_separation->GetYaxis()->SetTitleFont(22);
    h_dedx_separation->GetXaxis()->SetTitleOffset(1.2);
    h_dedx_separation->GetYaxis()->SetTitleOffset(1.1);
    h_dedx_separation->GetXaxis()->SetTitleSize(0.045);
    h_dedx_separation->GetYaxis()->SetTitleSize(0.045);
    h_dedx_separation->GetXaxis()->SetLabelFont(22);
    h_dedx_separation->GetYaxis()->SetLabelFont(22);
    
    h_dedx_separation->GetXaxis()->CenterTitle();
    h_dedx_separation->GetYaxis()->CenterTitle();   
    
    gPad->SetLeftMargin(0.13);
    gPad->SetBottomMargin(0.13);      
    

    double S_value[8];
    for( int i=0; i<8; i++)
    {
        double term1 = std::fabs( pion_dedx[i] - kaon_dedx[i] );
        double term2 = std::sqrt( pion_dedx_sigma[i]*pion_dedx_sigma[i] + kaon_dedx_sigma[i]*kaon_dedx_sigma[i] );

        S_value[i] = term1/term2;
    }            


    TGraph *gr5 = new TGraph(8, momentum, S_value); 
    gr5->SetMarkerStyle(20);
    gr5->SetMarkerColor(2);
    gr5->SetMarkerSize(1.25);   
    gr5->SetLineWidth(1);
    gr5->SetLineStyle(2);      
    gr5->SetLineColor(2);   
    gr5->Draw("SAME,LP");    

    
    MyC3->Print("./fig/dedx_plot/s_value.gif");

    
    return 0;
}
      
