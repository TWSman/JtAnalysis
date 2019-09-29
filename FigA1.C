#include "style.C"
Int_t N=0;
void FigA1(){
    style();
    gStyle->SetErrorX(0);
    TCanvas *cfig = new TCanvas("cfig", "cfig", 800, 1000);
    cfig->Divide(1, 2);
    cfig->GetPad(1)->SetFillColor(0);
    cfig->GetPad(2)->SetFillColor(0);
    cfig->GetPad(1)->SetTickx(1);
    cfig->GetPad(2)->SetTickx(1);
    cfig->GetPad(1)->SetTicky(0);
    cfig->GetPad(2)->SetTicky(0);
    cfig->GetPad(1)->SetPad(0.0, 0.375, 1.0, 1.0);
    cfig->GetPad(2)->SetPad(0.0, 0.0, 1.0, 0.375);
    cfig->GetPad(1)->SetTopMargin(0.02);
    cfig->GetPad(1)->SetRightMargin(0.02);
    cfig->GetPad(2)->SetRightMargin(0.02);
    cfig->GetPad(2)->SetTopMargin(0.00);
    cfig->GetPad(1)->SetBottomMargin(0.00);
    cfig->GetPad(2)->SetBottomMargin(0.25);
    cfig->cd(1);
    auto h = cfig->GetPad(1)->DrawFrame(0.85e-1,0.9e-3,3,2e8);
    h -> GetYaxis()->SetNdivisions(1003);
    h -> GetYaxis()->SetTickLength(0.01);
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(1.25);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetYaxis()->SetLabelOffset(-0.001);
    h->GetYaxis()->CenterTitle(true);
    // Please be consistent on the y label
    h->SetYTitle("#frac{1}{#it{N}_{jets}} #frac{1}{#it{j}_{T}} #frac{d#it{N}}{d#it{j}_{T}} (#it{c}^{ 2}/GeV^{2})");
    h->GetXaxis()->SetTitleSize(0.08);
    h->GetXaxis()->SetLabelSize(0.07);
    h->GetXaxis()->SetTitleSize(0.08);



    auto * f = TFile::Open("RootFiles/Fig6.root");
    vector<int> I = {3,4,5};
    vector<int> S = { 24, 25, 27, 28 };
    vector<int> C = { kBlack,kRed,kBlue,kGreen+3,kMagenta+2,kPink+10,kRed,kBlue,kGreen+3,kMagenta+2,kPink+10,kBlack,kRed,kBlue};
    vector<int> L = {100,80,60,40};
    vector<int> R = {150,100,80,60};

    int n = 0;

    TLegend *leg = new TLegend(0.189223,0.0471386,0.506266,0.296004, NULL, "brNDC");
    leg->SetFillColor(0);
    leg->SetTextAlign(12);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->SetTextFont(42);
    for (auto i : I){
        auto g = (TH1D*) gROOT -> FindObject(Form("jTSignalJetPt05_R0%d",i)) -> Clone(Form("g%d",i)); 
        g -> Scale (TMath::Power(10,n));
        g -> SetMarkerColor(C.at(n));
        g -> SetLineColor(C.at(n));
        g -> SetMarkerStyle(S.at(n));
        g -> SetMarkerSize(2);
        g -> GetXaxis()->SetRangeUser(0.85e-1,3);
        //g -> GetXaxis() -> SetRangeUser(0.1,3);
        g -> Draw("PZ same e1");
        leg -> AddEntry(g,Form("#it{R} = 0.%d (#times10^{%d}) ",i,n),"lpe2f");
        n++;
    }
    leg->Draw();

    leg = new TLegend(0.582707,0.676718,0.964912,0.904821, "", "brNDC");
    leg->SetFillColor(0);
    leg->SetTextAlign(12);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.045584);
    leg->AddEntry((TObject *)NULL, "pp #sqrt{#it{s}} = 5.02 TeV", "");
    leg->AddEntry((TObject *)NULL, "Full jets", "");
    leg->AddEntry((TObject *)NULL, "Anti-#it{k}_{T}", "");
    leg->AddEntry((TObject *)NULL, "60<#it{p}_{T, jet}<80 GeV/#it{c}", "");
    leg->Draw();

    leg = new TLegend(0.185464,0.853949,0.353383,0.944205, "", "brNDC");
    leg->SetFillColor(0);
    leg->SetTextAlign(12);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.065641);
    leg->AddEntry((TObject *)NULL, "PYTHIA 8", "");
    leg->Draw();

    cfig->GetPad(2)->cd();
    gPad -> SetLogx(true);
    gPad -> SetGridy(true);
    gPad -> SetGridx(true);
    style();
    h = cfig->DrawFrame(0.85e-1,0.4,3,4.8);
    h->SetXTitle("#it{j}_{T} (GeV/#it{c})"); //h2->GetXaxis()->SetTitleSize(0.07); h2->GetXaxis()->SetTitleOffset(1);h2->GetXaxis()->SetLabelSize(0.05);
    // Please be consistent on the y label
    //h->SetYTitle("Ratio (/40 < #it{p}_{T} < 60 GeV/#it{c})"); //h2->GetYaxis()->SetTitleSize(0.07); h2->GetYaxis()->SetTitleOffset(1);h2->GetYaxis()->SetLabelSize(0.05);
    h->SetYTitle("Ratio"); //h2->GetYaxis()->SetTitleSize(0.07); h2->GetYaxis()->SetTitleOffset(1);h2->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.13);
    h->GetYaxis()->SetTitleOffset(0.63);
    h->GetYaxis()->SetLabelSize(0.1);
    h->GetYaxis()->SetLabelFont(42);
    h->GetXaxis()->SetTitleSize(0.1);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.09);
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->CenterTitle(1);
    h->GetYaxis()->CenterTitle(1);
    h->Draw();

    TH1D *d = nullptr;
    n = 0;
    I = {3,4,5};
    for (auto i : I){
        auto gg = (TH1D*) gROOT -> FindObject(Form("jTSignalJetPt05_R0%d",i)) -> Clone(Form("rr%d",i)); 
        gg -> GetXaxis()->SetRangeUser(0.85e-1,3);
        if (n ==0 ) {
            d = (TH1D*) gg -> Clone();
        }
        auto r = (TH1D*) gg -> Clone();
        r->Divide (d);
        r->SetMarkerColor(C.at(n));
        r -> SetLineColor(C.at(n));
        r -> SetMarkerStyle(S.at(n));
        r -> SetMarkerSize(2);
        //g -> GetXaxis() -> SetRangeUser(0.1,3);
        if (n>0) r -> Draw("PZ same e1");
        n++;
    }
    cfig->GetPad(0)->SaveAs("../paperoriginal/newfigures/RcomparisonJT.pdf");
}
