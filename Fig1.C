#include "style.C"
TGraphErrors *CalcRatio(TGraphErrors* num, TGraphErrors* den);
Int_t N=0;
void Fig1(){
    style();
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
    auto h = cfig->GetPad(1)->DrawFrame(0.9e-1,0.9e-3,3,2e8);
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



    auto * f = TFile::Open("RootFiles/Fig1.root");
    vector<int> I = {3,2,1,0};
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
        auto g = (TGraphErrors*) gROOT -> FindObject(Form("jTSignalJetPt_Stat0%d",i)); 
        auto e = (TGraphErrors*) gROOT -> FindObject(Form("jTSignalJetPt_Syst0%d",i)); 
        for (int j = 0; j < g->GetN(); j++)
        {
            g->GetY()[j] *= TMath::Power(10,3-n);  //2.1;
            g->GetEY()[j] *= TMath::Power(10,3-n); //2.1;
            e->GetY()[j] *= TMath::Power(10,3-n);  //2.1;
            e->GetEY()[j] *= TMath::Power(10,3-n); //2.1;
        }
        g -> SetMarkerColor(C.at(n));
        g -> SetLineColor(C.at(n));
        g -> SetMarkerStyle(S.at(n));
        //g -> GetXaxis() -> SetRangeUser(0.1,3);
        g -> Draw("PZ same e1");
        e -> SetMarkerColor(C.at(n));
        e -> SetLineColor(C.at(n));
        e -> SetMarkerStyle(S.at(n));
        e -> SetFillStyle(1001);
        e -> SetFillColorAlpha(C.at(n),0.5);
        e -> Draw("PZ same e2");
        leg -> AddEntry(e,Form("%d < #it{p}_{T, jet} < %d GeV/#it{c} (#times10^{%d}) ",L.at(n),R.at(n),3-n),"lpe2f");
        n++;
    }
    leg->Draw();

    leg = new TLegend(0.511278,0.71434,0.892231,0.94226, "", "brNDC");
    leg->SetFillColor(0);
    leg->SetTextAlign(12);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.045584);
    leg->AddEntry((TObject *)NULL, "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "");
    leg->AddEntry((TObject *)NULL, "Full jets", "");
    leg->AddEntry((TObject *)NULL, "Anti-#it{k}_{T}, R=0.4", "");
    leg->AddEntry((TObject *)NULL, "|#it{#eta}_{jet}|<0.25", "");
    leg->Draw();

    leg = new TLegend(0.190476,0.82942,0.357143,0.920108, "", "brNDC");
    leg->SetFillColor(0);
    leg->SetTextAlign(12);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.107962);
    leg->AddEntry((TObject *)NULL, "ALICE", "");
    leg->Draw();

    cfig->GetPad(2)->cd();
    gPad -> SetLogx(true);
    gPad -> SetLogy(true);
    gPad -> SetGridy(true);
    gPad -> SetGridx(true);
    style();
    h = cfig->DrawFrame(0.9e-1,0.2,3,20);
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

    TGraphErrors *d;
    TGraphErrors *de;
    n = 0;
    I = {0,1,2,3};
    for (auto i : I){
        auto g = (TGraphErrors*) gROOT -> FindObject(Form("jTSignalJetPt_Stat0%d",i)); 
        auto e = (TGraphErrors*) gROOT -> FindObject(Form("jTSignalJetPt_Syst0%d",i));
        if (n ==0 ) {
            d = (TGraphErrors*) g -> Clone();
            de = (TGraphErrors*) e -> Clone();
        }
        auto r = CalcRatio(g,d);
        r->SetMarkerColor(C.at(3-n));
        r -> SetLineColor(C.at(3-n));
        r -> SetMarkerStyle(S.at(3-n));
        //g -> GetXaxis() -> SetRangeUser(0.1,3);
        if (n>0) r -> Draw("PZ same e1");
        auto re = CalcRatio(e,de);
        re -> SetMarkerColor(C.at(3-n));
        re -> SetLineColor(C.at(3-n));
        re -> SetMarkerStyle(S.at(3-n));
        re -> SetFillStyle(1001);
        re -> SetFillColorAlpha(C.at(3-n),0.2);
        if (n>0) re -> Draw("PZ same e3");

        n++;
    }
    cfig->GetPad(0)->SaveAs("../paperoriginal/newfigures/jTwithSystematics.pdf");
}
TGraphErrors *CalcRatio(TGraphErrors* num, TGraphErrors* den)
{
  TGraphErrors* ratio =  (TGraphErrors*) num -> Clone();
  ratio->SetName(Form("%s_%d",num->GetName(),N++));
  ratio->SetTitle(Form("%s_%d",num->GetName(),N++));
  ratio->SetLineColor(num->GetLineColor());
  ratio->SetLineStyle(num->GetLineStyle());
  ratio->SetLineWidth(num->GetLineWidth());
  ratio->SetMarkerStyle(num->GetMarkerStyle());
  ratio->SetFillColorAlpha(num->GetFillColor(),0.5);
  //ratio->SetFillStyle(0);
  //ratio->SetFillColor(0);
  //ratio->SetMarkerSize(0);

  for (Int_t i = 0; i < num->GetN(); i++) {
    double x, y;
    double xden, yden;
    num -> GetPoint(i,x,y);
    den -> GetPoint(i,xden, yden);
    if (yden !=0 ) ratio -> SetPoint(i,x,y/yden);

    double nume = num -> GetErrorY(i)/y;


    double dene = den -> GetErrorY(i)/yden;

    double esum = sqrt(nume*nume+dene*dene);

    ratio -> SetPointError (i,num->GetErrorX(i),esum*y/yden);
  }
  return ratio;
}
