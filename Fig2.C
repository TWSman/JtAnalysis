#include "style.C"
TGraphErrors *CalcRatio(TGraphErrors* num, TGraphErrors* den);
TGraphErrors *CalcRatio(TGraphErrors* num, TF1* den);
Int_t N=0;
void Fig2(){
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
    auto h = cfig->GetPad(1)->DrawFrame(0.9e-1,2e-3,3,2e5);
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



    auto * f = TFile::Open("RootFiles/Fig2.root");
    vector<int> I = {5};
    vector<int> S = { 24, 25, 27, 28 };
    vector<int> C = { kBlack,kRed,kBlue,kGreen+3,kMagenta+2,kPink+10,kRed,kBlue,kGreen+3,kMagenta+2,kPink+10,kBlack,kRed,kBlue};
    vector<int> L = {60};
    vector<int> R = {80};

    int n = 0;

    TLegend *leg = new TLegend(0.216792,0.634051,0.669173,0.811282, NULL, "brNDC");
    leg->SetFillColor(0);
    leg->SetTextAlign(12);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.0475034);
    leg->SetTextFont(42);
    TF1* gauss = nullptr;
    TF1* gamma = nullptr;
    for (auto i : I){
        auto g = (TGraphErrors*) gROOT -> FindObject(Form("jTSignalJetPt_Stat0%d",i)); 
        auto e = (TGraphErrors*) gROOT -> FindObject(Form("jTSignalJetPt_syst0%d",i)); 
        auto f = (TF1*) gROOT -> FindObject(Form("jTSignalJetPt_fit0%d",i)); 
        g -> SetMarkerColor(C.at(n));
        g -> SetLineColor(C.at(n));
        g -> SetMarkerStyle(S.at(n));
        g -> SetMarkerSize(2);
        //g -> GetXaxis() -> SetRangeUser(0.1,3);
        g -> Draw("PZ same e1");
        e -> SetMarkerColor(C.at(n));
        e -> SetLineColor(C.at(n));
        e -> SetMarkerStyle(S.at(n));
        e -> SetFillStyle(1001);
        e -> SetMarkerSize(2);
        e -> SetFillColorAlpha(C.at(n),0.5);
        e -> Draw("PZ same e2");
        leg -> AddEntry(e,"Data","lpe2f");
        gauss = new TF1("gausonly","gausn(0)",0,10);
        gauss->SetParameter(0,f->GetParameter(0));
        gauss->SetParameter(1,f->GetParameter(1));
        gauss->SetParameter(2,f->GetParameter(2));

        gamma = new TF1("gausonly","[0]* (pow([1],[2])/TMath::Gamma([2]))*exp(-[1]/x) * pow(x, -[2]-1)",0,10);
        gamma->SetParameter(0,f->GetParameter(3));
        gamma->SetParameter(1,f->GetParameter(4));
        gamma->SetParameter(2,f->GetParameter(5));
        f -> SetLineColor(C.at(n));
        f -> SetLineWidth(5);
        f -> Draw("same");
        gauss -> SetLineColor(C.at(n+1));
        gamma -> SetLineColor(C.at(n+2));
        gauss -> SetLineWidth(5);
        gamma -> SetLineWidth(5);
        gauss -> SetLineStyle(2);
        gamma -> SetLineStyle(9);
        gauss -> Draw("same");
        gamma -> Draw("same");

        leg -> AddEntry(f,"Total","l");
        leg -> AddEntry(gauss,"Narrow","l");
        leg -> AddEntry(gamma,"Wide","l");
        n++;
    }
    leg->Draw();

    leg = new TLegend(0.510025,0.647179,0.892231,0.942564, "", "brNDC");
    leg->SetFillColor(0);
    leg->SetTextAlign(12);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.045584);
    leg->AddEntry((TObject *)NULL, "p-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV", "");
    leg->AddEntry((TObject *)NULL, "Full jets", "");
    leg->AddEntry((TObject *)NULL, "Anti-#it{k}_{T}, R=0.4", "");
    leg->AddEntry((TObject *)NULL, "|#it{#eta}_{jet}|<0.25", "");
    leg -> AddEntry((TObject *)NULL,Form("%d < #it{p}_{T, jet} < %d GeV/#it{c}",L.at(0),R.at(0)),"");
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
    gPad -> SetGridy(true);
    gPad -> SetGridx(true);
    style();
    h = cfig->DrawFrame(0.9e-1,0.4,3,1.7);
    h->SetXTitle("#it{j}_{T} (GeV/#it{c})"); //h2->GetXaxis()->SetTitleSize(0.07); h2->GetXaxis()->SetTitleOffset(1);h2->GetXaxis()->SetLabelSize(0.05);
    // Please be consistent on the y label
    h->SetYTitle("Ratio (data/fit)"); //h2->GetYaxis()->SetTitleSize(0.07); h2->GetYaxis()->SetTitleOffset(1);h2->GetYaxis()->SetLabelSize(0.05);

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
    for (auto i : I){
        auto g = (TGraphErrors*) gROOT -> FindObject(Form("jTSignalJetPt_Stat0%d",i)); 
        auto e = (TGraphErrors*) gROOT -> FindObject(Form("jTSignalJetPt_syst0%d",i));
        auto f = (TF1*) gROOT -> FindObject(Form("jTSignalJetPt_fit0%d",i)); 
        if (n ==0 ) {
            d = (TGraphErrors*) g -> Clone();
            de = (TGraphErrors*) e -> Clone();
        }
        auto r = CalcRatio(g,f);
        r->SetMarkerColor(C.at(n));
        r -> SetLineColor(C.at(n));
        r -> SetMarkerStyle(S.at(n));
        //g -> GetXaxis() -> SetRangeUser(0.1,3);
        r -> Draw("PZ same e1");
        auto re = CalcRatio(e,f);
        re -> SetMarkerColor(C.at(n));
        re -> SetLineColor(C.at(n));
        re -> SetMarkerStyle(S.at(n));
        re -> SetFillStyle(1001);
        re -> SetFillColorAlpha(C.at(n),0.5);
        re -> Draw("PZ same e3");
        TLine * line = new TLine (0.9e-1,1,3,1);
        line -> SetLineColor(1);
        line -> SetLineStyle(5);
        line -> SetLineWidth(5);
        line -> Draw("same");

        n++;
    }
    cfig->GetPad(0)->SaveAs("../paperoriginal/newfigures/JtSignalFinalFitJetPt5.pdf");
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
TGraphErrors *CalcRatio(TGraphErrors* num, TF1* den)
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
    yden = den -> Eval(x);
    if (yden !=0 ) ratio -> SetPoint(i,x,y/yden);

    double nume = num -> GetErrorY(i)/y;

    ratio -> SetPointError (i,num->GetErrorX(i),nume*y/yden);
  }
  return ratio;
}

