#include "style.C"
#include "TGraphErrors.h"
TGraphErrors *CalcRatio(TGraphErrors* num, TGraphErrors* den);
TGraphErrors *CalcRatio(TGraphErrors* num, TF1* den);
TGraphErrors *CalcUnit(TGraphErrors* num);
TGraphErrors* modelfit (TString file, TString directory);
TGraphErrors *CalcRatioNumErrorOnly(TGraphErrors* num, TGraphErrors* den);
TAxis AxisVar(TString name, std::vector<Double_t> bin)
{
  TAxis axis(bin.size() - 1, &bin.front());
  axis.SetName(name);
  return axis;
}
Int_t globalcount = 0;

vector<double> jetptcenter = {6.5, 12.45, 23.34, 33.83, 46.75, 67.73, 88.01, 116.11, 194.61};
vector<int> C = {kBlack, kRed, kBlue, kGreen + 3, kMagenta + 2, kPink + 10, kRed, kBlue, kGreen + 3, kMagenta + 2, kPink + 10, kBlack, kRed, kBlue};
vector<int> LS = {1,2,9};
Int_t N=0;
void Fig3(){
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
    auto h = cfig->GetPad(1)->DrawFrame(38,-0.1,122,2.1);
    h->GetYaxis()->SetNdivisions(1003);
    h -> GetYaxis()->SetTickLength(0.01);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(1.25);
    h->GetYaxis()->SetLabelSize(0.06);
    h->GetYaxis()->SetLabelOffset(0.008);
    h->GetYaxis()->CenterTitle(true);

    // Please be consistent on the y label
    h->SetYTitle("#sqrt{#LT #it{j}^{ 2}_{T} #GT} (GeV/#it{c})");
    h->GetXaxis()->SetTitleSize(0.08);
    h->GetXaxis()->SetLabelSize(0.07);
    h->GetXaxis()->SetTitleSize(0.08);

    auto * f = TFile::Open("RootFiles/jtSystematics.root");
    vector<int> I = {0};
    vector<int> S = { 20, 21 };
    vector<double> jetptbins = { 5, 10, 20, 30, 40, 60, 80, 100, 150, 500 };
    TAxis Jetptbins = AxisVar("Jetptbins",jetptbins);

    int n = 0;

    TLegend *leg = new TLegend(0.216792,0.540081,0.492481,0.792713, NULL, "brNDC");
    leg->SetFillColor(0);
    leg->SetTextAlign(12);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.0475034);
    leg->SetTextFont(42);
    TF1* gauss = nullptr;
    TF1* gamma = nullptr;
    gStyle->SetErrorX(0);
    auto g = (TGraphErrors *)gROOT->FindObject(Form("gGausRMS0%d", 0));
    auto ge = (TGraphErrors *)gROOT->FindObject(Form("gGausRMS0%d_Systematics", 0));
    auto ig = (TGraphErrors *)gROOT->FindObject(Form("gGammaRMS0%d", 0));
    auto ige = (TGraphErrors *)gROOT->FindObject(Form("gGammaRMS0%d_Systematics", 0));
    vector<TGraphErrors *> G{g, ge, ig, ige};
    for (int j = 0; j < ge->GetN(); j++)
    {
      ge->SetPointError(j, ge->GetErrorX(j), ge->GetErrorY(j));
    }
    int nn = 0;
    for (auto g : G)
    {
      for (int j = 0; j < g->GetN(); j++)
      {
        double x, y;
        g->GetPoint(j, x, y);
        if (x < 40)
          g->SetPoint(j, 0, 0);
        else
        {
          Int_t foundbin = Jetptbins.FindBin(x) - 1;
          g->SetPoint(j, jetptcenter.at(foundbin), y);
          double yerr = g->GetErrorY(j);
          double xerr = g->GetErrorX(j);
          if (nn == 0 || nn == 2)
            g->SetPointError(j, 0, yerr);
          if (nn == 1)
          {
            g->SetPointError(j, 2, yerr);
            if (y > 0 && yerr / y > 0.13)
              g->SetPointError(j, g->GetErrorX(j), 0.13 * y);
            if (y > 0 && yerr / y < 0.09)
              g->SetPointError(j, g->GetErrorX(j), 0.09 * y);
          }
          if (nn == 3)
          {
            g->SetPointError(j, 2, yerr);
            if (y > 0 && yerr / y > 0.13)
              g->SetPointError(j, g->GetErrorX(j), 0.13 * y);
            if (y > 0 && yerr / y < 0.10)
              g->SetPointError(j, g->GetErrorX(j), 0.08 * y);
          }
        }
      }
      nn++;
    }

    g->SetMarkerColor(C.at(n));
    g->SetLineColor(C.at(n));
    g->SetMarkerStyle(S.at(n));
    g->Draw("PZ same e1 ");

    ge->SetMarkerColor(C.at(n));
    ge->SetLineColor(C.at(n));
    ge->SetMarkerStyle(S.at(n));
    ge->SetFillColorAlpha(C.at(n), 0.2);
    ge->SetFillStyle(1001);
    ge->Draw("P same e2 ");

    ig->SetMarkerColor(C.at(n));
    ig->SetLineColor(C.at(n));
    ig->SetMarkerStyle(S.at(n + 1));
    ig->Draw("PZ same e1 ");

    ige->SetMarkerColor(C.at(n));
    ige->SetLineColor(C.at(n));
    ige->SetMarkerStyle(S.at(n + 1));
    ige->SetFillColorAlpha(C.at(n), 0.2);
    ige->SetFillStyle(1001);
    ige->Draw("PZ same e2");
    leg->AddEntry(ige, "Wide", "pef");
    leg->AddEntry(ge, "Narrow", "pef");
    n++;
    TString file = "RootFiles/Grid_Monash.root";
    TString dir = "JCDijetBaseTask/jcdijet";
    auto monash = modelfit(file, dir);
    leg->AddEntry(monash,"PYTHIA8 Monash","lf");


    file = "RootFiles/Grid_Tune4c.root";
    auto tune4c = modelfit(file, dir);
    leg->AddEntry(tune4c,"PYTHIA8 4C","lf");

    file = "RootFiles/Herwig7.root";
    dir = "/JJetJt";
    auto herwig = modelfit(file, dir);
    leg->AddEntry(herwig,"Herwig 7.0","lf");
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
    style();
    h = cfig->DrawFrame(38, 0.4, 122, 1.7);
    h->SetXTitle("#it{p}_{T, jet} (GeV/#it{c})"); //h2->GetXaxis()->SetTitleSize(0.07); h2->GetXaxis()->SetTitleOffset(1);h2->GetXaxis()->SetLabelSize(0.05);
    // Please be consistent on the y label
    h->SetYTitle("Ratio (wide)"); //h2->GetYaxis()->SetTitleSize(0.07); h2->GetYaxis()->SetTitleOffset(1);h2->GetYaxis()->SetLabelSize(0.05);
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
    gPad->SetGridy(1);
    gPad->SetGridx(1);

    auto rig = CalcUnit(ig);
    auto rige = CalcUnit(ige);
    rig -> Draw("P1 same");
    rige -> GetXaxis()->SetRangeUser(40,120);
    rige -> Draw("P2 same ");
    auto rmonash = CalcRatioNumErrorOnly (monash,ige);
    auto rtune4c = CalcRatioNumErrorOnly (tune4c,ige);
    auto rherwig = CalcRatioNumErrorOnly (herwig,ige);
    rmonash -> Draw("l3 same");
    rtune4c -> Draw("l3 same");
    rherwig -> Draw("l3 same");
    cfig->GetPad(0)->SaveAs("../paperoriginal/newfigures/RMSWithSystematics_Pythia.pdf");

    //cfig->GetPad(0)->SaveAs("../paperoriginal/newfigures/JtSignalFinalFitJetPt5.pdf");
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
  ratio->SetFillColorAlpha(num->GetFillColor(),0.2);
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
TGraphErrors *CalcRatioNumErrorOnly(TGraphErrors* num, TGraphErrors* den)
{
  TGraphErrors* ratio =  (TGraphErrors*) num -> Clone();
  ratio->SetName(Form("%s_%d",num->GetName(),N++));
  ratio->SetTitle(Form("%s_%d",num->GetName(),N++));
  ratio->SetLineColor(num->GetLineColor());
  ratio->SetLineStyle(num->GetLineStyle());
  ratio->SetLineWidth(num->GetLineWidth());
  ratio->SetMarkerStyle(num->GetMarkerStyle());
  ratio->SetFillColorAlpha(num->GetFillColor(),0.2);
  //ratio->SetFillStyle(0);
  //ratio->SetFillColor(0);
  //ratio->SetMarkerSize(0);

  for (Int_t i = 0; i < num->GetN(); i++)
  {
    double x, y;
    num->GetPoint(i, x, y);
    double xden, yden;
    for (Int_t j = 0; j < den->GetN(); j++)
    {
      double xtemp, ytemp;
      den->GetPoint(j, xtemp, ytemp);
      cout << "xtemp" << xtemp << " " << xden << endl;
      if (fabs(xtemp-x)<1) {xden = xtemp; yden = ytemp;};
    }
    cout<< "x, xden"<<x<<" "<<xden<<endl;
    if (yden != 0)
      ratio->SetPoint(i, x, y / yden);

    double nume = num->GetErrorY(i) / y;



    ratio->SetPointError(i, num->GetErrorX(i), nume * y / yden);
  }
  return ratio;
}

TGraphErrors *CalcRatio(TGraph* num, TGraphErrors* den)
{
  TGraphErrors* ratio =  (TGraphErrors*) num -> Clone();
  ratio->SetName(Form("%s_%d",num->GetName(),N++));
  ratio->SetTitle(Form("%s_%d",num->GetName(),N++));
  ratio->SetLineColor(num->GetLineColor());
  ratio->SetLineStyle(num->GetLineStyle());
  ratio->SetLineWidth(num->GetLineWidth());
  ratio->SetMarkerStyle(num->GetMarkerStyle());
  ratio->SetFillColorAlpha(num->GetFillColor(),0.2);
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
  ratio->SetFillColorAlpha(num->GetFillColor(),0.2);
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
TGraphErrors *CalcUnit(TGraphErrors* num)
{
  TGraphErrors* ratio =  (TGraphErrors*) num -> Clone();
  ratio->SetName(Form("%s_%d",num->GetName(),N++));
  ratio->SetTitle(Form("%s_%d",num->GetName(),N++));
  ratio->SetLineColor(num->GetLineColor());
  ratio->SetLineStyle(num->GetLineStyle());
  ratio->SetLineWidth(num->GetLineWidth());
  ratio->SetMarkerStyle(num->GetMarkerStyle());
  ratio->SetFillColorAlpha(num->GetFillColor(),0.2);
  //ratio->SetFillStyle(0);
  //ratio->SetFillColor(0);
  //ratio->SetMarkerSize(0);

  for (Int_t i = 0; i < num->GetN(); i++) {
    double x, y;
    num -> GetPoint(i,x,y);
    ratio -> SetPoint(i,x,1);

    double nume = num -> GetErrorY(i)/y;

    ratio -> SetPointError (i,num->GetErrorX(i),nume*1);
  }
  return ratio;
}



TGraphErrors* modelfit (TString file, TString directory){

  TFile *f = TFile::Open(file.Data());
  vector<double> B3start = {7, 7, 5.06, 4.90, 10, 10, 7, 12.93, 12.93};
  vector<double> B5start = {1.5, 1.5, 4.5, 4.5, 1.5, 1.5, 2.88, 1.62, 1.62};
  vector<double> B4start = {3, 3, 9.91, 8.80, 3, 3, 5.62, 4.18, 4.18};
  Double_t x[4],ex[4];
  Double_t y[4],ey[4];
  Double_t y2[4],ey2[4];

  for (auto i = 4; i < 8; i++)
  {
    auto dir = (TDirectory *)f->GetDirectory(Form("%s/JetPtBin",directory.Data()));
    dir->cd();
    auto jetpt = (TH1D *)gROOT->FindObject(Form("JetPtBinNFin00JetPt0%d", i));
    double njets = jetpt->Integral();

    dir = (TDirectory *)f->GetDirectory(Form("%s/JetConeJtWeightBin",directory.Data()));
    dir->cd();
    auto incljt = (TH1D *)gROOT->FindObject(Form("JetConeJtWeightBinNFin00JetPt0%d", i));
    incljt->Scale(1. / njets, "width");
    //incljt -> Draw("PZsame");

    dir = (TDirectory *)f->GetDirectory(Form("%s/BgTrkNumberBin",directory.Data()));
    dir->cd();
    auto hnormbg = (TH1D *)gROOT->FindObject(Form("BgTrkNumberBinNFin00JetPt0%d", i));
    double normbg = hnormbg->Integral();

    dir = (TDirectory *)f->GetDirectory(Form("%s/BgJtWeightBin",directory.Data()));
    dir->cd();
    auto hbg = (TH1D *)gROOT->FindObject(Form("BgJtWeightBinNFin00JetPt0%d", i));
    hbg->Scale(1. / normbg, "width");
    hbg->SetLineColor(2);
    //hbg -> Draw("PZsame");
    incljt->Add(hbg, -1);
    //incljt -> Draw("PZ same");

    //auto f = new TF1("fit"," gausn(0) + [3]* (pow([4],[5])/TMath::Gamma([5]))*exp(-[4]/x) * pow(x, -[5]-1)",0,10);
    auto gaussfit = new TF1("gaussfit", "gausn", 0, 10);
    gaussfit->FixParameter(1, 0);
    gaussfit->SetParLimits(2, 0.1, 0.5);
    gaussfit->SetParLimits(0, 40, 150);
    incljt->Fit("gaussfit", "QN");

    auto invG = new TF1("invG", "[0]* (pow([1],[2])/TMath::Gamma([2]))*exp(-[1]/x) * pow(x, -[2]-1)", 0, 10);
    invG->SetParameter(0, B3start.at(i));
    if (i < 6)
      invG->SetParLimits(0, 0, 10);
    else
      invG->SetParLimits(0, 0, 25);

    invG->SetParameter(1, B5start.at(i));
    invG->SetParLimits(1, 0.90, 4.5);

    invG->SetParameter(2, B4start.at(i));
    invG->SetParLimits(2, 2.0, 15);

    incljt->Fit("invG", "QN", "", 1, 3);

    auto gaussfit3 = new TF1("gaussfit3", "gausn(0) + [3]* (pow([4],[5])/TMath::Gamma([5]))*exp(-[4]/x) * pow(x, -[5]-1)", 0, 10);
    gaussfit3->SetParameter(0, gaussfit->GetParameter(0));
    gaussfit3->SetParameter(1, gaussfit->GetParameter(1));
    gaussfit3->SetParameter(2, gaussfit->GetParameter(2));
    gaussfit3->SetParameter(3, invG->GetParameter(0));
    gaussfit3->SetParameter(4, invG->GetParameter(1));
    gaussfit3->SetParameter(5, invG->GetParameter(2));
    gaussfit3->FixParameter(1, 0);

    gaussfit3->SetParLimits(0, 40, 150);
    gaussfit3->SetParLimits(2, 0.1, 0.5);
    if (i < 6)
      gaussfit3->SetParLimits(3, 0, 10);
    else
      gaussfit3->SetParLimits(3, 0, 25);
    gaussfit3->SetParLimits(4, 0.9, 4.5);
    gaussfit3->SetParLimits(5, 2.0, 15);
    auto lastbin = incljt->FindLastBinAbove(1e-5);
    auto end = incljt->GetBinCenter(lastbin);
    incljt->Fit("gaussfit3", "QN", "", 0.01, end);
    //gaussfit3->Draw("same");

    auto B2 = gaussfit3->GetParameter(0);
    auto B1 = gaussfit3->GetParameter(2);
    auto B2e = gaussfit3->GetParError(0);
    auto B1e = gaussfit3->GetParError(2);
    auto gaussigma = TMath::Sqrt(2) * B1;
    auto gausyield = B2 * B1 / TMath::Sqrt(2 * TMath::Pi());
    auto gausyielde = TMath::Sqrt((B2 * B2 * B1e * B1e + B1 * B1 * B2e * B2e) / (2 * TMath::Pi()));
    auto gaussigmae = TMath::Sqrt(2) * B1e;
    auto constant = gaussfit3->GetParameter(3);
    auto constante = gaussfit3->GetParError(3);
    auto alpha = gaussfit3->GetParameter(5);
    auto alphae = gaussfit3->GetParError(5);

    auto beta = gaussfit3->GetParameter(4);
    auto betae = gaussfit3->GetParError(4);
    auto peak = beta / (alpha + 1);
    auto peake = TMath::Sqrt(pow(betae / (alpha + 1), 2) + pow(alphae * beta / pow(alpha + 1, 2), 2));
    auto gammaRMS = beta / sqrt((alpha - 2) * (alpha - 3));
    auto gammaYield = constant * beta / (alpha - 1);
    auto gammaYielde = sqrt(pow(beta * constante / (alpha - 1), 2) + pow(constant * beta * alphae / pow(alpha - 1, 2), 2) + pow(constant * betae / (alpha - 1), 2));
    auto gammaRMSe = sqrt(pow((5 - 2 * alpha) * beta * alphae / pow(2 * ((alpha - 2) * (alpha - 3)), 1.5), 2) + pow(betae / sqrt((alpha - 2) * (alpha - 3)), 2));

    x[i - 4] = jetptcenter.at(i);
    cout<<"x4 :"<<x[i-4]<<endl;
    y[i - 4] = gaussigma;
    y2[i - 4] = gammaRMS;
    ex[i - 4] = 0;
    ey[i-4] = gaussigmae;
    ey2[i-4] = gammaRMSe;
  }
  auto model = new TGraphErrors(4, x, y,ex,ey);
  auto model2 = new TGraphErrors(4, x, y2,ex,ey2);

  model->SetLineColor(C.at(globalcount+1));
  model2->SetLineColor(C.at(globalcount+1));
  model2->SetFillColorAlpha(C.at(globalcount+1),0.2);
  model2->SetFillStyle(1001);
  model->SetFillColorAlpha(C.at(globalcount+1),0.2);
  model->SetFillStyle(1001);
  model->SetLineWidth(5);
  model2->SetLineWidth(5);
  model->SetLineStyle(LS.at(globalcount));
  model2->SetLineStyle(LS.at(globalcount));
  model->Draw("le3 same");
  model2->Draw("le3 same");
  globalcount++;
  return model2;
}

