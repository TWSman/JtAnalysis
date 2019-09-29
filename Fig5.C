#include "style.C"
#include "TGraphErrors.h"
TGraphErrors *CalcRatio(TGraphErrors* num, TGraphErrors* den);
TGraphErrors *CalcRatio(TGraphErrors* num, TF1* den);
TGraphErrors *CalcUnit(TGraphErrors* num);
void modelfit (TString file, TString directory, Int_t reference, TGraphErrors & wide, TGraphErrors &narrow);
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
vector<int> LS2 = {3,5,7};
Int_t N=0;
void Fig5(){
    style();
    TCanvas *cfig = new TCanvas("cfig", "cfig", 800, 1000);
    cfig->Divide(1, 3);
    cfig->GetPad(1)->SetFillColor(0);
    cfig->GetPad(2)->SetFillColor(0);
    cfig->GetPad(1)->SetTickx(1);
    cfig->GetPad(2)->SetTickx(1);
    cfig->GetPad(1)->SetTicky(0);
    cfig->GetPad(2)->SetTicky(0);
    cfig->GetPad(1)->SetPad(0.0, 0.5, 1.0, 1.0);
    cfig->GetPad(2)->SetPad(0.0, 0.3, 1.0, 0.5);
    cfig->GetPad(3)->SetPad(0.0, 0.0, 1.0, 0.3);
    cfig->GetPad(1)->SetTopMargin(0.02);
    cfig->GetPad(2)->SetTopMargin(0.00);
    cfig->GetPad(3)->SetTopMargin(0.00);
    cfig->GetPad(1)->SetRightMargin(0.02);
    cfig->GetPad(2)->SetRightMargin(0.02);
    cfig->GetPad(3)->SetRightMargin(0.02);
    cfig->GetPad(2)->SetTopMargin(0.00);
    cfig->GetPad(1)->SetBottomMargin(0.00);
    cfig->GetPad(2)->SetBottomMargin(0.0);
    cfig->GetPad(3)->SetBottomMargin(0.315);
    cfig->cd(1);
    auto h = cfig->GetPad(1)->DrawFrame(30,-0.1,122,2.1);
    h->GetYaxis()->SetNdivisions(1003);
    h -> GetYaxis()->SetTickLength(0.01);
    h->GetYaxis()->SetTitleSize(0.08);
    h->GetYaxis()->SetTitleOffset(0.93);
    h->GetYaxis()->SetLabelSize(0.08);
    h->GetYaxis()->SetLabelOffset(0.013);
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

    TLegend *leg = new TLegend(0.469925,0.645538,0.744361,0.898256, NULL, "brNDC");
    leg->SetFillColor(0);
    leg->SetTextAlign(12);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.0410256);
    leg->SetTextFont(42);
    TF1* gauss = nullptr;
    TF1* gamma = nullptr;
    gStyle->SetErrorX(0);

    TString file = "rootFiles/legotrain_CF_pPb_2305_20190109_LHC13bcde_minimal.root";
    TString dir = "AliJJetJtTask_kEMCEJE/AliJJetJtHistManager";
    TGraphErrors wide, narrow;
    modelfit(file, dir, 1, wide, narrow);
    TGraphErrors wide2, narrow2;
    modelfit(file, dir, 2,wide2,narrow2);
    leg->AddEntry(&wide,"Wide, Jet-axis reference","lf");
    leg->AddEntry(&wide2,"Wide, Leading track reference","lf");
    leg->AddEntry(&narrow,"Narrow, Jet-axis reference","lf");
    leg->AddEntry(&narrow2,"Narrow, Leading track reference","lf");
    //leg->AddEntry(trkreference,"Leading track reference","lf");
    leg -> Draw();

    leg = new TLegend(0.182957,0.630769,0.387218,0.858872, "", "brNDC");
    leg->SetFillColor(0);
    leg->SetTextAlign(12);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.045584);
    leg->AddEntry((TObject *)NULL, "#sqrt{#it{s}} = 5.02 TeV", "");
    leg->AddEntry((TObject *)NULL, "Full jets", "");
    leg->AddEntry((TObject *)NULL, "Anti-#it{k}_{T}, R=0.4", "");
    leg->Draw();

    leg = new TLegend(0.186717,0.862154,0.354637,0.95241, "", "brNDC");
    leg->SetFillColor(0);
    leg->SetTextAlign(12);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.065641);
    leg->AddEntry((TObject *)NULL, "PYTHIA 6", "");
    leg->Draw();


    
    cfig->GetPad(2)->cd();
    style();
    h = cfig->DrawFrame(30, 0.4, 122, 1.7);
    h->SetXTitle("#it{p}_{T, jet} (GeV/#it{c})"); //h2->GetXaxis()->SetTitleSize(0.07); h2->GetXaxis()->SetTitleOffset(1);h2->GetXaxis()->SetLabelSize(0.05);
    // Please be consistent on the y label
    h->SetYTitle("Ratio (W)"); //h2->GetYaxis()->SetTitleSize(0.07); h2->GetYaxis()->SetTitleOffset(1);h2->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.25);
    h->GetYaxis()->SetTitleOffset(0.32);
    h->GetYaxis()->SetLabelSize(0.18);
    h->GetYaxis()->SetLabelOffset(0.01);
    h->GetYaxis()->SetLabelFont(42);
    h->GetXaxis()->SetTitleSize(0.1);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.09);
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->CenterTitle(1);
    h->GetYaxis()->CenterTitle(1);
    h->Draw();
    gPad->SetGridy(1);

    auto ratio = CalcRatioNumErrorOnly (&wide2,&wide);
    auto ratio2 = CalcRatioNumErrorOnly (&narrow2,&narrow);

    ratio -> Draw("l3 same");
    auto center = CalcUnit(&wide);
    center -> Draw("3 same");
    leg = new TLegend(0.196742,0.771282,0.962406,0.968205, NULL, "brNDC");
    leg->SetFillColor(0);
    leg->SetTextAlign(12);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.0820513);
    leg->SetTextFont(42);
    leg->AddEntry(ratio,"Wide, leading track / jet axis","lf");
    leg->AddEntry(ratio2,"Narrow, leading track / jet axis","lf");
    //leg->Draw();
    
    cfig->GetPad(3)->cd();
    style();
    h = cfig->DrawFrame(30, 0.4, 122, 1.7);
    h->SetXTitle("#it{p}_{T, jet} (GeV/#it{c})"); //h2->GetXaxis()->SetTitleSize(0.07); h2->GetXaxis()->SetTitleOffset(1);h2->GetXaxis()->SetLabelSize(0.05);
    // Please be consistent on the y label
    h->SetYTitle("Ratio (N)"); //h2->GetYaxis()->SetTitleSize(0.07); h2->GetYaxis()->SetTitleOffset(1);h2->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.17);
    h->GetYaxis()->SetTitleOffset(0.46);
    h->GetYaxis()->SetLabelSize(0.12);
    h->GetYaxis()->SetLabelOffset(0.012);
    h->GetYaxis()->SetLabelFont(42);
    h->GetXaxis()->SetTitleSize(0.13);
    h->GetXaxis()->SetTitleOffset(1.03);
    h->GetXaxis()->SetLabelSize(0.12);
    h->GetXaxis()->SetLabelOffset(0.0);
    h->GetXaxis()->SetLabelFont(42);
    h->GetXaxis()->CenterTitle(1);
    h->GetYaxis()->CenterTitle(1);
    h->Draw();
    gPad->SetGridy(1);
    ratio2 -> Draw("l3 same");

    auto center2 = CalcUnit(&narrow);
    center2 -> Draw("3 same");

    cfig->GetPad(0)->SaveAs("../paperoriginal/newfigures/MixedFullJetsR04JetConeJtLeadingRefPtFrom2To6.pdf");
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



void modelfit (TString file, TString directory, Int_t reference, TGraphErrors & wide, TGraphErrors &narrow){

  TFile *f = TFile::Open(file.Data());
  vector<double> B3start = {7, 7, 5.06, 4.90, 10, 10, 7, 12.93, 12.93};
  vector<double> B5start = {1.5, 1.5, 4.5, 4.5, 1.5, 1.5, 2.88, 1.62, 1.62};
  vector<double> B4start = {3, 3, 9.91, 8.80, 3, 3, 5.62, 4.18, 4.18};
  Double_t x[5],ex[5];
  Double_t y[5],ey[5];
  Double_t y2[5],ey2[5];

  for (auto i = 3; i < 8; i++)
  {
    auto dir = (TDirectory *)f->GetDirectory(Form("%s/JetPtBin",directory.Data()));
    dir->cd();
    auto jetpt = (TH1D *)gROOT->FindObject(Form("JetPtBinNFin00JetPt0%d", i));
    double njets = jetpt->Integral();

    dir = (TDirectory *)f->GetDirectory(Form("%s/JetConeJtWeightBin",directory.Data()));
    if (reference == 2 )
      dir = (TDirectory *)f->GetDirectory(Form("%s/JetConeJtWeightLeadingRefBin", directory.Data()));
    dir->cd();
    auto incljt = (TH1D *)gROOT->FindObject(Form("JetConeJtWeightBinNFin00JetPt0%d", i));
    if (reference == 2) {
      incljt = (TH1D *)gROOT->FindObject(Form("JetConeJtWeightLeadingRefBinNFin00JetPt0%dXlong00", i));
      auto incljt2 = (TH1D *)gROOT->FindObject(Form("JetConeJtWeightLeadingRefBinNFin00JetPt0%dXlong01", i));
      auto incljt3 = (TH1D *)gROOT->FindObject(Form("JetConeJtWeightLeadingRefBinNFin00JetPt0%dXlong02", i));
      auto incljt4 = (TH1D *)gROOT->FindObject(Form("JetConeJtWeightLeadingRefBinNFin00JetPt0%dXlong03", i));
      incljt -> Add(incljt2);
      incljt -> Add(incljt3);
      incljt -> Add(incljt4);
    }
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

    x[i - 3] = jetptcenter.at(i);
    y[i - 3] = gaussigma;
    y2[i - 3] = gammaRMS;
    ex[i - 3] = 0;
    ey[i-3] = gaussigmae;
    ey2[i-3] = gammaRMSe;
  }
  auto model = new TGraphErrors(5, x, y,ex,ey);
  auto model2 = new TGraphErrors(5, x, y2,ex,ey2);

  model->SetLineColor(C.at(globalcount+1));
  model2->SetLineColor(C.at(globalcount+3));
  model2->SetFillColorAlpha(C.at(globalcount+3),0.2);
  model2->SetFillStyle(1001);
  model->SetFillColorAlpha(C.at(globalcount+1),0.2);
  model->SetFillStyle(1001);
  model->SetLineWidth(5);
  model2->SetLineWidth(5);
  model->SetLineStyle(LS.at(globalcount));
  model2->SetLineStyle(LS2.at(globalcount));
  model->Draw("le3 same");
  model2->Draw("le3 same");
  globalcount++;
  narrow = *model;
  wide = *model2;
}

