{//BOF
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetStatColor(0);
  gStyle->SetTitleH(0.15);
  gStyle->SetTitleW(0.2);
  gStyle->SetPadColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  gStyle->SetOptTitle(0);
  gStyle->SetPadColor(10);
  gStyle->SetPadBorderMode(0);
  gStyle->SetDrawBorder(0);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.08);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetTitleFont(42,"XYZ");
  // gStyle->SetTitleSize(0.035,"YZ");
  gStyle->SetTitleSize(0.065,"X");
  //gStyle->SetTitleOffset(1.8,"YZ");
  gStyle->SetTitleOffset(1.2,"X");
  gStyle->SetLabelSize(0.07,"XYZ");
  gStyle->SetLabelFont(42,"XYZ");
  //gStyle->SetLabelOffset(0.01,"YZ");
  gStyle->SetLabelOffset(0.01,"X");
  //gStyle->SetOptLogx(0);
  //gStyle->SetOptLogy(1); 


  TCanvas* canv = new TCanvas("canv","my canvas",2000,1200);
  canv->SetFillColor(10);
  TPad *YnamePad = new TPad("pad1","pad1",0.0,0.0,0.1,0.99);
  YnamePad->SetFillColor(10);
  YnamePad->Draw();
  TPad *pad2 = new TPad("pad2","pad2",0.1,0.0,0.99,0.99);
  pad2->SetFillColor(10);
  pad2->Divide(4,4);	
  pad2->Draw();
  cout<<"ok"<<endl;
  TString hname;
  TH2F *hFrame[16];
  for(Int_t i=0;i<16;i++){
    hname=Form("hist%d",i+1);
    hFrame[i] = new TH2F(hname,hname,100,0.,50,100,-22,1);
    hFrame[i]->GetXaxis()->SetTitle("A");
    hFrame[i]->GetXaxis()->CenterTitle();
     //hFrame[i]->GetYaxis()->SetTickLength(0.022+i*0.003);
      //hFrame[i]->GetYaxis()->SetLabelSize(0.07-i*0.012);
      //hFrame[i]->GetYaxis()->SetLabelOffset(0.01+i*0.001);
    hFrame[i]->GetXaxis()->SetNdivisions(510);
    hFrame[i]->GetYaxis()->SetNdivisions(505);
      //hFrame[i]->GetXaxis()->SetTickLength(0.045-i*0.005);
  }
   
  TString title;
  TLatex Yname; 
  YnamePad->cd();
  title="E/A/MeV";
  Yname.SetTextSize(0.2);
  Yname.SetTextAlign(42);
  Yname.SetTextAngle(90.);
  Yname.DrawLatex(0.8,0.5,title.Data());
  
  ifstream infile[10];
  string strfile;

  //nndc_data
  const Int_t maxA1=61;
  const Int_t maxZ1=27;
  Double_t BE_ZA1[maxZ1][maxA1],X1[maxZ1][maxA1],xx1[maxZ1][maxA1],yy1[maxZ1][maxA1];
  Int_t Z1,A1;
  Double_t pp11,pp21,pp31;
  Double_t npoint1[maxZ1];
  
  for(Int_t i=0;i<maxZ1;i++){
    npoint1[i]=0;
    for(Int_t j=0;j<maxA1;j++){
      BE_ZA1[i][j]=0;
      X1[i][j]=0;
    }
  }
  const int cdata1=453;
  infile[1].open("Binding_energy_nndc.dat");
  // getline(infile[0],strfile);
  
  for(Int_t ii=0;ii<cdata1;ii++){
    infile[1]>>Z1>>A1>>pp11>>pp21;
    X1[Z1][A1]=A1;    
    BE_ZA1[Z1][A1]=pp11/A1;
    npoint1[Z1]++;
    //cout<<"Z="<<Z1<<" A="<< A1<<" pp1="<<pp11/A1<<"  "<<npoint1[Z1]<<endl;
  }
    
  for(Int_t i=0;i<maxZ1;i++){
    int nn1=0;
    for(Int_t j=0;j<maxA1;j++){
      if(BE_ZA1[i][j]!=0 && X1[i][j]!=0){
	xx1[i][nn1]=X1[i][j];
	yy1[i][nn1]=BE_ZA1[i][j];
	nn1++;
	//cout<<X1[i][j]<<endl;
      }
    }
  }
    
  TGraph* gr_ads1[27];
  gr_ads1[0] = new TGraph(1,xx1[0],yy1[0]);
  for(int iz1=1;iz1<16;iz1++){
    //cout<< npoint[iz]<<endl;
    gr_ads1[iz1] = new TGraph(npoint1[iz1],xx1[iz1],yy1[iz1]);
    gr_ads1[iz1]->SetMarkerStyle(20);
    gr_ads1[iz1]->SetMarkerColor(1);
    gr_ads1[iz1]->SetMarkerSize(1);
    gr_ads1[iz1]->SetLineStyle(1);
    gr_ads1[iz1]->SetLineColor(1);
    gr_ads1[iz1]->SetLineWidth(1);
    // gr_ads[iz]->Draw("LP");
  } 
  
  //35Ca-t1000-b=1_5wan
  // const Int_t maxA=56;
  //const Int_t maxZ=27;
  const Int_t maxA=61;
  const Int_t maxZ=30;
  TString Zname[maxZ]={"n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","X","XX","XXX"};
  Double_t BE_ZA[maxZ][maxA],X[maxZ][maxA],xx[maxZ][maxA],yy[maxZ][maxA];
  Int_t Z,A;
  Double_t pp1,pp2,pp3;
  Double_t npoint[maxZ];
  
  for(Int_t i=0;i<maxZ;i++){
    npoint[i]=0;
    for(Int_t j=0;j<maxA;j++){
      BE_ZA[i][j]=0;
      X[i][j]=0;
    }
  }
  
  //const int cdata=344;
  const int cdata=256;
  
  infile[0].open("BE_CoMDold_35Ca_t3000_b04.dat");//5th
  for(Int_t i=0;i<cdata;i++){
    infile[0]>>Z>>A>>pp1>>pp2;
    X[Z][A]=A;    
    BE_ZA[Z][A]=pp1/A;//for LowBE
    npoint[Z]++;
   }
  
  for(Int_t i=0;i<maxZ;i++){
    int nn=0;
    for(Int_t j=0;j<maxA;j++){
      if(BE_ZA[i][j]!=0 && X[i][j]!=0){
	xx[i][nn]=X[i][j];
	yy[i][nn]=BE_ZA[i][j];
	nn++;
	//cout<<X[i][j]<<endl;
      }
    }
  }
    
  TGraph* gr_ads[27];
  for(int iz=0;iz<16;iz++){
    gr_ads[iz] = new TGraph(npoint[iz],xx[iz],yy[iz]);
    gr_ads[iz]->SetMarkerStyle(24);
    gr_ads[iz]->SetMarkerColor(2);
    gr_ads[iz]->SetMarkerSize(1);
    gr_ads[iz]->SetLineStyle(1);
    gr_ads[iz]->SetLineColor(2);
    gr_ads[iz]->SetLineWidth(1);
    // gr_ads[iz]->Draw("LP");
  } 


  const Int_t maxA3=61;
  const Int_t maxZ3=30;
  Double_t BE_ZA3[maxZ3][maxA3],X3[maxZ3][maxA3],xx3[maxZ3][maxA3],yy3[maxZ3][maxA3];
  Int_t Z3,A3;
  Double_t pp13,pp23,pp33;
  Double_t npoint3[maxZ3];
  
  for(Int_t i=0;i<maxZ3;i++){
    npoint3[i]=0;
    for(Int_t j=0;j<maxA3;j++){
      BE_ZA3[i][j]=0;
      X3[i][j]=0;
    }
  }
  //const int cdata3=256;
  const int cdata3=187;
  //infile[3].open("BE_CoMDold_35Ca_t3000_b04_Hist.dat");//5th
  //infile[3].open("BE_CoMDoldMC_35Ca_t3000_b04.dat");//5th,added mass correction
  infile[3].open("BE_CoMDoldMassCor_35Ca_t3000_b04.dat");//5th,added mass correction
  for(Int_t i=0;i<cdata3;i++){
    infile[3]>>Z3>>A3>>pp13>>pp23;
    X3[Z3][A3]=A3;    
    BE_ZA3[Z3][A3]=pp13/A3;//for LowBE
    //BE_ZA3[Z3][A3]=pp13;//for HistBE
    npoint3[Z3]++;
   }
  
  for(Int_t i=0;i<maxZ3;i++){
    int nn3=0;
    for(Int_t j=0;j<maxA3;j++){
      if(BE_ZA3[i][j]!=0 && X3[i][j]!=0){
	xx3[i][nn3]=X3[i][j];
	yy3[i][nn3]=BE_ZA3[i][j];
	nn3++;
      }
    }
  }
    
  TGraph* gr_ads3[27];
  for(int iz=0;iz<16;iz++){
    gr_ads3[iz] = new TGraph(npoint3[iz],xx3[iz],yy3[iz]);
    gr_ads3[iz]->SetMarkerStyle(26);
    gr_ads3[iz]->SetMarkerColor(4);
    gr_ads3[iz]->SetMarkerSize(1);
    gr_ads3[iz]->SetLineStyle(1);
    gr_ads3[iz]->SetLineColor(4);
    gr_ads3[iz]->SetLineWidth(1);
    // gr_ads[iz]->Draw("LP");
  } 

  
 

  TLegend* legend=new TLegend(0.25,0.4,0.45,0.8);
  legend->SetTextFont(42);
  legend->SetFillColor(10);
  legend->SetTextSize(0.05);
  legend->SetBorderSize(0);
  legend->AddEntry(gr_ads1[1],"Exp.(NNDC)","PL");
  legend->AddEntry(gr_ads[1],"CoMD_original: 40Ca+40Ca_b=0-4fm,t=3000fm/c,LowBE","PL");
  legend->AddEntry(gr_ads3[1],"CoMD_original: 40Ca+40Ca_b=0-4fm,t=3000fm/c,LowBE_MassCrrection","PL");

  
  //pad2->cd(1);
  //legend->Draw();

  TLatex Chargename; 

  for(Int_t j=0;j<16;j++){
    pad2->cd(j+1);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.05);
    hFrame[j]->Draw();  
    gr_ads1[j]->Draw("PL");
    gr_ads[j]->Draw("PL");
    gr_ads3[j]->Draw("PL");
         
    Chargename.SetTextSize(0.1);
    Chargename.SetTextAlign(42);
    Chargename.DrawLatex(10,-3,Zname[j]);
    if(j==0)legend->Draw();
  }
 
}//EOF
