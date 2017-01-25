#include <iostream>
#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TPaveStats.h"
#include "TSystem.h"
#include "TMath.h"
#include "TF1.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2D.h"
#include <string>
// basic file operations
#include <iostream>
#include <fstream>
#include "TStyle.h"
#include "TPad.h"
#include "TLatex.h"
#include "TPaveLabel.h"

void DrawStat(TH1* gr1)
{
  gr1->Draw();
  gPad->Update();
  gr1->SetName(gr1->GetTitle());
	TPaveStats *st = ((TPaveStats*)(gr1->GetListOfFunctions()->FindObject("stats")));
  if (st) 
	{
  	st->SetTextColor(gr1->GetLineColor());
    st->SetX1NDC(0.64); st->SetX2NDC(0.99);
    st->SetY1NDC(0.75); st->SetY2NDC(0.99);
  }
}

 // define a function with 3 parameters
Double_t fitf(Double_t *x,Double_t *par) 
{
  Double_t arg = 0.;
  if (par[2]!=0) arg  = (x[0] - par[1])/par[2];
  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;
}


const double length_strip=100;/*cm*/
const double velocity=2.0/3*29.9792;/*cm.ns-1*/

std::vector<TH1F*>DiffTimes;
TH1F* DiffTime=new TH1F("Diff_Time","Diff_Time",400,-20e-9,20e-9);
TH1F* triggerr=new TH1F("Diff_trigger","Diff_trigger",400,-20,20);
TH2F* Position=new TH2F("Position","Position",32,0,32,400,-200,200);
TH2F* StripCoincidence=new TH2F("StripCoincidence","StripCoincidence",32,0,32,1,0,50);
class hit
{
  public:
    hit()=delete;
    ~hit(){};
    hit(ULong64_t& abci,Double_t& ltd,UInt_t& bci,int& id,unsigned int& even,unsigned int& tim,unsigned& mezzanin,double& rtd,UInt_t& stri,double diff)
    {
      abcid=abci;
      ltdc=ltd;
      bcid=bci;
      idx=id;
      event=even;
      time=tim;
      mezzanine=mezzanin;
      rtdc=rtd;
      strip=stri;
      diff=diff_time;
      if(ChannelTDC1toStrip.size()>strip)
      {
          if(mezzanine==1)strip_1=ChannelTDC1toStrip[strip];
          else strip_1=ChannelTDC2toStrip[strip];
      }
    } 
    double diff_time;
    ULong64_t abcid;
    Double_t ltdc;
    UInt_t bcid;
    int idx;
    unsigned int event;
    unsigned int time;
    unsigned int mezzanine;
    UInt_t strip;
    double rtdc;
    UInt_t strip_1;
    static std::vector<int> ChannelTDC1toStrip;
    static std::vector<int>ChannelTDC2toStrip;
    static int number1;
    static int number2;
};

class trigger
{
  public:
    trigger()=delete;
    ~trigger(){};
    trigger(ULong64_t& abci,Double_t& ltd,UInt_t& bci,int& id,unsigned int& even,unsigned int& tim,unsigned& mezzanin,double& rtd)
    {
      if(mezzanin==1)trigger::number1++;
      else trigger::number2++;
      abcid=abci;
      ltdc=ltd;
      bcid=bci;
      idx=id;
      event=even;
      time=tim;
      mezzanine=mezzanin;
      rtdc=rtd;
    } 
    ULong64_t abcid;
    Double_t ltdc;
    UInt_t bcid;
    int idx;
    unsigned int event;
    unsigned int time;
    unsigned int mezzanine;
    double rtdc;
    static int number1;
    static int number2;
};

int trigger::number1=0;
int trigger::number2=0;

class Selector
{
  public:
    Selector()=delete;
    Selector(TTree* data):dataTree(data)
    {
      dataTree->SetBranchAddress("abcid", &abcid);
      dataTree->SetBranchAddress("ltdc", &ltdc);
      dataTree->SetBranchAddress("bcid", &bcid);
      dataTree->SetBranchAddress("idx", &idx);
      dataTree->SetBranchAddress("event", &event);
      dataTree->SetBranchAddress("time", &time);
      dataTree->SetBranchAddress("mezzanine", &mezzanine);
      dataTree->SetBranchAddress("strip", &strip);
      dataTree->SetBranchAddress("rtdc", &rtdc);
    };
    void Fill(unsigned int i)
    {
      dataTree->GetEntry(i);
      if(old_idx!=idx)
      {
        Analyse();
        hits.clear();
        triggers.clear();
        hitstriggered.clear();
      }
      if(strip==16)triggers.push_back(trigger(abcid,ltdc,bcid,idx,event,time,mezzanine,rtdc));
      else hits.push_back(hit(abcid,ltdc,bcid,idx,event,time,mezzanine,rtdc,strip,0));
      old_idx=idx;
    }
    ~Selector()
    {
    };
    
    void Analyse()
    {
      for(unsigned int i=0;i!=triggers.size();++i)
      {
        bool is_touch1=false;
        bool is_touch2=false;
        for(unsigned int j=0;j!=hits.size();++j)
        {
          if(hits[j].abcid==triggers[i].abcid&&fabs(hits[j].bcid-triggers[i].bcid)<=5) 
          {
            if(hits[j].mezzanine==1)is_touch1=true; 
            else is_touch2=true;
            hits[j].diff_time=hits[j].rtdc*1.0e-9-triggers[j].rtdc*1.e-9;
            hitstriggered.push_back(hits[j]);
          }
          //if(is_touch1==true)hit::number1++;
          //if(is_touch2==true)hit::number2++;
        }
	if(is_touch1==true)hit::number1++;
        if(is_touch2==true)hit::number2++;
      }
      Select();
    }
    void Select()
    {
      bool is_touch1=false;
      bool is_touch2=false;
      for(unsigned int i=0;i!=hitstriggered.size();++i)
      {
        for(unsigned int j=i;j!=hitstriggered.size();++j)
        {
          if(j==i)continue;
          if(hitstriggered[i].strip_1==hitstriggered[j].strip_1&&hitstriggered[i].strip==myfriend[hitstriggered[j].strip])
          {
            if(hitstriggered[i].mezzanine==1) is_touch1=true;
            else is_touch2=true;
              double dt=0.0;
              if(hitstriggered[j].strip%2==0) dt=(hitstriggered[i].rtdc-hitstriggered[j].rtdc)*1.e-9;
              else dt=(hitstriggered[j].rtdc-hitstriggered[i].rtdc)*1.e-9;;
              std::cout<<dt<<std::endl;
              Position->Fill(hitstriggered[j].strip_1,-100+velocity*0.5*dt);
              StripCoincidence->Fill(hitstriggered[j].strip_1,25);
              DiffTime->Fill(dt);
              DiffTimes[hitstriggered[j].strip_1]->Fill(dt);
              triggerr->Fill(100./velocity-(hitstriggered[i].rtdc+hitstriggered[j].rtdc)*1.e-9);
          }
        }
      }
	if(is_touch1==true)Selector::number1++;
        if(is_touch2==true)Selector::number2++;
    }
    ULong64_t abcid;
    Double_t ltdc;
    UInt_t bcid;
    int idx;
    int old_idx;
    unsigned int event;
    unsigned int time;
    unsigned int mezzanine;
    UInt_t strip;
    double rtdc;
    static std::vector<int>myfriend;
    TTree* dataTree;
    std::vector<hit> hits;
    std::vector<hit> hitstriggered;
    std::vector<trigger>triggers;
    static int number1;
    static int number2;
};

int hit::number1=0;
int hit::number2=0;
int Selector::number1=0;
int Selector::number2=0;
std::vector<int> hit::ChannelTDC1toStrip={14,14,12,12,8 ,8 ,6 ,6 ,4 ,4 ,2 ,2 ,0 ,0 };
std::vector<int> hit::ChannelTDC2toStrip={16,16,18,18,20,20,22,22,24,24,28,28,30,30};
std::vector<int> Selector::myfriend={1,0,3,2,5,4,7,6,9,8,11,10,13,12};
  
  
int main(int argc , char *argv[])
{
  gStyle->SetOptFit(1111);
  std::ofstream myfile;
  myfile.open ("Results_efficiency.txt",std::ios_base::app);
  for(unsigned int i=0;i!=31;++i)
  {
    DiffTimes.push_back(new TH1F(("Diff_Time_"+std::to_string(i)).c_str(),("Diff_Time_"+std::to_string(i)).c_str(),400,-20e-9,20e-9));
    DiffTimes[i]->SetTitle(("T_{1}-T_{2} for channel "+std::to_string(i)).c_str());
    DiffTimes[i]->GetXaxis()->SetTitle("T_{1}-T_{2} (ns)");
    DiffTimes[i]->GetYaxis()->SetTitle("#");
  }
  DiffTime->SetTitle("T_{1}-T_{2}");
  DiffTime->GetXaxis()->SetTitle("T_{1}-T_{2} (ns)");
  DiffTime->GetYaxis()->SetTitle("#");
  TStopwatch ti;
  ti.Start();
  std::string filename="";
  if(argc==2)
  {
		filename=argv[1];
  }
  else
  {
  	std::cout<<"Please give me a filename !!"<<std::endl;
  	std::exit(1);
  }
  TFile dataFile(filename.c_str());
  if (dataFile.IsOpen() != true) 
  {
  	std::cout <<"Impossible to read " <<filename<< std::endl;
    std::exit(1);
  }
  TTree *dataTree = (TTree *)dataFile.Get("events");
  if (!dataTree) 
  {
  	std::cout << "Impossible to read TTree RAWData"<< std::endl;
    delete dataTree;
    std::exit(1);
  }
  std::size_t found = filename.find_last_of("/\\");
  std::string namee=filename.substr(found+1);
  TFile f(("Result_"+namee).c_str(),"recreate");
  Selector selector(dataTree);
  for (unsigned int i = 0; i < dataTree->GetEntries(); i++) 
  {
    if(i%50000==0)std::cout<<"Entry : "<<i<<std::endl;
    selector.Fill(i);
  }
  dataFile.Close();
  //delete dataTree;
  for(unsigned int i=0;i!=DiffTimes.size();++i)
  {
    if(DiffTimes[i]->GetEntries()!=0)
    {
    	TF1 *func = new TF1("fit",fitf,DiffTimes[i]->GetMean()-DiffTimes[i]->GetRMS(),DiffTimes[i]->GetMean()+DiffTimes[i]->GetRMS(),3);
    	// set the parameters to the mean and RMS of the histogram
    	func->SetParameters(1.0,DiffTimes[i]->GetMean()-DiffTimes[i]->GetRMS(),DiffTimes[i]->GetRMS()/2.0);
        //func->SetParameters(5.0,DiffTimes[i]->GetMean(),DiffTimes[i]->GetRMS());
    	// give the parameters meaningful names
    	func->SetParNames ("Constant","Mean_value","Sigma");
    	// call TH1::Fit with the name of the TF1 object
    	DiffTimes[i]->Fit("fit");
      DrawStat(DiffTimes[i]);
      double sigma= func->GetParameter("Sigma")*1.0/sqrt(2);
      std::string text="#sigma = "+std::to_string(sigma);
      TPaveLabel *fLabel   = new TPaveLabel(0.14,0.70,0.36,0.80,text.c_str(),"NDC");
      fLabel->Draw();
    	DiffTimes[i]->Write();
			delete func;
   }
   delete DiffTimes[i];
  }
  DiffTimes.clear();
  TF1 *func = new TF1("fit",fitf,DiffTime->GetMean()-DiffTime->GetRMS(),DiffTime->GetMean()+DiffTime->GetRMS(),3);
  func->SetParameters(1.0,DiffTime->GetMean(),DiffTime->GetRMS(),1.0,DiffTime->GetMean(),DiffTime->GetRMS(),DiffTime->GetMean());
  //func->SetParameters(4.0,DiffTime->GetMean(),DiffTime->GetRMS());
  func->SetParNames ("Constant","Mean_value","Sigma");
  DiffTime->Fit("fit");
  DrawStat(DiffTime);
  StripCoincidence->Write();
  DiffTime->Write();
  Position->Write();
  triggerr->Write();
  f.Close();
  ti.Stop();
  delete func;
  std::cout<<"Efficiency : "<<hit::number1*100.0/trigger::number1<<"  "<<hit::number2*100.0/trigger::number2<<std::endl;
  std::cout<<"Efficiency both side (%) : "<<Selector::number1*100.0/hit::number1<<"  "<<Selector::number2*100.0/hit::number2<<std::endl;
  myfile <<namee<<" TDC1: "<<hit::number1*100.0/trigger::number1<<" ; "<<Selector::number1*100.0/hit::number1<<"; TDC2: "<<hit::number2*100.0/trigger::number2<<"  "<<Selector::number2*100.0/hit::number2<<std::endl;
  myfile.close();
  std::cout << " Time :" << ti.RealTime() << "  " << ti.CpuTime() << std::endl;
}

