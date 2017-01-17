#include <iostream>
#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TMath.h"
#include "TF1.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2D.h"
#include <string>

 // define a function with 3 parameters
Double_t fitf(Double_t *x,Double_t *par) 
{
  Double_t arg = 0;
  if (par[2]!=0) arg = (x[0] - par[1])/par[2];
  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;
}


const double length_strip=100;/*cm*/
const double velocity=2.0/3*29.9792;/*cm.ns-1*/
TH2F* Position=new TH2F("Position","Position",32,0,32,100,0,100);
TH2F* StripCoincidence=new TH2F("StripCoincidence","StripCoincidence",32,0,32,1,0,50);
std::vector<TH1F*>DiffTimes;
TH1F* DiffTime=new TH1F("Diff_Time","Diff_Time",400,-20e-9,20e-9);
TH1F* triggerr=new TH1F("Diff_trigger","Diff_trigger",400,-20,20);
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
      number++;
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
    static int number;
};

int trigger::number=0;

class Selector
{
  public:
    Selector()=delete;
    Selector(TTree* data,double time_trig1,double time_trig2):dataTree(data)
    {
      time_from_trigger={time_trig1,time_trig2};
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
        for(unsigned int j=0;j!=hits.size();++j)
        {
          if(hits[j].abcid==triggers[i].abcid&&fabs(hits[j].bcid-triggers[i].bcid)<=2) 
          {
            if(hits[j].mezzanine==1) hit::number1++;
            else hit::number2++;
            hits[j].diff_time=hits[j].rtdc*1.0e-9-triggers[j].rtdc*1.e-9;
            hitstriggered.push_back(hits[j]);
          }
        }
      }
      Select();
    }
    void Select()
    {
      for(unsigned int i=0;i!=hitstriggered.size();++i)
      {
        for(unsigned int j=i;j!=hitstriggered.size();++j)
        {
          if(j<=i)continue;
          if(hitstriggered[i].strip_1==hitstriggered[j].strip_1&&hitstriggered[i].strip==myfriend[hitstriggered[j].strip])
          {
            if(hitstriggered[i].mezzanine==1) Selector::number1++;
            else Selector::number2++;
           //if(velocity*(hitstriggered[i].diff_time+hitstriggered[j].diff_time)<105&&velocity*(hitstriggered[i].diff_time+hitstriggered[j].diff_time)>95)
           //{
              Position->Fill(hitstriggered[j].strip_1,50-velocity*0.5e-9*(hitstriggered[i].rtdc-hitstriggered[j].rtdc));
              StripCoincidence->Fill(hitstriggered[j].strip_1,25);
              DiffTime->Fill((hitstriggered[i].rtdc-hitstriggered[j].rtdc)*1.e-9);
              std::cout<<hitstriggered[j].strip_1<<"  "<<(hitstriggered[i].rtdc-hitstriggered[j].rtdc)*1.e-9<<std::endl;
              std::cout<<hitstriggered[i].diff_time+hitstriggered[j].diff_time<<"  "<<velocity*(hitstriggered[i].diff_time+hitstriggered[j].diff_time)<<std::endl<<std::endl;
              DiffTimes[hitstriggered[j].strip_1]->Fill((hitstriggered[i].rtdc-hitstriggered[j].rtdc)*1.e-9);
              std::cout<<100./velocity-(hitstriggered[i].rtdc+hitstriggered[j].rtdc)*1.e-9<<std::endl;
              triggerr->Fill(100./velocity-(hitstriggered[i].rtdc+hitstriggered[j].rtdc)*1.e-9);
           //}
          }
        }
      }
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
    std::array<double,2> time_from_trigger;
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
  for(unsigned int i=0;i!=31;++i)
  {
    DiffTimes.push_back(new TH1F(("Diff_Time_"+std::to_string(i)).c_str(),("Diff_Time_"+std::to_string(i)).c_str(),400,-20e-9,20e-9));
  }
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
  TFile f(("Result"+std::string(argv[1])).c_str(),"recreate");
  Selector selector(dataTree,100.,100.);
  for (unsigned int i = 0; i < dataTree->GetEntries(); i++) 
  {
    //std::cout<<"Entry : "<<i<<std::endl;
    selector.Fill(i);
  }
  dataFile.Close();
  //delete dataTree;
  for(unsigned int i=0;i!=DiffTimes.size();++i)
  {
    TF1 *func = new TF1("fit",fitf,DiffTimes[i]->GetMean()-DiffTimes[i]->GetRMS(),DiffTimes[i]->GetMean()+DiffTimes[i]->GetRMS(),3);
    // set the parameters to the mean and RMS of the histogram
    func->SetParameters(1.0,DiffTimes[i]->GetMean(),DiffTimes[i]->GetRMS());
    // give the parameters meaningful names
    func->SetParNames ("Constant","Mean_value","Sigma");
    // call TH1::Fit with the name of the TF1 object
    DiffTimes[i]->Fit("fit");
    DiffTimes[i]->Write();
  }
  TF1 *func = new TF1("fit",fitf,DiffTime->GetMean()-DiffTime->GetRMS(),DiffTime->GetMean()+DiffTime->GetRMS(),3);
  func->SetParameters(1.0,DiffTime->GetMean(),DiffTime->GetRMS());
  func->SetParNames ("Constant","Mean_value","Sigma");
  DiffTime->Fit("fit");
  StripCoincidence->Write();
  DiffTime->Write();
  Position->Write();
  triggerr->Write();
  f.Close();
  ti.Stop();
  std::cout<<"Efficiency : "<<hit::number1*1.0/trigger::number<<"  "<<hit::number2*1.0/trigger::number<<std::endl;
  std::cout<<"Efficiency : "<<Selector::number1*1.0/trigger::number<<"  "<<Selector::number2*1.0/trigger::number<<std::endl;
  std::cout << " Time :" << ti.RealTime() << "  " << ti.CpuTime() << std::endl;
}

