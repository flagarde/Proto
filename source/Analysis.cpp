//g++ -std=c++11 -I $ROOTSYS/include -L/home/couet/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGpad -lTree -lRint -lMatrix -lMathCore main.cpp -o test

#include <iostream>
#include <fstream>
#include <string>
#include <cstdint>
#include <array>
#include <sstream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2D.h"
#include <unistd.h>
#include <thread>
#include "TStreamerInfo.h"
#include "TCanvas.h"
#include <map>
using namespace std;

int main(int argc , char *argv[])
{
  TStopwatch ti;
  ti.Start();
  const int ch_thr=18;
  gROOT->ProcessLine("#include<vector>");
	TFile f("Results.root","recreate");
  std::vector<int> ChannelTDC1toStrip={14,14,12,12,8 ,8 ,6 ,6 ,4 ,4 ,2 ,2 ,0 ,0 };
  std::vector<int> ChannelTDC2toStrip={16,16,18,18,20,20,22,22,24,24,28,28,30,30};
  std::vector<int>myfriend={1,0,3,2,5,4,7,6,9,8,11,10,13,12};
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
  TTree *dataTree = (TTree *)dataFile.Get("RawData");
  if (!dataTree) 
  {
  	std::cout << "Impossible to read TTree RAWData"<< std::endl;
    delete dataTree;
    std::exit(1);
  }
  std::vector<float>* Timestamp1 = new std::vector<float>; 
  std::vector<float>* Timestamp2 = new std::vector<float>;
  std::vector<float>* Channel1 = new std::vector<float>; 
  std::vector<float>* Channel2 = new std::vector<float>;  
  int eventNbr1;
  int eventNbr2;
  Timestamp1->clear();
  Timestamp2->clear();
  Channel1->clear();
  Channel2->clear();
  dataTree->SetBranchAddress("TimeStamp_TDC1", &Timestamp1);
  dataTree->SetBranchAddress("TimeStamp_TDC2", &Timestamp2);
  dataTree->SetBranchAddress("Channel_TDC1", &Channel1);
  dataTree->SetBranchAddress("Channel_TDC2", &Channel2);
  dataTree->SetBranchAddress("Event_NbrTDC1", &eventNbr1);
  dataTree->SetBranchAddress("Event_NbrTDC2", &eventNbr2);
  unsigned int nEntries = dataTree->GetEntries();
  TH1F* Hits=new TH1F("Hits","Hits",32,0,32);
  TH1F* gTimestamp1=new TH1F("Timestamps_trigger","Timestamps Triggers",30000,0,15);
  TH1F* gTimestamp2=new TH1F("Timestamps_trigger","Timestamps Triggers",30000,0,15);
  TH1F* diff=new TH1F("diff","diff",60000,-15,15);
  TH1F* diff1=new TH1F("diff1","diff1",60000,-15,15);
  TH1F* diff2=new TH1F("diff2","diff2",60000,-15,15);
  TH1F* gTimestampnearTrigger1=new TH1F("Timestamps_near_trigger1","Timestamps_near_trigger1",30000,0,15);
  TH1F* gTimestampnearTrigger2=new TH1F("Timestamps_near_trigger2","Timestamps_near_trigger2",30000,0,15);
  TH2F* gvision=new TH2F("selected1","selected1",32,0,32,30000,0,15);
  TH2F* real=new TH2F("real","real",32,0,32,30000,0,20);
  TH2F* gchannel=new TH2F("channel_near_trigger1","channel_near_trigger1",32,0,32,1,0,20);
  for (unsigned int i = 0; i < nEntries; i++) 
  {
  	std::vector<float> trigger1;
  	std::vector<float> trigger2;
  	std::map<float,float>Closetotrigger1;
  	std::map<float,float>Closetotrigger2;
  	trigger1.clear();
  	trigger2.clear();
  	Closetotrigger1.clear();
  	Closetotrigger2.clear();
  	dataTree->GetEntry(i);
    for (unsigned int j = 0; j !=Channel1->size(); ++j) 
    {
    	if(Channel1->at(j)<ChannelTDC1toStrip.size())
    	{
    	 //std::cout<<Channel1->at(j)<<"  "<<ChannelTDC1toStrip[Channel1->at(j)]<<std::endl;
    	 Hits->Fill(ChannelTDC1toStrip[Channel1->at(j)],0.5);
    	}
    	if(Channel1->at(j)==ch_thr)
    	{
    		trigger1.push_back(Timestamp1->at(j));
    		gTimestamp1->Fill(Timestamp1->at(j)*1e-12);
      }
    }
    for (unsigned int j = 0; j !=Channel2->size(); ++j) 
    {
    	if(Channel2->at(j)<ChannelTDC2toStrip.size())
    	{
    	 //std::cout<<Channel2->at(j)<<"  "<<ChannelTDC2toStrip[Channel2->at(j)]<<std::endl;
    	 Hits->Fill(ChannelTDC2toStrip[Channel2->at(j)],0.5);
    	}
    	if(Channel2->at(j)==ch_thr)
    	{
    		trigger2.push_back(Timestamp2->at(j));
    		gTimestamp2->Fill(Timestamp2->at(j)*1e-12);
      }
    }
    //std::cout<<trigger1.size()<<"  "<<trigger2.size()<<std::endl;
    for(unsigned int i=0;i!=trigger1.size();++i)
    {
    	for(unsigned int j=0;j!=trigger2.size();++j)
    	{
    		//std::cout<<trigger1[i]-trigger2[j]<<std::endl;
    		diff->Fill((trigger1[i]-trigger2[j])*1e-12,1);
    	}
    }
    for (unsigned int j = 0; j !=Channel1->size(); ++j) 
    {
      for(unsigned int i=0;i!=trigger1.size();++i)
      {
        if(fabs(Timestamp1->at(j)-trigger1[i])*1e-12<=100.0&&Channel1->at(j)<ChannelTDC1toStrip.size()&&Channel1->at(j)!=ch_thr)
        {
          Closetotrigger1[Timestamp1->at(j)]=Channel1->at(j);
          //std::cout<<Channel1->at(j)<<"  "<<Timestamp1->at(j)<<"  "<<trigger1[i]<<std::endl;
        }
        if(Channel1->at(j)<ChannelTDC1toStrip.size()&&Channel1->at(j)!=ch_thr)diff1->Fill((Timestamp1->at(j)-trigger1[i])*1e-9);
      }
    }
    std::map<float,float>realstrip;
    for(std::map<float,float>::iterator it=Closetotrigger1.begin();it!=Closetotrigger1.end();++it)
    {
      std::set<float>difff;
      std::map<float,float>::iterator pp=it;
      pp++;
      for(std::map<float,float>::iterator itt=pp;itt!=Closetotrigger1.end();++itt)
      {
        
        if(itt->second==myfriend[int(it->second)])difff.insert((itt->first-it->first));
      }
      std::cout<<"Values for real strip "<< ChannelTDC2toStrip[it->second]<<" : "<<std::endl;
      for(std::set<float>::iterator y=difff.begin();y!=difff.end();++y)
      {
        if(y==difff.begin())realstrip[ChannelTDC2toStrip[it->second]]=*y;
        std::cout<<*y<<std::endl;
      }
      std::set<float>::iterator ppp=difff.begin();
      
    }
    for(std::map<float,float>::iterator m=realstrip.begin();m!=realstrip.end();++m)
    {
    
      real->Fill(m->first,(m->second*1e-12*299792458000)/3+30);
      std::cout<<m->first<<"  "<<m->second*1.0e-12<<"  "<<(m->second*1.0e-12*299792458000)/3+30<<std::endl;
    }
    for (unsigned int j = 0; j !=Channel2->size(); ++j) 
    {
      for(unsigned int i=0;i!=trigger2.size();++i)
      {
        if(fabs(Timestamp2->at(j)-trigger2[i])*1e-12<=100.0&&Channel2->at(j)<ChannelTDC2toStrip.size()&&Channel2->at(j)!=ch_thr)Closetotrigger2[Timestamp2->at(j)]=Channel2->at(j);
        if(Channel2->at(j)<ChannelTDC2toStrip.size()&&Channel2->at(j)!=ch_thr)diff2->Fill((Timestamp2->at(j)-trigger2[i])*1e-12);
      }
    }
    for(std::map<float,float>::iterator it=Closetotrigger1.begin();it!=Closetotrigger1.end();++it)
    {
      if(Channel1->at(it->second)<ChannelTDC1toStrip.size())
      {
        gTimestampnearTrigger1->Fill(it->first*1e-12);
        gchannel->Fill(ChannelTDC1toStrip[Channel1->at(it->second)],1,0.5);
        gvision->Fill(ChannelTDC1toStrip[Channel1->at(it->second)],it->first*1e-12,0.5);
      }
    }
    for(std::map<float,float>::iterator it=Closetotrigger2.begin();it!=Closetotrigger2.end();++it)
    {
      if(Channel2->at(it->second)<ChannelTDC2toStrip.size())
      {
        gTimestampnearTrigger2->Fill(it->first*1e-12);
        gchannel->Fill(ChannelTDC2toStrip[Channel2->at(it->second)],1,0.5);
        //std::cout<<ChannelTDC2toStrip[Channel2->at(it->second)]<<"  "<<Channel2->at(it->second)<<"  "<<ChannelTDC2toStrip.size()<<std::endl;
        gvision->Fill(ChannelTDC2toStrip[Channel2->at(it->second)],it->first*1e-12,0.5);
      }
      //std::cout<<ChannelTDC2toStrip[Channel2->at(it->second)]<<"  "<<Channel2->at(it->second)<<std::endl;
    }
  }
	//dataFile.Close();
	f.cd();
	Hits->Write();
	real->Write();
	gTimestamp1->Write();
	gTimestamp2->Write();
	gTimestampnearTrigger2->Write();
	gchannel->Write();
	gTimestampnearTrigger1->Write();
	diff->Write();
	diff1->Write();
	diff2->Write();
	gvision->Write();
	f.Close();
	//delete Hits1;
	ti.Stop();
	//delete dataTree;
  std::cout << " Time :" << ti.RealTime() << "  " << ti.CpuTime() << std::endl;
	return 0;
}
