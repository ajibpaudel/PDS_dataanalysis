//===================Function filter===============================================================
void TV1D_denoise(vector<float>& input, vector<float>& output, const int width, const float lambda){
  if (width>0) {                /*to avoid invalid memory access to input[0]*/
    int k=0, k0=0;            /*k: current sample location, k0: beginning of current segment*/
    float umin=lambda, umax=-lambda;    /*u is the dual variable*/
    float vmin=input[0]-lambda, vmax=input[0]+lambda;    /*bounds for the segment's value*/
    int kplus=0, kminus=0;     /*last positions where umax=-lambda, umin=lambda, respectively*/
    const float twolambda=2.0*lambda;    /*auxiliary variable*/
    const float minlambda=-lambda;        /*auxiliary variable*/
    for (;;) {                /*simple loop, the exit test is inside*/
      while (k==width-1) {    /*we use the right boundary condition*/
	if (umin<0.0) {            /*vmin is too high -> negative jump necessary*/
	  do output[k0++]=vmin; while (k0<=kminus);
	  umax=(vmin=input[kminus=k=k0])+(umin=lambda)-vmax;
	} else if (umax>0.0) {    /*vmax is too low -> positive jump necessary*/
	  do output[k0++]=vmax; while (k0<=kplus);
	  umin=(vmax=input[kplus=k=k0])+(umax=minlambda)-vmin;
	} else {
	  vmin+=umin/(k-k0+1);
	  do output[k0++]=vmin; while(k0<=k);
	  return;
	}
      }
      if ((umin+=input[k+1]-vmin)<minlambda) {        /*negative jump necessary*/
	do output[k0++]=vmin; while (k0<=kminus);
	vmax=(vmin=input[kplus=kminus=k=k0])+twolambda;
	umin=lambda; umax=minlambda;
      } else if ((umax+=input[k+1]-vmax)>lambda) {    /*positive jump necessary*/
	do output[k0++]=vmax; while (k0<=kplus);
	vmin=(vmax=input[kplus=kminus=k=k0])-twolambda;
	umin=lambda; umax=minlambda;
      } else {     /*no jump necessary, we continue*/
	k++;
	if (umin>=lambda) {        /*update of vmin*/
	  vmin+=(umin-lambda)/((kminus=k)-k0+1);
	  umin=lambda;
	}
	if (umax<=minlambda) {    /*update of vmax*/
	  vmax+=(umax+lambda)/((kplus=k)-k0+1);
	  umax=minlambda;
	}
      }
    }
  }
}
//==================================================================================================




void PDSsignalanalysis(){
  TFile *f=new TFile("dvdm-2kfeedback-100ohm-adapter-45v.root","RECREATE");
  
  TH1D *Total_charge=new TH1D("Total_charge","Integrated charge;Charge [V*ns];entries",1120,-0.4*4,2.4*4);
  TH2D *time_vs_signal=new TH2D("time_vs_signal",";time [ticks];signal[V]",4000,0,4000,2000,-0.1,0.1);
  TH1D *amplitude =new TH1D("amplitude","",1000,-0.05,0.05);
  TProfile *avg_SPE=new TProfile("avg_SPE","Mean SPE;time [ns];signal [V]",4000,0,16000);


  std::vector<float> ondama(4000);
  std::vector<float> ondaden(4000);
  std::vector<float> input(4000);
  std::vector<float> output(4000);
   Float_t mobileAVG,avg,c0,c1,c2;
   double mu,sigma,sigma2,sv1,sv2,base,baseCC,mediaPre,sum,sum2,base2;
   int index;
   int inum=0;
   double ampl=0;//amplitude

   
  
  TGraph *gr[4000];
  TGraph *gr_denoised[4000];
  double times[4000]={0.0};
  double waveform[4000]={0.0};
  int entry=0;
  std::vector<float> timev(4000);
  for(int i=0;i<4000;i++){
    times[i]=i;
    timev[i]=i;
  }
  
  //  ifstream file("test_save_led6p8.txt");
  // ifstream file("july8_led6p8_ch3_36V.txt");
  // ifstream file("july8_led6p8_ch2_36V_1.txt");
  //ifstream file("test_save_led10p0-35v.txt");
  ifstream file("dvdm-2kfeedback-100ohm-adapter-45v.txt");
  int line=0;
  vector<double> waveformy, time;
  double y, x, sizev;
  int entries=0;
  while (!file.eof()){
    file>>y>>x;
    line++;
    waveformy.push_back(y);
    if(line<1000) continue;
    if(x>1.0){
      sizev=waveformy.size();
      for(int i=0;i<1000;i++){
	waveform[999-i]=waveformy[sizev-1-i];
      }
      for(int i=1000;i<4000;i++){
	file>>y>>x;
	waveform[i]=y;
      }
      waveformy.clear();
      entry++;


      /////////////////Charge analysis///////////////////

      //////////////////////////////////////////////////
      //========= Denoising ============================================================
      mobileAVG=4.0;  avg=0.;  c0=0.;  c1=0.;  c2=0.;
      for(Int_t n=0; n<4000 ; n++){if(n>=mobileAVG && n<4000-mobileAVG){
	  for(Int_t i1=n-mobileAVG; i1<=n+mobileAVG; i1++){avg=avg+waveform[i1]; c0=c0+1;}
	  avg=avg/c0; ondama[n]=avg; avg=0; c0=0;}
	else{if(n<mobileAVG){
	    for(Int_t i1=0; i1<=n+mobileAVG; i1++){avg=avg+waveform[i1]; c1=c1+1;}
	    avg=avg/c1; ondama[n]=avg; avg=0; c1=0;}
	  else if(n>=4000-mobileAVG){
	    for(Int_t i1=n-mobileAVG; i1<4000; i1++){avg=avg+waveform[i1]; c2=c2+1;}
	    avg=avg/c2; ondama[n]=avg; avg=0; c2=0;}}}
            
      for(Int_t i1=0; i1<4000; i1++){input[i1]=ondama[i1]; output[i1]=input[i1];}
            
      TV1D_denoise(input,output,4000,0.005);
      for(Int_t i1=0; i1<4000; i1++){
	ondaden[i1]=output[i1];
      }
            
      //===========================BASELINE Histo========================================
            
      //===========================BASELINE Histo========================================
      base=0;
      base2=0;
      
      TH1F *basehelp= new TH1F("basehelp","basehelp",500,-0.01,0.01);//channel 4
      for(int ib=0; ib<800; ib++){basehelp->Fill(ondaden[ib]);}
      /* basebinmax = basehelp->GetMaximumBin();
	 base = basehelp->GetXaxis()->GetBinCenter(basebinmax);*/
      base=basehelp->GetMean();//try using mean for baseline
      
      basehelp->Delete();
      for(int i=0;i<4000;i++){
	waveform[i]=waveform[i]-base;
	time_vs_signal->Fill(i,waveform[i]);
      }

       for(int i=0; i<4000; i++){
	ondaden[i]=ondaden[i]-base;
      }
       
      if(entry<1000){
	gr[entry]=new TGraph(4000,times,waveform);
	gr[entry]->Write(Form("gr_raw_%d",entry));
	gr_denoised[entry]=new TGraph(timev.size(),&timev[0],&ondaden[0]);
	gr_denoised[entry]->Write(Form("gr_denoised_%d",entry));
      }
      
     
      
      //=================== Integral ====================================================
      sum=0.; sum2=-4000.; sigma2=-10;ampl=-100;index=0;
      for(int i=0;i<4000;i++){
	if(i>1000 && i<1600){//noise
	  sum+=ondaden[i];
	  if(waveform[i]>ampl) ampl=waveform[i];
	}

      }
     
      //=================================================================================
      //=================================================================================




      ////////////////////////////////////////////////
      ///////////////////////////////////////////////

      Total_charge->Fill(sum*4);
      if(sum*4>0.618-0.094 && sum*4<0.618+0.094){
	for(int i=0;i<4000;i++){
	  avg_SPE->Fill(i*4,waveform[i]);
	}
      }
      
      amplitude->Fill(ampl);
    }
  }

  double timesspe[4000];
  double sigspe[4000];
  for(int i=0;i<4000;i++){
    timesspe[i]=4*i;
    sigspe[i]=avg_SPE->GetBinContent(i+1);
  }
  TGraph *gr_spe=new TGraph(4000,timesspe,sigspe);
  gr_spe->SetTitle("Average SPE signal;time [ns];signal [V]");
  gr_spe->Write();
  Total_charge->Write();
  amplitude->Write();
  time_vs_signal->Write();
  Total_charge->Draw();
  avg_SPE->Write();
  f->Close();
  
 
}
