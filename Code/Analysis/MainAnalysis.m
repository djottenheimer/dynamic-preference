%% get relevant info for each trial
load('RAWCh');
%load('RAWCh_extra');

global Dura Tm BSIZE
BSIZE=0.01; %Do not change
Dura=[-5 15]; Tm=Dura(1):BSIZE:Dura(2);

R.Tm=Tm;

%windows for events
Window{1}=[-2 0]; %pre-cue
Window{2}=[0 0.75]; %post-cue
Window{3}=[-3 -2]; %pre-LP, relative to delivery
Window{4}=[-2 -1.25]; %post-LP, relative to delivery
Window{5}=[-1.25 0]; %pre-RD
Window{6}=[0.75 1.95]; %post-RD
Baseline=[-6 -1]; %for PSTHs

R.Events={'SucF1';'SucF2';'SucF3';'SucF4';
          'WatF1';'WatF2';'WatF3';'WatF4';
          'CueS1';'CueS2';'CueS3';'CueS4';
          'CueW1';'CueW2';'CueW3';'CueW4';
          'CueC1';'CueC2';'CueC3';'CueC4';};

%Smoothing PSTHs
%smoothing filter for whole PSTH
PSTHsmoothbins=25; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',8); %used to be 6.6
PSTHfilterweights=pdf(halfnormal,0:PSTHsmoothbins);

%smoothing filter for individual trials
trialsmoothbins=10; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',3); %std=3.98
trialfilterweights=pdf(halfnormal,0:trialsmoothbins);


%smoothing filter for licking PSTH
licksmoothbins=50; %number of previous bins used to smooth
halfnormal=makedist('HalfNormal','mu',0,'sigma',25); %std=3.98
lickfilterweights=pdf(halfnormal,0:licksmoothbins);

%preference smoothing filter
sigma=5;
smoothbins=10; %number of trials before and after current using in smoothing
normal=makedist('Normal','mu',0,'sigma',sigma);
filterweights=pdf(normal,-smoothbins:smoothbins);


NN=0;
R.GLM=[];
R.Hz=[];
%%
tic
for i=1:length(RAW)
    %events
    SucroseF=strmatch('SucroseF',RAW(i).Einfo(:,2),'exact');
    SucroseC=strmatch('SucroseC',RAW(i).Einfo(:,2),'exact');
    WaterF=strmatch('WaterF',RAW(i).Einfo(:,2),'exact');
    WaterC=strmatch('WaterC',RAW(i).Einfo(:,2),'exact');
    Cues=strmatch('CueAll',RAW(i).Einfo(:,2),'exact');
    CueS=strmatch('CueS',RAW(i).Einfo(:,2),'exact');
    CueC=strmatch('CueC',RAW(i).Einfo(:,2),'exact');
    CueW=strmatch('CueW',RAW(i).Einfo(:,2),'exact');
    LeverR=strmatch('LeverR',RAW(i).Einfo(:,2),'exact');
    LeverS=strmatch('LeverS',RAW(i).Einfo(:,2),'exact');
    LeverW=strmatch('LeverW',RAW(i).Einfo(:,2),'exact');
    Delivery=strmatch('Delivery',RAW(i).Einfo(:,2),'exact');
    PortEntry=strmatch('PortEntry',RAW(i).Einfo(:,2),'exact');
    PortExit=strmatch('PortExit',RAW(i).Einfo(:,2),'exact');
    Licks=strmatch('Licks',RAW(i).Einfo(:,2),'exact');
    
    Ev{1}=Cues;
    Ev{2}=Cues;
    Ev{3}=Delivery;
    Ev{4}=Delivery;
    Ev{5}=Delivery;
    Ev{6}=Delivery;

    TrialTbl=[];
    TrialTbl=table();

    Trials1=cat(2,RAW(i).Erast{SucroseF,1},ones(length(RAW(i).Erast{SucroseF,1}),1));
    Trials2=cat(2,RAW(i).Erast{SucroseC,1},2*ones(length(RAW(i).Erast{SucroseC,1}),1));
    Trials3=cat(2,RAW(i).Erast{WaterF,1},3*ones(length(RAW(i).Erast{WaterF,1}),1));
    Trials4=cat(2,RAW(i).Erast{WaterC,1},4*ones(length(RAW(i).Erast{WaterC,1}),1));

    Trials=cat(1,Trials1,Trials2,Trials3,Trials4);
    [Trials(:,1),ind]=sort(Trials(:,1));
    Trials(:,2)=Trials(ind,2);


    TrialTbl.Outcome = Trials(:,2)==1 | Trials(:,2)==2; %1 if sucrose
    TrialTbl.Choice = Trials(:,2)==2 | Trials(:,2)==4; %1 if choice trial

    TrialNo=1:length(Trials);
    
    %find trials where rat was in port for delivery
    InPort=[];
    for trial=1:length(RAW(i).Erast{Delivery})
        InPort(trial,1)=sum(RAW(i).Erast{Delivery}(trial)>RAW(i).Erast{PortEntry} & RAW(i).Erast{Delivery}(trial)<RAW(i).Erast{PortExit});
    end
    TrialTbl.InPort=InPort;
    %preference throughout session


    Choices = TrialTbl.Outcome(TrialTbl.Choice==1 & TrialTbl.InPort);
    ChoiceTimes = Trials(TrialTbl.Choice==1 & TrialTbl.InPort,1);
    Preference=[];
    for l=1:length(Choices)
       Preference(1,l)=sum(Choices(l-min([l-1 smoothbins]):l+min([length(Choices)-l smoothbins]),1)'.*fliplr(filterweights(1+max([0 smoothbins+l-length(Choices)]):smoothbins+min([l smoothbins+1]))))/sum(filterweights(1+max([0 smoothbins+l-length(Choices)]):smoothbins+min([l smoothbins+1])));                  
    end


    TrialPref=[];
    for trial=1:length(Trials)
        [~,closest_choice]=min(abs(ChoiceTimes-Trials(trial,1)));
        TrialPref(trial,1)=Preference(closest_choice);
    end
    
    %licks per trial
    LickTimes=RAW(i).Erast{Licks,1};
    for trial=1:length(Trials)
        LicksPerTrial(trial,1)=sum(LickTimes>Trials(trial,1) & LickTimes<(Trials(trial,1)+13));
    end
    
    TrialTbl.Preference = TrialPref;
    TrialTbl.TrialNo = TrialNo';
    TrialTbl.RDTime = Trials(:,1);
    R.TrialInfo{i,1}=TrialTbl;
    
    Break1=length(Trials)/4;
    Break2=length(Trials)/2;
    Break3=3*length(Trials)/4;
    

    %for reward analysis, make sure rat is in port during delivery
    Q1=TrialTbl.TrialNo<Break1 & TrialTbl.InPort;
    Q2=TrialTbl.TrialNo>=Break1 & TrialTbl.TrialNo<Break2 & TrialTbl.InPort;
    Q3=TrialTbl.TrialNo>=Break2 & TrialTbl.TrialNo<Break3 & TrialTbl.InPort;
    Q4=TrialTbl.TrialNo>=Break3 & TrialTbl.InPort;
    
    %for cue analysis, doesn't matter
    Q1cue=TrialTbl.TrialNo<Break1;
    Q2cue=TrialTbl.TrialNo>=Break1 & TrialTbl.TrialNo<Break2;
    Q3cue=TrialTbl.TrialNo>=Break2 & TrialTbl.TrialNo<Break3;
    Q4cue=TrialTbl.TrialNo>=Break3;

    FS=TrialTbl.Outcome==1 & TrialTbl.Choice==0;
    FW=TrialTbl.Outcome==0 & TrialTbl.Choice==0;
    CS=TrialTbl.Outcome==1 & TrialTbl.Choice==1;
    CW=TrialTbl.Outcome==0 & TrialTbl.Choice==1;
    Ch=TrialTbl.Choice==1;
    
    EvTimestamps{1,1}=RAW(i).Erast{SucroseF,1}(Q1(FS));
    EvTimestamps{2,1}=RAW(i).Erast{SucroseF,1}(Q2(FS));
    EvTimestamps{3,1}=RAW(i).Erast{SucroseF,1}(Q3(FS));
    EvTimestamps{4,1}=RAW(i).Erast{SucroseF,1}(Q4(FS));
    EvTimestamps{5,1}=RAW(i).Erast{WaterF,1}(Q1(FW));
    EvTimestamps{6,1}=RAW(i).Erast{WaterF,1}(Q2(FW));
    EvTimestamps{7,1}=RAW(i).Erast{WaterF,1}(Q3(FW));
    EvTimestamps{8,1}=RAW(i).Erast{WaterF,1}(Q4(FW));
    EvTimestamps{9,1}=RAW(i).Erast{CueS,1}(Q1cue(FS));
    EvTimestamps{10,1}=RAW(i).Erast{CueS,1}(Q2cue(FS));
    EvTimestamps{11,1}=RAW(i).Erast{CueS,1}(Q3cue(FS));
    EvTimestamps{12,1}=RAW(i).Erast{CueS,1}(Q4cue(FS));
    EvTimestamps{13,1}=RAW(i).Erast{CueW,1}(Q1cue(FW));
    EvTimestamps{14,1}=RAW(i).Erast{CueW,1}(Q2cue(FW));
    EvTimestamps{15,1}=RAW(i).Erast{CueW,1}(Q3cue(FW));
    EvTimestamps{16,1}=RAW(i).Erast{CueW,1}(Q4cue(FW));
    EvTimestamps{17,1}=RAW(i).Erast{CueC,1}(Q1cue(Ch));
    EvTimestamps{18,1}=RAW(i).Erast{CueC,1}(Q2cue(Ch));
    EvTimestamps{19,1}=RAW(i).Erast{CueC,1}(Q3cue(Ch));
    EvTimestamps{20,1}=RAW(i).Erast{CueC,1}(Q4cue(Ch));
    

%%



    for j=1:length(RAW(i).Nrast) %Number of neurons within sessions
        NN=NN+1;
        
        
        R.Rat(NN,1)=RAW(i).Rat;
        R.Ses(NN,1)=RAW(i).Ses;
        R.Typ(NN,1)=RAW(i).Ses>10;
        

        %get trial by trial firing rate for all trials
        for event=1:length(Window)
            rewspk=[];
            [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{Ev{event}},Window{event},{2});% makes trial by trial rasters for event
            for y= 1:EVn
                rewspk(y,1)=sum(EVcell{1,y}>Window{event}(1));
            end       
            R.Hz{i,event}(:,j)=rewspk(1:length(Trials))/(Window{event}(1,2)-Window{event}(1,1)); %convert to firing rate

            for trialtype=1:2
                if trialtype==1 Sel=TrialTbl.Choice==0; end
                if trialtype==2 Sel=TrialTbl.Choice==1; end
                tbl=[TrialTbl(Sel,:) table(R.Hz{i,event}(Sel,j),'VariableNames',{'Hz'})];
                modelspec = 'Hz ~ Outcome*TrialNo'; %this used to be preference instead of trial number
                mdl = fitglm(tbl,modelspec,'Distribution','poisson');

                R.GLM{trialtype,1}.Coefficients.Outcome(NN,event)=mdl.Coefficients.Estimate(2);
                R.GLM{trialtype,1}.Coefficients.Preference(NN,event)=mdl.Coefficients.Estimate(3);
                R.GLM{trialtype,1}.Coefficients.OutByPref(NN,event)=mdl.Coefficients.Estimate(4);
                R.GLM{trialtype,1}.PVal.Outcome(NN,event)=mdl.Coefficients.pValue(2);
                R.GLM{trialtype,1}.PVal.Preference(NN,event)=mdl.Coefficients.pValue(3);
                R.GLM{trialtype,1}.PVal.OutByPref(NN,event)=mdl.Coefficients.pValue(4);
                R.GLM{trialtype,1}.R2(NN,event)=mdl.Rsquared.Ordinary;
            end
        end
        
        %get mean baseline firing for all reward trials
        [Bcell1,B1n]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{Cues},Baseline,{2});% makes trial by trial rasters for baseline
        for y= 1:B1n
            basespk(1,y)=sum(Bcell1{1,y}>Baseline(1));
        end
        
        FiringRate=basespk/(Baseline(1,2)-Baseline(1,1));
        Bmean=nanmean(FiringRate);
        Bstd=nanstd(FiringRate);
        
        R.Bmean(NN,1)=Bmean;
        R.Bstd(NN,1)=Bstd;
        for event=1:length(R.Events)
            if length(EvTimestamps{event})>=3
                
                %make PSTHs
                [PSR1,N1]=MakePSR04(RAW(i).Nrast(j),EvTimestamps{event},Dura,{2});% makes collpased rasters. PSR1 is a cell(neurons)
                smoothedtrials=[];
                binned=[];
                for trial=1:length(PSR1)
                    binned(trial,:)=histcounts(PSR1{trial},[Tm Tm(end)+(Tm(end)-Tm(end-1))]);
                    for l=1:length(Tm)
                        smoothedtrials(trial,l)=sum(binned(trial,l-min([l-1 trialsmoothbins]):l).*fliplr(trialfilterweights(1:min([l trialsmoothbins+1]))))/sum(trialfilterweights(1:min([l trialsmoothbins+1])));
                    end
                end
                PTH1=mean(smoothedtrials,1)/BSIZE;
                PTH1nosmooth=mean(binned,1)/BSIZE;
                
                PTH1smooth=[];
                for l=1:length(Tm)
                    PTH1smooth(1,l)=sum(PTH1(1,l-min([l-1 PSTHsmoothbins]):l).*fliplr(PSTHfilterweights(1:min([l PSTHsmoothbins+1]))))/sum(PSTHfilterweights(1:min([l PSTHsmoothbins+1])));
                end
                R.Ev(event).PSTHraw(NN,:)=PTH1smooth;
                R.Ev(event).PSTHnosmooth(NN,1:length(Tm))=PTH1nosmooth;

                %normalize already smoothed activity
                for l=1:length(PTH1smooth)
                    R.Ev(event).PSTHz(NN,l)=(PTH1smooth(l)-Bmean)/Bstd;
                    R.Ev(event).PSTHnosmoothz(NN,l)=(PTH1nosmooth(l)-Bmean)/Bstd;
                end
                
                %behavioral measures
                if j==1 %only do this once per session
                    
                    %licking
                    Licks=strmatch('Licks',RAW(i).Einfo(:,2),'exact');
                     %make PSTHs
                    [PSR1,N1]=MakePSR04(RAW(i).Erast(Licks),EvTimestamps{event},Dura,{2});% makes collpased rasters. PSR1 is a cell(neurons)
                    smoothedtrials=[];
                    binned=[];
                    bigbinned=[];
                    bigbins=1:5:length(Tm); %make the bins 5x as big
                    bigbins=[bigbins length(Tm)+1];
                    for trial=1:length(PSR1)
                        binned(trial,:)=histcounts(PSR1{trial},[Tm Tm(end)+(Tm(end)-Tm(end-1))]);
                        for l=1:length(Tm)
                            bigbinned(trial,l)=mean(binned(trial,max(bigbins(bigbins<=l)):min(bigbins(bigbins>l))-1));
                            smoothedtrials(trial,l)=sum(bigbinned(trial,l-min([l-1 trialsmoothbins]):l).*fliplr(trialfilterweights(1:min([l trialsmoothbins+1]))))/sum(trialfilterweights(1:min([l trialsmoothbins+1])));
                        end
                    end
                    PTH1=mean(smoothedtrials,1)/BSIZE;
                    PTH1nosmooth=mean(binned,1)/BSIZE;

                    PTH1smooth=[];
                    for l=1:length(Tm)
                        PTH1smooth(1,l)=sum(PTH1(1,l-min([l-1 licksmoothbins]):l).*fliplr(lickfilterweights(1:min([l licksmoothbins+1]))))/sum(lickfilterweights(1:min([l licksmoothbins+1])));
                    end
                    R.Lick(event).PSTHraw(i,:)=PTH1smooth;
                    R.Lick(event).PSTHnosmooth(i,1:length(Tm))=PTH1nosmooth;

                    
                    
                    %port entry
                    %make PSTHs
                    inport=[];
                    for trial=1:length(EvTimestamps{event})
                        for time=1:length(Tm)
                            inport(trial,time)=sum(EvTimestamps{event}(trial,1)+Tm(time)>RAW(i).Erast{PortEntry} & EvTimestamps{event}(trial,1)+Tm(time)<RAW(i).Erast{PortExit});
                        end
                    end
                    PTH1=mean(inport,1);
                    R.Port(event).PSTHnosmooth(i,:)=PTH1;
                    
                    PTH1smooth=[];
                    for l=1:length(Tm)
                        PTH1smooth(1,l)=sum(PTH1(1,l-min([l-1 PSTHsmoothbins]):l).*fliplr(PSTHfilterweights(1:min([l PSTHsmoothbins+1]))))/sum(PSTHfilterweights(1:min([l PSTHsmoothbins+1])));
                    end
                    R.Port(event).PSTHraw(i,:)=PTH1smooth;
                end
                
                
            else
                R.Ev(event).PSTHraw(NN,1:length(Tm))=NaN;
                R.Ev(event).PSTHz(NN,1:length(Tm))=NaN;
                R.Ev(event).PSTHz(NN,1:length(Tm))=NaN;
                R.Ev(event).PSTHnosmoothz(NN,1:length(Tm))=NaN;
                R.Lick(event).PSTHraw(i,:)=NaN;
                R.Lick(event).PSTHnosmooth(i,1:length(Tm))=NaN;
                R.Lick(event).PSTHz(i,1:length(Tm))=NaN;
                R.Lick(event).PSTHnosmoothz(i,1:length(Tm))=NaN;
                R.Port(event).PSTHraw(i,:)=NaN;
            end
        end


        fprintf('Neuron ID # %d\n',NN);
    end

end

%save('R_choice.mat','R')

toc