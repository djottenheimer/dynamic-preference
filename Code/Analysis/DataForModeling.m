%% get relevant info for each trial
clear all;
load('RAWCh');
%load('RAWCh_extra'); %for extra sessions
%%
clear NeuronInfo;

%windows for events
Window{1}=[-2 0]; %pre-cue
Window{2}=[0 0.75]; %post-cue
Window{3}=[-3 -2]; %pre-LP, relative to delivery
Window{4}=[-2 -1.25]; %post-LP, relative to delivery
Window{5}=[-1.25 0]; %pre-RD
Window{6}=[0.75 1.95]; %post-RD

Evs={'pre-cue';'post-cue';'pre-act';'post-act';
          'pre-RD';'post-RD'};


%get all the data necessary for model fitting for each neuron
NN=0;
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
 
    %licks per trial
    LickTimes=RAW(i).Erast{Licks,1};
    for trial=1:length(Trials)
        LicksPerTrial(trial,1)=sum(LickTimes>Trials(trial,1) & LickTimes<(Trials(trial,1)+13));
    end
    
    
    for j=1:length(RAW(i).Nrast) %Number of neurons within sessions
        NN=NN+1;
        
        
        NeuronInfo(NN).Rat=RAW(i).Rat;
        NeuronInfo(NN).Ses=RAW(i).Ses;
        NeuronInfo(NN).Typ=RAW(i).Ses>10;
        NeuronInfo(NN).Outcome=TrialTbl.Outcome;
        NeuronInfo(NN).Choice=TrialTbl.Choice;
        NeuronInfo(NN).Licks = LicksPerTrial;
        NeuronInfo(NN).InPort = InPort;
        NeuronInfo(NN).Window = Window;
        NeuronInfo(NN).Evs = Evs;
        
        
        
        %get trial by trial firing rate for all trials
        for event=1:length(Window)
            rewspk=[];
            [EVcell,EVn]=MakePSR04(RAW(i).Nrast(j),RAW(i).Erast{Ev{event}},Window{event},{2});% makes trial by trial rasters for event
            for y= 1:EVn
                rewspk(y,1)=sum(EVcell{1,y}>Window{event}(1));
            end
            NeuronInfo(NN).Spikes(:,event)=rewspk(1:length(Trials));
            NeuronInfo(NN).Hz(:,event)=rewspk(1:length(Trials))/(Window{event}(1,2)-Window{event}(1,1)); %convert to firing rate
            
        end
        
        
        
        
    end
    
    fprintf('Session ID # %d\n',i);
    
end