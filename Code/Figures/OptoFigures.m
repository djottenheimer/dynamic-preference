%get cue and reward delivery info from the medpc file
clear all;

gfp=[4 6 10 11 13 14 15 20 21];
chr=[1 2 7 8 9 16 17 18 22 23 24];

%file location
address=['/Users/david/Documents/GitHub/dynamic-preference/Data/OptoFiles'];

%files=dir([address,'\\*!*']);
files=dir([address,'//*!*']); %mac

%get day of training from date
for i=1:length(files)
    date(i)=str2num(files(i).name([8 10:11]));
end
dates=unique(date);

for k=1:length(files)
    filename=fullfile(files(k).folder,files(k).name);
    file=fopen(filename);
    L=textscan(file,'%s','Delimiter',':');
    fclose('all');
    x=str2num(files(k).name(30:31)); %x is rat #
    y=str2num(files(k).name([8 10:11])); %y is date
    y=find(dates==y);
    
    %find start of A, start of C, and end of C
    Hstrt=find(strcmp('H',L{1,1})); %Sucrose cue onset
    Istrt=find(strcmp('I',L{1,1})); %Maltodextrin cue onset
    Jstrt=find(strcmp('J',L{1,1})); %Choice
    Kstrt=find(strcmp('K',L{1,1})); %Sucrose Delivery
    Lstrt=find(strcmp('L',L{1,1})); %Maltodextrin Delivery
    Mstrt=find(strcmp('M',L{1,1})); %LLP (maltodextrin)
    Nstrt=find(strcmp('N',L{1,1})); %RLP (sucrose)
    Ostrt=find(strcmp('O',L{1,1})); %Port entry times
    Pstrt=find(strcmp('P',L{1,1})); %Port durations
    Qstrt=find(strcmp('Q',L{1,1})); %Laser stimulations
    Tstrt=find(strcmp('X',L{1,1})); %List (marks end of laser)
    
    %Sucrose trial times
    suct=(L{1,1}(Hstrt+2:2:Istrt-1));
    sucrosettime=[];
    for i=1:length(suct)
        medtext=textscan(char(suct(i,1)),'%f','Delimiter',' ');
        stime=medtext{1,1}(:,1);
        sucrosettime=(cat(1,sucrosettime,stime(isnan(stime)==0)));
    end
    sucrosettime(sucrosettime==0)=[]; %delete 0s
    
    %deals with session where commutator fell off
    if x==8 & y==13
        sucrosettime(sucrosettime>2580)=[];
    end
    
    sucrosetimes{x,y}=sucrosettime;
    
    %Maltodextrin trial times
    maltot=(L{1,1}(Istrt+2:2:Jstrt-1));
    maltotime=[];
    for i=1:length(maltot)
        medtext=textscan(char(maltot(i,1)),'%f','Delimiter',' ');
        mtime=medtext{1,1}(:,1);
        maltotime=(cat(1,maltotime,mtime(isnan(mtime)==0)));
    end
    maltotime(maltotime==0)=[]; %delete 0s
    
    %deals with session where commutator fell off
    if x==8 & y==13
        maltotime(maltotime>2580)=[];
    end
    
    maltotimes{x,y}=maltotime;
    
    %Choice trial times
    choicet=(L{1,1}(Jstrt+2:2:Kstrt-1));
    choicetime=[];
    for i=1:length(choicet)
        medtext=textscan(char(choicet(i,1)),'%f','Delimiter',' ');
        ctime=medtext{1,1}(:,1);
        choicetime=(cat(1,choicetime,ctime(isnan(ctime)==0)));
    end
    choicetime(choicetime==0)=[]; %delete 0s
    
    if x==8 & y==13
        choicetime(choicetime>2580)=[];
    end
    
    %Sucrose delivery times
    sucrosedel=(L{1,1}(Kstrt+2:2:Lstrt-1));
    sucrosedeltime=[];
    for i=1:length(sucrosedel)
        medtext=textscan(char(sucrosedel(i,1)),'%f','Delimiter',' ');
        sdtime=medtext{1,1}(:,1);
        sucrosedeltime=(cat(1,sucrosedeltime,sdtime(isnan(sdtime)==0)));
    end
    sucrosedeltime(sucrosedeltime==0)=[]; %delete 0s
    
    %Maltodextrin delivery times
    maltdel=(L{1,1}(Lstrt+2:2:Mstrt-1));
    maldeltime=[];
    for i=1:length(maltdel)
        medtext=textscan(char(maltdel(i,1)),'%f','Delimiter',' ');
        mdtime=medtext{1,1}(:,1);
        maldeltime=(cat(1,maldeltime,mdtime(isnan(mdtime)==0)));
    end
    maldeltime(maldeltime==0)=[]; %delete 0s
    
    %Port entry times
    PEt=(L{1,1}(Ostrt+2:2:Pstrt-1));
    petimes=[];
    for i=1:length(PEt)
        medtext=textscan(char(PEt(i,1)),'%f','Delimiter',' ');
        ptimes=medtext{1,1}(:,1);
        petimes=(cat(1,petimes,ptimes(isnan(ptimes)==0)));
    end
    petimes(petimes==0)=[]; %delete 0s
    
    %Left lever presses
    LLP=(L{1,1}(Mstrt+2:2:Nstrt-1));
    LLPtime=[];
    for i=1:length(LLP)
        medtext=textscan(char(LLP(i,1)),'%f','Delimiter',' ');
        ltime=medtext{1,1}(:,1);
        LLPtime=(cat(1,LLPtime,ltime(isnan(ltime)==0)));
    end
    LLPtime(LLPtime==0)=[]; %delete 0s
    
    LLPtimes{x,y}=LLPtime;
    
    %Right lever presses
    RLP=(L{1,1}(Nstrt+2:2:Ostrt-1));
    RLPtime=[];
    for i=1:length(RLP)
        medtext=textscan(char(RLP(i,1)),'%f','Delimiter',' ');
        rtime=medtext{1,1}(:,1);
        RLPtime=(cat(1,RLPtime,rtime(isnan(rtime)==0)));
    end
    
    RLPtime(RLPtime==0)=[]; %delete 0s
    
    RLPtimes{x,y}=RLPtime;
    
    
    
    %Sucrose or maltodextrin?
    choice=[];
    for i=1:length(choicetime)
        malt=find(maldeltime(:,1) > choicetime(i,1),1);
        sucrose=find(sucrosedeltime(:,1) > choicetime(i,1),1);
        malttime=maldeltime(malt);
        sucrosetime=sucrosedeltime(sucrose);
        if malttime < sucrosetime
            choice(i,1)=0;
        elseif (isempty(sucrosetime) & isempty(malttime)==0)==1
            choice(i,1)=0;
        elseif malttime > sucrosetime
            choice(i,1)=1;
        elseif isempty(sucrosetime)==0 & isempty(malttime)
            choice(i,1)=1;
        elseif isempty(malttime) && isempty(sucrosetime)
            choice(i,1)=NaN;
        end
        
        
    end
    
    if y>=12 %12 is first stim session
        
        
        %Laser stims
        ST=(L{1,1}(Qstrt+2:2:Tstrt-1));
        STtime=[];
        for i=1:length(ST)
            medtext=textscan(char(ST(i,1)),'%f','Delimiter',' ');
            stime=medtext{1,1}(:,1);
            STtime=(cat(1,STtime,stime(isnan(stime)==0)));
        end
        
        STtime(STtime==0)=[]; %delete 0s
        %create array with all trials and info about them
        TrialTbl=[];
        TrialTbl=table();
        
        Trials1=cat(2,sucrosettime,ones(length(sucrosettime),1));
        Trials2=cat(2,maltotime,2*ones(length(maltotime),1));
        Trials3=cat(2,choicetime(isnan(choice)==0),3*ones(sum(isnan(choice)==0),1));
        
        Trials=cat(1,Trials1,Trials2,Trials3);
        [Trials(:,1),ind]=sort(Trials(:,1));
        Trials(:,2)=Trials(ind,2);
        
        
        TrialTbl.Outcome = Trials(:,2)==1; %1 if sucrose
        TrialTbl.Outcome(Trials(:,2)==3) = choice(isnan(choice)==0);
        TrialTbl.Choice = Trials(:,2)==3; %1 if choice trial
        
        %effect of stim on choice
        stimmed=zeros(length(Trials),1);
        poststimmed=zeros(length(Trials),1);
        for stim=1:length(STtime)
            [~,entry]=max(Trials(Trials(:,1)<STtime(stim),1));
            stimmed(entry)=1;
            if entry<length(Trials)
                poststimmed(entry+1)=1;
            end
        end
        
        %calculate latency
        alllptimes = sort(cat(1,RLPtime,LLPtime));
        latency=[];
        for trial = 1:length(Trials)
            if sum(alllptimes>Trials(trial,1))>0
                latency(trial,1) = min(alllptimes(alllptimes>Trials(trial,1))) - Trials(trial,1);
            end
        end
        alllatencies{x,y} = latency;
        stimtrls{x,y} = stimmed;
        poststimtrls{x,y} = poststimmed;
        
        TrialTbl.Stimmed = stimmed;
        
        
    end
    
    
    %Converting choice time into minutes
    %     choicetime=choicetime/60;
    
    %Deleting rows with NaN values
    choice_array=[];
    choice_array=cat(2,choicetime,choice);
    choice_array(any(isnan(choice_array),2),:) = [];
    
    %Numbering choice trials
    choicetrial=[];
    for i=1:length(choice_array)
        choicetrial(i,1)=i;
    end
    
    %Making choice array
    Choice{x,y}=cat(2,choicetrial,choice_array);
    
    
    %get total number of presses
    
    totalpresses(x,y)=sum(isnan(sucrosedeltime)==0)+sum(isnan(maldeltime)==0);
    presspermin(x,y)=60*totalpresses(x,y)/petimes(end);
    
    totalpes(x,y)=length(petimes);
    pepermin(x,y)=60*totalpes(x,y)/petimes(end);
    
    %latency to press sucrose
    sucroselatency=[];
    for trial = 1:length(sucrosettime)
        if sum(RLPtime>sucrosettime(trial))>0
            sucroselatency(trial,1) = min(RLPtime(RLPtime>sucrosettime(trial)))-sucrosettime(trial);
        end
    end
    
    %latency to press maltodextrin
    maltodextrinlatency=[];
    for trial = 1:length(maltotime)
        if sum(LLPtime>maltotime(trial))>0
            maltodextrinlatency(trial,1) = min(LLPtime(LLPtime>maltotime(trial)))-maltotime(trial);
        end
    end
    
    suclatencies{x,y} = sucroselatency;
    mallatencies{x,y} = maltodextrinlatency;
    
    
    
    
end

included=ones(length(Choice(:,1)),1);
included=logical(included);

for x=1:length(Choice(:,1))
    for y=1:size(Choice,2)
        if isempty(Choice{x,y}) | isnan(Choice{x,y})
            choicea(x,y)=NaN;
            Choice{x,y}=NaN;
            presspermin(x,y)=NaN;
            pepermin(x,y)=NaN;
        else
            choicea(x,y)=mean((Choice{x,y}(:,3)));
        end
    end
end
choiceavg=nanmean(choicea);


%% performance for days of interest

sessions=[12 13 14 15 16 17 18];
colors{1,1}=[0.5 0.5 0.5];
colors{2,1}=[0.2 0.4 0.8];

figure;
subplot(1,1,1)
hold on;
for i=1:length(sessions)
    scatter(rand([length(gfp) 1])/2+(i-1)*3+0.65,choicea(gfp,sessions(i)),55,colors{1},'filled');
    scatter(rand([length(chr) 1])/2+(i-1)*3+1.35,choicea(chr,sessions(i)),55,colors{2},'filled');
end

plot([0 3*(length(sessions)+1)],[0.5 0.5],':','color','k');
ylabel('Preference');
axis([0.5 length(sessions)*3+0.5 0 1.02]);
xticks(1:3:length(sessions)*3);
xticklabels({'Baseline','Test','Recovery 1','Recovery 2','Recovery 3','Recovery 4','Recovery 5'});
xtickangle(45);
yticks([0 0.5 1]);
yticklabels([-1;0;1]);
legend('YFP','ChR2');

%stats
for i=1:length(sessions)
    [~,stats_vs_control(i,1)]=kstest2(choicea(gfp,sessions(i)),choicea(chr,sessions(i)));
    [~,stats_vs_baseline(i,1)]=kstest2(choicea(chr,sessions(1)),choicea(chr,sessions(i)));
end

data=cat(1,choicea(gfp,sessions),choicea(chr,sessions));
day=[];
group=[];
ratname=[];
for i=1:length(sessions)
    day=cat(2,day,i*ones(length(gfp)+length(chr),1));
    group=cat(2,group,cat(1,zeros(length(gfp),1),ones(length(chr),1)));
    ratname=cat(2,ratname,cat(1,gfp',chr'));
end

[~,test_tbl,test_stats]=anovan(data(:),{group(:),day(:)},'varnames',{'group','day'},'model','full','display','off');
%c= multcompare(test_stats,'dimension',[1,2]);

%% preference on session of interest
rats=length(Choice(:,1));

session=13;
groups{1,1}=gfp;
groups{2,1}=chr;

%smooth the preference like ephys
%preference smoothing filter
sigma=5;
smoothbins=10; %number of trials before and after current using in smoothing
normal=makedist('Normal','mu',0,'sigma',sigma);
filterweights=pdf(normal,-smoothbins:smoothbins);

Preference=NaN(length(Choice(:,session)),150);
for rat=1:length(Choice(:,session))
    Choices=Choice{rat,session}(:,3);
    numchoices(rat,1)=length(Choices);
    for l=1:length(Choices)
       Preference(rat,l)=sum(Choices(l-min([l-1 smoothbins]):l+min([length(Choices)-l smoothbins]),1)'.*fliplr(filterweights(1+max([0 smoothbins+l-length(Choices)]):smoothbins+min([l smoothbins+1]))))/sum(filterweights(1+max([0 smoothbins+l-length(Choices)]):smoothbins+min([l smoothbins+1])));                  
    end
end


figure;
n=0;
a={};
numtrials=min(numchoices(chr));
indtrials=max(numchoices);
for group=1:2
    subplot(1,3,3);
    hold on

    a{group}=plot(1:numtrials,nanmean(Preference(groups{group,1},1:numtrials)),'linewidth',2,'color',colors{group});
    up=nanmean(Preference(groups{group,1},1:numtrials))+nanste(Preference(groups{group,1},1:numtrials),1);
    down=nanmean(Preference(groups{group,1},1:numtrials))-nanste(Preference(groups{group,1},1:numtrials),1);
    patch([1:numtrials,numtrials:-1:1],[up,down(end:-1:1)],colors{group},'EdgeColor','none');alpha(0.3);

    plot([0 numtrials+1],[0.5 0.5],':','color','k');
    xlabel('Choice trials');
    axis([0 numtrials+1 0 1]);
    title('Group mean+/-SEM');
    yticks([0 0.5 1]);
    yticklabels([-1;0;1]);
    if  group==2 legend([a{:}],'Control','ChR2','location','southwest'); end
    
    subplot(1,3,1:2);
    hold on
    
    for i=1:rats
        if sum(groups{group,1}==i)==1
            n=n+1;
            plot(1:indtrials,Preference(i,1:indtrials),'linewidth',2,'color',colors{group});
        end
    end
    

    plot([0 indtrials+1],[0.5 0.5],':','color','k');
    xlabel('Choice trials');
    ylabel('Preference')
    yticks([0 0.5 1]);
    yticklabels([-1;0;1]);
    axis([0 indtrials 0 1]);
    title('Individual rats');
end

%% Latency for session of interest
rats=length(Choice(:,1));

session=13;
groups{1,1}=gfp;
groups{2,1}=chr;


SucLat=NaN(length(Choice(:,session)),4);
for rat=1:length(Choice(:,session))
    Latencies=log(suclatencies{rat,session});
    trials=1:length(Latencies);
    Break1=length(trials)/4;
    Break2=length(trials)/2;
    Break3=3*length(trials)/4;
    
    SucLat(rat,1) = mean(Latencies(trials<Break1));
    SucLat(rat,2) = mean(Latencies(trials>=Break1 & trials<Break2));
    SucLat(rat,3) = mean(Latencies(trials>=Break2 & trials<Break3));
    SucLat(rat,4) = mean(Latencies(trials>=Break3));
end

MalLat=NaN(length(Choice(:,session)),4);
for rat=1:length(Choice(:,session))
    Latencies=log(mallatencies{rat,session});
    trials=1:length(Latencies);
    Break1=length(trials)/4;
    Break2=length(trials)/2;
    Break3=3*length(trials)/4;
    
    MalLat(rat,1) = mean(Latencies(trials<Break1));
    MalLat(rat,2) = mean(Latencies(trials>=Break1 & trials<Break2));
    MalLat(rat,3) = mean(Latencies(trials>=Break2 & trials<Break3));
    MalLat(rat,4) = mean(Latencies(trials>=Break3));
end

LatDiff=MalLat - SucLat;


figure;
n=0;
a={};

for group=1:2
    subplot(1,3,3);
    hold on

    a{group}=plot(1:4,nanmean(LatDiff(groups{group,1},1:4)),'linewidth',2,'color',colors{group});
    up=nanmean(LatDiff(groups{group,1},1:4))+nanste(LatDiff(groups{group,1},1:4),1);
    down=nanmean(LatDiff(groups{group,1},1:4))-nanste(LatDiff(groups{group,1},1:4),1);
    patch([1:4,4:-1:1],[up,down(end:-1:1)],colors{group},'EdgeColor','none');alpha(0.3);

    plot([0 4+1],[0 0],':','color','k');
    xlabel('Quarter');

    axis([0 5 -3 3]);
    title('Group mean+/-SEM');

    if  group==2 legend([a{:}],'Control','ChR2','location','southwest'); end
    xticks(1:4);
    
    subplot(1,3,2);
    hold on
    
    for i=1:rats
        if sum(groups{group,1}==i)==1
            n=n+1;
            plot(1:4,LatDiff(i,1:4),'linewidth',2,'color',colors{group});
        end
    end
    

    plot([0 4+1],[0 0],':','color','k');
    xlabel('Quarter');

    ylabel({'Maltodextrin - sucrose log(s)'});

    axis([0.5 4.5 -3 4.5]);
    title({'Effect of stim on sucrose';'versus maltodextrin press latency'});
    xticks(1:4);
end


%stats
signrank(LatDiff(chr,1),LatDiff(chr,4));
signrank(LatDiff(gfp,1),LatDiff(gfp,4));
ranksum(LatDiff(chr,1),LatDiff(gfp,1));
ranksum(LatDiff(chr,4),LatDiff(gfp,4));

%latency to press post stimulation
LatencyPostStim=NaN(length(Choice(:,session)),2);
for rat=1:length(Choice(:,session))
    LatencyPostStim(rat,2)=mean(log(alllatencies{rat,session}(poststimtrls{rat,session}(1:length(alllatencies{rat,session}))==1)));
    LatencyPostStim(rat,1)=mean(log(alllatencies{rat,session}(poststimtrls{rat,session}(1:length(alllatencies{rat,session}))==0)));
end
LatDiffStim=LatencyPostStim(:,2)-LatencyPostStim(:,1);

subplot(1,3,1);
hold on;
scatter(rand([length(gfp) 1])/2+0.65,LatDiffStim(gfp),55,colors{1},'filled');
scatter(rand([length(chr) 1])/2+1.35,LatDiffStim(chr),55,colors{2},'filled');
title({'Effect of stim on';'next trial press latency'});
ylabel({'Post-stim - post-no-stim log(s)'});
plot([0 2.4],[0 0],':','color','k');
xticks([0.9 1.6]);
xticklabels({'GFP','Chr2'});
axis([0.4 2.1 -3 3]);

%stats
ranksum(LatDiffStim(chr),LatDiffStim(gfp));
%% self-stim
gfp=[4 6 10 11 13 14 15 21];
chr=[1 2 7 8 9 16 17 18 22 23 24];

%file location
address=['/Users/david/Documents/GitHub/dynamic-preference/Data/ICSSFiles'];
stim_files=dir([address,'//*12-12*']); %mac

for k=1:length(stim_files)
    filename=fullfile(stim_files(k).folder,stim_files(k).name);
    file=fopen(filename);
    L=textscan(file,'%s','Delimiter',':');
    fclose('all');
    x=str2num(files(k).name(30:31)); %x is rat #
    

    %find start of A, start of C, and end of C
    Nstrt=find(strcmp('N',L{1,1})); %important values
    
    stim_pokes(x,1) = str2num(L{1,1}{Nstrt+2});
    cont_pokes(x,1) = str2num(L{1,1}{Nstrt+4});
    stims(x,1) = str2num(L{1,1}{Nstrt+10});



end

all_pokes = [cont_pokes stim_pokes];

figure;
subplot(1,3,1);
hold on;
plot([1 2],all_pokes(gfp,:),'color',[0.4 0.4 0.4]);
errorbar(1,nanmean(all_pokes(gfp,1)),nanste(all_pokes(gfp,1),1),'color',[0.4 0.4 0.4],'marker','o','linewidth',1.5);
errorbar(2,nanmean(all_pokes(gfp,2)),nanste(all_pokes(gfp,2),1),'color',[0 0.3 1],'marker','o','linewidth',1.5);
ylabel('Nosepokes');
xticks([1 2]);
xticklabels({'No laser','Laser'});
xtickangle(45);
axis([0.5 2.5 0 100]);
title('YFP');

pval=signrank(all_pokes(gfp,1),all_pokes(gfp,2));
text(1.2,100,['p = ' num2str(round(pval,2,'significant'))]);


subplot(1,3,2);
hold on;
plot([1 2],all_pokes(chr,:),'color',[0.4 0.4 0.4]);
errorbar(1,nanmean(all_pokes(chr,1)),nanste(all_pokes(chr,1),1),'color',[0.4 0.4 0.4],'marker','o','markerfacecolor',[0.4 0.4 0.4],'linewidth',1.5);
errorbar(2,nanmean(all_pokes(chr,2)),nanste(all_pokes(chr,2),1),'color',[0 0.3 1],'marker','o','markerfacecolor',[0 0.3 1],'linewidth',1.5);
ylabel('Nosepokes');
xticks([1 2]);
xticklabels({'No laser','Laser'});
xtickangle(45);
axis([0.5 2.5 0 4500]);
title('ChR2');

pval=signrank(all_pokes(chr,1),all_pokes(chr,2));
text(1.2,375,['p = ' num2str(round(pval,2,'significant'))]);


subplot(1,3,3)
hold on;
scatter(rand([length(gfp) 1])/2+0.65,stims(gfp),55,colors{2});
scatter(rand([length(chr) 1])/2+1.35,stims(chr),55,colors{2},'filled');

axis([0.5 2 0 2500]);pval=ranksum(stims(gfp,1),stims(chr,1));
text(1.2,375,['p = ' num2str(round(pval,2,'significant'))]);
xticks([1 1.5]);
xticklabels({'YFP','ChR2'});
ylabel('Stimulations earned');
xtickangle(45);
