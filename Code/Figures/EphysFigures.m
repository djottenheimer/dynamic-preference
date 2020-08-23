%%
load('R_choice.mat')
load('RAWCh.mat')

%% highlighting cue and reward
predictable=[0.9 0.9 0.9];
unpredictable=[0.1 0.1 0.1];
figure;

props=[];
for trial=1 %forced or choice

    
    Evs={'Cue','Reward'};

    %Outcome
    subplot(2,3,1+(trial-1)*3)
    hold on
    props(1:2,1)=(sum(R.GLM{trial,1}.PVal.Outcome(:,2:4:end)<0.05&R.Typ==0 ,1)/sum(R.Typ==0))'; %task 1
    props(1:2,2)=(sum(R.GLM{trial,1}.PVal.Outcome(:,2:4:end)<0.05&R.Typ,1)/sum(R.Typ))'; %task 2
    b=bar(props); 
    b(1).FaceColor=predictable;
    b(2).FaceColor=unpredictable;
    

    %stats
    pval=[];
    chi2=[];
    for i=1:2
        [~,chi2(1,i),pval(1,i)]=crosstab(R.GLM{trial,1}.PVal.Outcome(:,2+4*(i-1))<0.05,R.Typ);
    end    
    plot(1:2,(pval(1,:)<0.05)-0.4,'*','color','k');
    
    xticks(1:2);
    xticklabels(Evs);
    xtickangle(45);
    axis([0.5 2.5 0 1])
    legend('Predictive task','Non-predictive task')
    title('Outcome')
    if trial==1 ylabel({'Forced trials';'Fraction of the population'}); end
    if trial==2 ylabel({'Choice trials';'Fraction of the population'}); end

    %Preference
    subplot(2,3,2+(trial-1)*3)
    hold on
    props(1:2,1)=(sum(R.GLM{trial,1}.PVal.Preference(:,2:4:end)<0.05&R.Typ==0,1)/sum(R.Typ==0))'; %task 1
    props(1:2,2)=(sum(R.GLM{trial,1}.PVal.Preference(:,2:4:end)<0.05&R.Typ,1)/sum(R.Typ))'; %task 2
    b=bar(props); 
    b(1).FaceColor=predictable;
    b(2).FaceColor=unpredictable;
    
    %stats
    for i=1:2
        [~,chi2(2,i),pval(2,i)]=crosstab(R.GLM{trial,1}.PVal.Preference(:,2+4*(i-1))<0.05,R.Typ);
    end    
    plot(1:2,(pval(2,:)<0.05)-0.3,'*','color','k');
    
    
    xticks(1:6);
    xticklabels(Evs);
    xtickangle(45);
    axis([0.5 2.5 0 1])
    title('Time')

    %Outcome by preference
    subplot(2,3,3+(trial-1)*3)
    hold on
    props(1:2,1)=(sum(R.GLM{trial,1}.PVal.OutByPref(:,2:4:end)<0.05&R.Typ==0&R.GLM{trial,1}.Coefficients.OutByPref(:,2:4:end)>0,1)/sum(R.Typ==0))'; %task 1
    props(1:2,2)=(sum(R.GLM{trial,1}.PVal.OutByPref(:,2:4:end)<0.05&R.Typ&R.GLM{trial,1}.Coefficients.OutByPref(:,2:4:end)>0,1)/sum(R.Typ))'; %task 2
    b=bar(props); 
    b(1).FaceColor=predictable;
    b(2).FaceColor=unpredictable;
    
    %stats
    for i=1:2
        [~,chi2(3,i),pval(3,i)]=crosstab(R.GLM{trial,1}.PVal.OutByPref(:,2+4*(i-1))<0.05&R.GLM{trial,1}.Coefficients.OutByPref(:,2+4*(i-1))>0,R.Typ);
    end    
    plot(1:2,(pval(3,:)<0.05)-0.3,'*','color','k');
    
    
    xticks(1:2);
    xticklabels(Evs);
    xtickangle(45);
    axis([0.5 2.5 0 1])
    title('Out X Time')
end

%% individual rat points

predictable=[0.9 0.9 0.9];
unpredictable=[0.1 0.1 0.1];
figure;

rats=unique(R.Rat);
ratprops={};
for trial=1 %forced or choice
    
    Evs={'Cue','Reward'};

    %Outcome
    subplot(2,3,1+(trial-1)*3)
    hold on
    for rat=1:length(rats)
        Sel=R.Rat==rats(rat);
        for task=0:1
            ratprops{task+1}(rat,1:2) = sum(R.GLM{trial,1}.PVal.Outcome(:,2:4:end)<0.05&R.Typ==task&Sel ,1)/sum(R.Typ==task & Sel);
        end
    end
    
    plot(1:2,[ratprops{1}(:,1) ratprops{2}(:,1)],'color',[0 0 0]);
    plot(3:4,[ratprops{1}(:,2) ratprops{2}(:,2)],'color',[0 0 0]);
    axis([0 5 0 1]);
    xticks([1 2 3 4]);
    xticklabels({'S.C.','U.O.','S.C.','U.O'});
    text(1,0.8,'Cue');
    text(3,0.8,'Reward');
    title('Outcome');
    

    %stats
    [~,cue_pval]=ttest(ratprops{1}(:,1),ratprops{2}(:,1));
    text(1,0.9,['p=' num2str(round(cue_pval,2,'significant'))]);
    [~,reward_pval]=ttest(ratprops{1}(:,2),ratprops{2}(:,2));
    text(3,0.9,['p=' num2str(round(reward_pval,2,'significant'))]);


    %Preference
    subplot(2,3,2+(trial-1)*3)
    hold on
    for rat=1:length(rats)
        Sel=R.Rat==rats(rat);
        for task=0:1
            ratprops{task+1}(rat,1:2) = sum(R.GLM{trial,1}.PVal.Preference(:,2:4:end)<0.05&R.Typ==task&Sel ,1)/sum(R.Typ==task & Sel);
        end
    end
    
    plot(1:2,[ratprops{1}(:,1) ratprops{2}(:,1)],'color',[0 0 0]);
    plot(3:4,[ratprops{1}(:,2) ratprops{2}(:,2)],'color',[0 0 0]);
    axis([0 5 0 1]);
    xticks([1 2 3 4]);
    xticklabels({'S.C.','U.O.','S.C.','U.O'});
    text(1,0.8,'Cue');
    text(3,0.8,'Reward');
    title('Preference');
    

    %stats
    [~,cue_pval]=ttest(ratprops{1}(:,1),ratprops{2}(:,1));
    text(1,0.9,['p=' num2str(round(cue_pval,2,'significant'))]);
    [~,reward_pval]=ttest(ratprops{1}(:,2),ratprops{2}(:,2));
    text(3,0.9,['p=' num2str(round(reward_pval,2,'significant'))]);

    %Outcome by preference
    subplot(2,3,3+(trial-1)*3)
    hold on
    for rat=1:length(rats)
        Sel=R.Rat==rats(rat);
        for task=0:1
            ratprops{task+1}(rat,1:2) = sum(R.GLM{trial,1}.PVal.OutByPref(:,2:4:end)<0.05&R.GLM{trial,1}.Coefficients.OutByPref(:,2:4:end)>0&R.Typ==task&Sel ,1)/sum(R.Typ==task & Sel);
        end
    end
    
    plot(1:2,[ratprops{1}(:,1) ratprops{2}(:,1)],'color',[0 0 0]);
    plot(3:4,[ratprops{1}(:,2) ratprops{2}(:,2)],'color',[0 0 0]);
    axis([0 5 0 1]);
    xticks([1 2 3 4]);
    xticklabels({'S.C.','U.O.','S.C.','U.O'});
    text(1,0.8,'Cue');
    text(3,0.8,'Reward');
    title('Out X Time');
    

    %stats
    [~,cue_pval]=ttest(ratprops{1}(:,1),ratprops{2}(:,1));
    text(1,0.9,['p=' num2str(round(cue_pval,2,'significant'))]);
    [~,reward_pval]=ttest(ratprops{1}(:,2),ratprops{2}(:,2));
    text(3,0.9,['p=' num2str(round(reward_pval,2,'significant'))]);
end

%% R2
figure;
WhichNeurons=R.GLM{1,1}.PVal.OutByPref(:,2)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,2)>0;
subplot(1,2,1);
histogram(R.GLM{1,1}.R2(WhichNeurons&R.Typ==0,2),0:0.05:1);
title('Specific Cues, Out X Time at Cue');
xlabel('R2');
ylabel('# neurons');
text(0.6,6,['Median = ' num2str(round(median(R.GLM{1,1}.R2(WhichNeurons&R.Typ==0,2)),2))]);

WhichNeurons=R.GLM{1,1}.PVal.OutByPref(:,6)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,6)>0;
subplot(1,2,2);
histogram(R.GLM{1,1}.R2(WhichNeurons&R.Typ==1,6),0:0.05:1);
title('Uncertain Outcome, Out X Time at Reward');
xlabel('R2');
ylabel('# neurons');
text(0.75,15,['Median = ' num2str(round(median(R.GLM{1,1}.R2(WhichNeurons&R.Typ==1,6)),2))]);

%% behavior over time, all sessions

figure;
Sel=cat(1,RAW.Ses)>0;

%setting up parameters
Xaxis=[-3 13];
inttime=find(R.Tm>=Xaxis(1) & R.Tm<=Xaxis(2)); %+2 because it's aligned to lever press/port entry
PSTHtime=R.Tm(inttime); %+2 because it's aligned to lever press/port entry

Colors={[0.2 0.1 0],[0 0 0.1];
        [0.4 0.1 0],[0.2 0.1 0.4];
        [0.8 0.3 0],[0.3 0.2 0.7];
        [1 0.5 0],[0.5 0.3 0.9]};
titles={'Sucrose licking, forced','Water licking, forced'};
for i=1:2

    subplot(2,2,i)
    hold on

    for j=1:4
        %plot each quadrant
        psth1=nanmean(R.Lick((i-1)*4+j).PSTHraw(Sel,inttime),1); 
        sem1=nanste(R.Lick((i-1)*4+j).PSTHraw(Sel,inttime),1); %calculate standard error of the mean
        upE=psth1+sem1;
        downE=psth1-sem1;


        %plotting
        p{j}=plot(PSTHtime,psth1,'Color',Colors{j,i},'linewidth',1);
        patch([PSTHtime,PSTHtime(end:-1:1)],[upE,downE(end:-1:1)],Colors{j,i},'EdgeColor','none');alpha(0.5);
    end
    
    
    plot([0 0],[-3 8],':','color','k','linewidth',0.75);
    axis([Xaxis(1) Xaxis(2) 0 8]);
    ylabel('Licks/s');
    xlabel('Seconds from reward delivery');
    legend([p{:}],'Q1','Q2','Q3','Q4');
    title(titles{i})


end

titles={'Sucrose port occupancy','Water port occupancy'};
for i=1:2

    subplot(2,2,i+2)
    hold on

    for j=1:4
        %plot each quadrant
        psth1=nanmean(R.Port((i-1)*4+j).PSTHraw(Sel,inttime),1); 
        sem1=nanste(R.Port((i-1)*4+j).PSTHraw(Sel,inttime),1); %calculate standard error of the mean
        upE=psth1+sem1;
        downE=psth1-sem1;


        %plotting
        p{j}=plot(PSTHtime,psth1,'Color',Colors{j,i},'linewidth',1);
        patch([PSTHtime,PSTHtime(end:-1:1)],[upE,downE(end:-1:1)],Colors{j,i},'EdgeColor','none');alpha(0.5);
    end

    plot([-4 5],[0 0],':','color','k','linewidth',0.75);
    plot([0 0],[-3 8],':','color','k','linewidth',0.75);
    axis([Xaxis(1) Xaxis(2) 0 1.03]);
    ylabel('Probability');
    xlabel('Seconds from reward delivery');
    legend([p{:}],'Q1','Q2','Q3','Q4');
    title(titles{i})


end

%% example session
session=2;

sucrose=[0.8 0.3 0];
water=[0.3 0.2 0.7];
choice=[0.1 0.6 0.4];
figure;


subplot(1,1,1)
hold on;

tickwidth=1;

pl={};
choices=find(R.TrialInfo{session,1}.Choice==1&R.TrialInfo{session,1}.Outcome==1);
for i=1:length(choices)
    pl{1}=plot([choices(i) choices(i)],[0.6 1],'color',sucrose,'linewidth',tickwidth);
end

choices=find(R.TrialInfo{session,1}.Choice==1&R.TrialInfo{session,1}.Outcome==0);
for i=1:length(choices)
    pl{2}=plot([choices(i) choices(i)],[-1 -0.6],'color',water,'linewidth',tickwidth);
end



forced=find(R.TrialInfo{session,1}.Choice==0&R.TrialInfo{session,1}.Outcome==0);
for i=1:length(forced)
    pl{3}=plot([forced(i) forced(i)],[-1 -.8],'color','k','linewidth',tickwidth);
end

forced=find(R.TrialInfo{session,1}.Choice==0&R.TrialInfo{session,1}.Outcome==1);
for i=1:length(forced)
    plot([forced(i) forced(i)],[0.8 1],'color','k','linewidth',tickwidth);
end

plot([0 length(R.TrialInfo{session,1}.Preference)],[0 0],':','color','k','linewidth',0.75);
pl{4}=plot(R.TrialInfo{session,1}.Preference*2-1,'color',choice,'linewidth',2);

legend([pl{:}],'Sucrose choice','Water choice','Forced trial','Smoothed preference','location','southeast');

yticks([-1 -0.5 0 0.5 1]);
ylabel('Preference');
xlabel('Trial');
title('Example session');

%% quadrant summary for behavior

%choice behavior
sections=4;
Preference=[];
for i=1:length(RAW)
    Trials=R.TrialInfo{i,1}.TrialNo;
    FS=R.TrialInfo{i,1}.Outcome==1 & R.TrialInfo{i,1}.Choice==0 & R.TrialInfo{i,1}.InPort;
    FW=R.TrialInfo{i,1}.Outcome==0 & R.TrialInfo{i,1}.Choice==0 & R.TrialInfo{i,1}.InPort;
    Ch=R.TrialInfo{i,1}.Choice==1;

    breakpoint(1)=0;
    section={};
    for j=1:sections
        breakpoint(j+1)=j*Trials(end)/sections;
        section{j}=Trials>breakpoint(j) & Trials<=breakpoint(j+1);
        Preference(i,j)=mean(R.TrialInfo{i,1}.Outcome(section{j}&Ch));
    end
end

%normalize it from -1 to 1
Preference=Preference*2-1;

figure;
for task=1:2
    %Sel=R.GLM{1,1}.PVal.OutByPref(:,6)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,6)>0;
    if task==1 Sel=cat(1,RAW.Ses)<10; end
    if task==2 Sel=cat(1,RAW.Ses)>10; end

    %setting up parameters
    Colors={[0.2 0.1 0],[0 0 0.1],[0 0.2 0];
            [0.4 0.1 0],[0.2 0.1 0.4],[0 0.4 0.3];
            [0.8 0.3 0],[0.3 0.2 0.7],[0.1 0.6 0.4];
            [1 0.5 0],[0.5 0.3 0.9],[0.2 0.8 0.6]};
    titles={'Average preference','Individual preference'};
    



        
    subplot(2,2,1+(task-1)*2)
    hold on

    prefup=mean(Preference(Sel,:))+nanste(Preference(Sel,:),1);
    prefdown=mean(Preference(Sel,:))-nanste(Preference(Sel,:),1);

    plot(1:4,mean(Preference(Sel,:)),'Color',Colors{3,3},'linewidth',1);
    patch([1:4,4:-1:1],[prefup,prefdown(end:-1:1)],Colors{3,3},'EdgeColor','none');alpha(0.5);
    plot([0 5],[0 0],':','color','k','linewidth',0.75);

    axis([0 5 -1 1]);
    if task==2 xlabel('Quadrant'); end
    if task==1 title(titles{1}); end
    ylabel('Fraction sucrose choices');
    xticks(1:4);

    subplot(2,2,2+(task-1)*2)
    hold on

    plot(1:4,Preference(Sel,:),'Color',Colors{3,3},'linewidth',1);
    plot([0 5],[0 0],':','color','k','linewidth',0.75);
    
    if task==1 title(titles{2}); end
    if task==2 xlabel('Quadrant'); end
    axis([0 5 -1 1]);
    xticks(1:4);

   
    
    

end

%combined behavior
figure;

for window=1:2
    if window==1 bin=[0 13]; end
    if window==2 bin=[0.75 1.95]; end

    %setting up parameters
    Colors={[0.2 0.1 0],[0 0 0.1],[0 0.2 0];
            [0.4 0.1 0],[0.2 0.1 0.4],[0 0.4 0.3];
            [0.8 0.3 0],[0.3 0.2 0.7],[0.1 0.6 0.4];
            [1 0.5 0],[0.5 0.3 0.9],[0.2 0.8 0.6]};
    titles={'Licks','Port Occupancy'};
    
    
    %time periods we're plotting
    toi=find(R.Tm>=bin(1) & R.Tm<=bin(2));
    a=[];
    for i=1:2 %which reward



        for j=1:4 %which quadrant
            %get mean for each quadrant           
            lickmean(i,j)=nanmean(nanmean(R.Lick((i-1)*4+j).PSTHnosmooth(Sel,toi),2));
            lickste(i,j)=nanste(nanmean(R.Lick((i-1)*4+j).PSTHnosmooth(Sel,toi),2),1);
            pemean(i,j)=nanmean(nanmean(R.Port((i-1)*4+j).PSTHraw(Sel,toi),2));
            peste(i,j)=nanste(nanmean(R.Port((i-1)*4+j).PSTHraw(Sel,toi),2),1);
        end
        
        subplot(2,2,1+(window-1)*2)
        hold on
        
        lickup=lickmean(i,:)+lickste(i,:);
        lickdown=lickmean(i,:)-lickste(i,:);
        
        a{i}=plot(1:4,lickmean(i,:),'Color',Colors{3,i},'linewidth',1);
        patch([1:4,4:-1:1],[lickup,lickdown(end:-1:1)],Colors{3,i},'EdgeColor','none');alpha(0.5);
        ylabel({[num2str(bin)];'Licks/s'});
        if window==2 xlabel('Quadrant'); end
        if window==1 title(titles{1}); end
        if window==1 legend([a{:}],'Sucrose','Water'); end
        axis([0 5 0 4+window*2]);
        xticks(1:4);


        subplot(2,2,2+(window-1)*2)
        hold on
        
        peup=pemean(i,:)+peste(i,:);
        pedown=pemean(i,:)-peste(i,:);
        
        plot(1:4,pemean(i,:),'Color',Colors{3,i},'linewidth',1);
        patch([1:4,4:-1:1],[peup,pedown(end:-1:1)],Colors{3,i},'EdgeColor','none');alpha(0.5);    
        if window==2 xlabel('Quadrant'); end
        if window==1 title(titles{2}); end
        ylabel('Probability');
        axis([0 5 0 1]);
        xticks(1:4);

        
        
    end
    
    

end


%% example neuron raster -- cue
figure;

WhichNeurons=R.GLM{1,1}.PVal.OutByPref(:,2)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,2)>0 & R.Typ==0;
NeuronList=find(WhichNeurons);

alph=0.3;
%raster
window=[-0.5 1.5];
sample_neuron_number=43; 
AnalysisWindow=[0 0.75];

sessrat=R.Rat(sample_neuron_number);
sessneurons=find(R.Rat==sessrat&R.Typ==0);
sessneuron=sessneurons==sample_neuron_number;
session=find(unique(R.Rat)==sessrat)*2-1;

subplot(1,2,1);
hold on;


RD=strcmp('CueS', RAW(session).Einfo(:,2));
RDtimes=RAW(session).Erast{RD,1};

patch([AnalysisWindow(1) AnalysisWindow(1) AnalysisWindow(2) AnalysisWindow(2)],[0 length(RDtimes)*11 length(RDtimes)*11 0],[.8 .3 0],'edgecolor','none');alpha(alph);
PlotRaster(RAW(session).Nrast{sessneuron,1},RDtimes,window);
ylabel('Trials, sorted by session progress');
yticks([]);
xticks(-1:4);
xticklabels(-1:2);
xlabel('Seconds from cue onset');
xlim(window);
title('Sucrose cue');

subplot(1,2,2);
hold on;

RD=strcmp('CueW', RAW(session).Einfo(:,2));
RDtimes=RAW(session).Erast{RD,1};

patch([AnalysisWindow(1) AnalysisWindow(1) AnalysisWindow(2) AnalysisWindow(2)],[0 length(RDtimes)*11 length(RDtimes)*11 0],[.3 .2 .7],'edgecolor','none');alpha(alph);
PlotRaster(RAW(session).Nrast{sessneuron,1},RDtimes,window);
ylabel('');
yticks([]);
xticks(-1:2);
xticklabels(-1:2);
xlabel('');
xlim(window);
title('Water cue');

%% firing over time -- cues

WhichNeurons=R.GLM{1,1}.PVal.OutByPref(:,2)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,2)>0;
%WhichNeurons=WhichNeurons==0; %to get the other neurons
AnalysisWindow=[0 0.75];
alph=0.3;

figure;
for task=1

    if task==1 Sel=WhichNeurons & R.Typ==0; end
    if task==2 Sel=WhichNeurons & R.Typ; end

    %setting up parameters
    Xaxis=[-1 2];
    inttime=find(R.Tm>=Xaxis(1) & R.Tm<=Xaxis(2));
    PSTHtime=R.Tm(inttime); 

    Colors={[0.2 0.1 0],[0 0 0.1],[0 0.2 0];
            [0.4 0.1 0],[0.2 0.1 0.4],[0 0.4 0.3];
            [0.8 0.3 0],[0.3 0.2 0.7],[0.1 0.6 0.4];
            [1 0.5 0],[0.5 0.3 0.9],[0.2 0.8 0.6]};
    titles={'Sucrose cue','Water cue','Choice cue'};
    for i=1:2

        subplot(2,3,i+(task-1)*3)
        hold on
        
        patch([AnalysisWindow(1) AnalysisWindow(1) AnalysisWindow(2) AnalysisWindow(2)],[-3 11 11 -3],[.8 .8 .8],'edgecolor','none');alpha(alph);

        pl={};
        for j=1:4
            %plot each quadrant
            psth1=nanmean(R.Ev(8+(i-1)*4+j).PSTHz(Sel,inttime),1); 
            sem1=nanste(R.Ev(8+(i-1)*4+j).PSTHz(Sel,inttime),1); %calculate standard error of the mean
            upE=psth1+sem1;
            downE=psth1-sem1;


            %plotting
            pl{j}=plot(PSTHtime,psth1,'Color',Colors{j,i},'linewidth',1);
            patch([PSTHtime,PSTHtime(end:-1:1)],[upE,downE(end:-1:1)],Colors{j,i},'EdgeColor','none');alpha(0.5);
        end

        plot([-4 5],[0 0],':','color','k','linewidth',0.75);
        plot([0 0],[-3 11],':','color','k','linewidth',0.75);
        
        
        axis([Xaxis(1) Xaxis(2) -3 11]);
        ylabel('Mean firing (z-score)');
        xlabel('Seconds from cue onset');
        legend([pl{:}],'Q1','Q2','Q3','Q4');
        title(titles{i})
        
        
    end

end



% quadrant summary for cue
axes=[0 5 -2 4];
a=[];
for task=1
    %Sel=R.GLM{1,1}.PVal.OutByPref(:,6)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,6)>0;
    if task==1 Sel=WhichNeurons & R.Typ==0; end
    if task==2 Sel=WhichNeurons & R.Typ; end

    %setting up parameters
    Colors={[0.2 0.1 0],[0 0 0.1],[0 0.2 0];
            [0.4 0.1 0],[0.2 0.1 0.4],[0 0.4 0.3];
            [0.8 0.3 0],[0.3 0.2 0.7],[0.1 0.6 0.4];
            [1 0.5 0],[0.5 0.3 0.9],[0.2 0.8 0.6]};
    titles={'Binned firing','Action','Reward'};
    
    %time periods we're plotting
    cuebin=[0 0.75];
    cuetoi=find(R.Tm>=cuebin(1) & R.Tm<=cuebin(2));
 
    for i=1:2 %which reward



        for j=1:4 %which quadrant
            %get mean for each quadrant
            cuemean(i,j)=nanmean(nanmean(R.Ev(8+(i-1)*4+j).PSTHnosmoothz(Sel,cuetoi),2));
            cueste(i,j)=nanste(nanmean(R.Ev(8+(i-1)*4+j).PSTHnosmoothz(Sel,cuetoi),2),1);
        end
        
        subplot(2,3,3)
        hold on
        
        cueup=cuemean(i,:)+cueste(i,:);
        cuedown=cuemean(i,:)-cueste(i,:);
        
        a{i}=plot(1:4,cuemean(i,:),'Color',Colors{3,i},'linewidth',1);
        patch([1:4,4:-1:1],[cueup,cuedown(end:-1:1)],Colors{3,i},'EdgeColor','none');alpha(0.5);
        if task==1 ylabel({'Population mean (z-score)'}); end
        if task==2 ylabel({'Unpredictable';'Population mean (z-score)'}); end
        if task==1 xlabel('Quarter'); end
        if task==1 title(titles{1}); end
        
        axis(axes);
        plot([0 5],[0 0],':','color','k','linewidth',0.75);
        xticks(1:4);
        if task==1 legend([a{:}],'Sucrose','Water'); end

        
    end

    

end

%% firing over time for time neurons -- cues

WhichNeurons=R.GLM{1,1}.PVal.Preference(:,2)<0.05 & R.GLM{1,1}.Coefficients.Preference(:,2)<0;

AnalysisWindow=[0 0.75];
alph=0.3;

figure;
for task=1:2

    if task==1 Sel=WhichNeurons & R.Typ==0; end
    if task==2 Sel=WhichNeurons & R.Typ; end

    %setting up parameters
    Xaxis=[-1 2];
    inttime=find(R.Tm>=Xaxis(1) & R.Tm<=Xaxis(2));
    PSTHtime=R.Tm(inttime); 

    Colors={[0.2 0.1 0],[0 0 0.1],[0 0.2 0];
            [0.4 0.1 0],[0.2 0.1 0.4],[0 0.4 0.3];
            [0.8 0.3 0],[0.3 0.2 0.7],[0.1 0.6 0.4];
            [1 0.5 0],[0.5 0.3 0.9],[0.2 0.8 0.6]};
    titles={'Sucrose cue','Water cue','Choice cue'};
    for i=1:3

        subplot(2,4,i+(task-1)*4)
        hold on
        
        patch([AnalysisWindow(1) AnalysisWindow(1) AnalysisWindow(2) AnalysisWindow(2)],[-3 11 11 -3],[.8 .8 .8],'edgecolor','none');alpha(alph);

        pl={};
        for j=1:4
            %plot each quadrant
            psth1=nanmean(R.Ev(8+(i-1)*4+j).PSTHz(Sel,inttime),1); 
            sem1=nanste(R.Ev(8+(i-1)*4+j).PSTHz(Sel,inttime),1); %calculate standard error of the mean
            upE=psth1+sem1;
            downE=psth1-sem1;


            %plotting
            pl{j}=plot(PSTHtime,psth1,'Color',Colors{j,i},'linewidth',1);
            patch([PSTHtime,PSTHtime(end:-1:1)],[upE,downE(end:-1:1)],Colors{j,i},'EdgeColor','none');alpha(0.5);
        end

        plot([-4 5],[0 0],':','color','k','linewidth',0.75);
        plot([0 0],[-3 11],':','color','k','linewidth',0.75);
        
        
        axis([Xaxis(1) Xaxis(2) -3 11]);
        ylabel('Mean firing (z-score)');
        xlabel('Seconds from cue onset');
        legend([pl{:}],'Q1','Q2','Q3','Q4');
        title(titles{i})
        
        
    end

end



% quadrant summary for cue
axes=[0 5 -2 4];
a=[];
for task=1:2
    %Sel=R.GLM{1,1}.PVal.OutByPref(:,6)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,6)>0;
    if task==1 Sel=WhichNeurons & R.Typ==0; end
    if task==2 Sel=WhichNeurons & R.Typ; end

    %setting up parameters
    Colors={[0.2 0.1 0],[0 0 0.1],[0 0.2 0];
            [0.4 0.1 0],[0.2 0.1 0.4],[0 0.4 0.3];
            [0.8 0.3 0],[0.3 0.2 0.7],[0.1 0.6 0.4];
            [1 0.5 0],[0.5 0.3 0.9],[0.2 0.8 0.6]};
    titles={'Binned firing','Action','Reward'};
    
    %time periods we're plotting
    cuebin=[0 0.75];
    cuetoi=find(R.Tm>=cuebin(1) & R.Tm<=cuebin(2));
    actbin=[-2 -1.25]; %relative to lever press
    acttoi=find(R.Tm>=actbin(1) & R.Tm<=actbin(2));    
    rewbin=[0 2.75]; %relative to lever press
    rewtoi=find(R.Tm>=rewbin(1) & R.Tm<=rewbin(2));    
    for i=1:3 %which reward



        for j=1:4 %which quadrant
            %get mean for each quadrant
            cuemean(i,j)=nanmean(nanmean(R.Ev(8+(i-1)*4+j).PSTHnosmoothz(Sel,cuetoi),2));
            cueste(i,j)=nanste(nanmean(R.Ev(8+(i-1)*4+j).PSTHnosmoothz(Sel,cuetoi),2),1);
        end
        
        subplot(2,4,4+(task-1)*4)
        hold on
        
        cueup=cuemean(i,:)+cueste(i,:);
        cuedown=cuemean(i,:)-cueste(i,:);
        
        a{i}=plot(1:4,cuemean(i,:),'Color',Colors{3,i},'linewidth',1);
        patch([1:4,4:-1:1],[cueup,cuedown(end:-1:1)],Colors{3,i},'EdgeColor','none');alpha(0.5);
        if task==1 ylabel({'Population mean (z-score)'}); end
        if task==2 ylabel({'Unpredictable';'Population mean (z-score)'}); end
        if task==1 xlabel('Quadrant'); end
        if task==1 title(titles{1}); end
        
        axis(axes);
        plot([0 5],[0 0],':','color','k','linewidth',0.75);
        xticks(1:4);
        if task==1 legend([a{:}],'Sucrose','Water'); end

        
    end

    

end

%% latencies
figure;
WhichNeurons=R.GLM{1,1}.PVal.OutByPref(:,2)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,2)>0;
included_sessions = find(cat(1,RAW.Ses)<10);
ExampleRat = 2;

Colors={[0.2 0.1 0],[0 0 0.1],[0 0.2 0];
        [0.4 0.1 0],[0.2 0.1 0.4],[0 0.4 0.3];
        [0.8 0.3 0],[0.3 0.2 0.7],[0.1 0.6 0.4];
        [1 0.5 0],[0.5 0.3 0.9],[0.2 0.8 0.6]};
        
NN = 0;
a={};
b={};
pval=[];
rho=[];
for sess_num = 1:length(included_sessions)
    session = included_sessions(sess_num);
    Cues = strmatch('CueAll',RAW(session).Einfo(:,2),'exact');
    CueTimes = RAW(session).Erast{Cues,1};
    Presses = strmatch('LeverR',RAW(session).Einfo(:,2),'exact');
    PressTimes = RAW(session).Erast{Presses,1};
    
    Latency = PressTimes - CueTimes(1:length(PressTimes));
    Latencies{sess_num,1} = Latency;
    TrialNo = 1:length(Latency);
    
    SucCue = R.TrialInfo{session,1}.Outcome == 1 & R.TrialInfo{session,1}.Choice == 0;
    WatCue = R.TrialInfo{session,1}.Outcome == 0 & R.TrialInfo{session,1}.Choice == 0;
    
    Break1=length(Latency)/4;
    Break2=length(Latency)/2;
    Break3=3*length(Latency)/4;
    
    SucQuads(sess_num,1) = mean(log(Latency(TrialNo<Break1 & SucCue')));
    SucQuads(sess_num,2) = mean(log(Latency(TrialNo>=Break1 & TrialNo<Break2 & SucCue')));
    SucQuads(sess_num,3) = mean(log(Latency(TrialNo>=Break2 & TrialNo<Break3 & SucCue')));
    SucQuads(sess_num,4) = mean(log(Latency(TrialNo>=Break3 & SucCue')));

    WatQuads(sess_num,1) = mean(log(Latency(TrialNo<Break1 & WatCue')));
    WatQuads(sess_num,2) = mean(log(Latency(TrialNo>=Break1 & TrialNo<Break2 & WatCue')));
    WatQuads(sess_num,3) = mean(log(Latency(TrialNo>=Break2 & TrialNo<Break3 & WatCue')));
    WatQuads(sess_num,4) = mean(log(Latency(TrialNo>=Break3 & WatCue')));    
    

    if sess_num == ExampleRat
        subplot(2,3,1);
        hold on;
        scatter(TrialNo(SucCue),log(Latency(SucCue)),12,Colors{3,1});
        axis([0 TrialNo(end) -2 7]);
        ylabel('log(Latency)');
        xlabel('Trial no.');
        title('Forced sucrose trials');
        xticks(20:20:200);
        
        %add a line :(
        Fit = polyfit(TrialNo(SucCue),log(Latency(SucCue))',1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
        plot(TrialNo(SucCue),polyval(Fit,TrialNo(SucCue)),'color','k')

        subplot(2,3,2);
        hold on;
        scatter(TrialNo(WatCue),log(Latency(WatCue)),12,Colors{3,2});    
        axis([0 TrialNo(end) -2 7]);
        xlabel('Trial no.');
        title('Forced water trials');
        xticks(20:20:200);
        
        %add a line :(
        Fit = polyfit(TrialNo(WatCue),log(Latency(WatCue))',1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
        plot(TrialNo(WatCue),polyval(Fit,TrialNo(WatCue)),'color','k')
        
    end
    
    %correlation between spikes and latency
    for neuron = 1:length(R.Hz{session,2}(1,:))
        NN = NN+1;
        [rho(NN,1),pval(NN,1)] = corr(Latency,R.Hz{session,2}(:,neuron),'type','spearman');      
        
    end
end


subplot(2,3,3)
hold on

sucmean=mean(SucQuads,1);
sucste=nanste(SucQuads,1);
sucup=sucmean+sucste;
sucdown=sucmean-sucste;

watmean=mean(WatQuads,1);
watste=nanste(WatQuads,1);
watup=watmean+watste;
watdown=watmean-watste;

a{1} = plot(1:4,sucmean,'Color',Colors{3,1},'linewidth',1);
patch([1:4,4:-1:1],[sucup,sucdown(end:-1:1)],Colors{3,1},'EdgeColor','none');alpha(0.5);

a{2} = plot(1:4,watmean,'Color',Colors{3,2},'linewidth',1);
patch([1:4,4:-1:1],[watup,watdown(end:-1:1)],Colors{3,2},'EdgeColor','none');alpha(0.5);

xticks(1:4);
xlabel('Quadrants');
ylabel('log(Latency)');
legend([a{:}],'Sucrose','Water','location','southeast');
axis([0 5 0 4]);

%example neuron
neuron=40;

SucCue = R.TrialInfo{1,1}.Outcome == 1 & R.TrialInfo{1,1}.Choice == 0;
WatCue = R.TrialInfo{1,1}.Outcome == 0 & R.TrialInfo{1,1}.Choice == 0;
ChoiceCue = R.TrialInfo{1,1}.Choice == 1;

subplot(2,3,5);
hold on;
scatter(R.Hz{1,2}(SucCue,neuron),log(Latencies{1,1}(SucCue)),12,Colors{3,1});
scatter(R.Hz{1,2}(WatCue,neuron),log(Latencies{1,1}(WatCue)),12,Colors{3,2});
scatter(R.Hz{1,2}(ChoiceCue,neuron),log(Latencies{1,1}(ChoiceCue)),12,Colors{3,3});

xlabel('Firing rate (Hz)');
ylabel('log(Latency)');
title('Example neuron');
text(12,4.5,['\rho = ' num2str(round(rho(neuron),2))]);
text(12,3.5,['p = ' num2str(round(pval(neuron),8))]);

%summary of correlation coefficients
subplot(2,3,6);
hold on;

plot([-1 1],[0.5 0.5],'color','k','linewidth',1);
plot([0 0],[0 01],'color','k','linewidth',1);

[cdf,x] = ecdf(rho(WhichNeurons(R.Typ==0)));
b{1} = plot(x,cdf,'linewidth',1.5,'color',[0 0.6 0.8]);
plot([mean(rho(WhichNeurons(R.Typ==0))) mean(rho(WhichNeurons(R.Typ==0)))],[0 1],'color',[0 0.6 0.8],'linewidth',1.5);

[cdf,x] = ecdf(rho(WhichNeurons(R.Typ==0)==0));
b{2} = plot(x,cdf,'linewidth',1.5,'color',[0.6 0.6 0.6]);
plot([mean(rho(WhichNeurons(R.Typ==0)==0)) mean(rho(WhichNeurons(R.Typ==0)==0))],[0 1],'color',[0.6 0.6 0.6],'linewidth',1.5);

legend([b{:}],'Out X Time','Others','location','southeast');

axis([-1 1 0 1]);
xlabel('Spearman''s \rho');
ylabel('Cumulative fraction of population');

%test compared to other neurons
int_pval = ranksum(rho(WhichNeurons(R.Typ==0)),rho(WhichNeurons(R.Typ==0)==0));
text(0.2,0.7,['p = ' num2str(round(int_pval,8))]);

%fraction negatively correlated
frac = sum(rho(WhichNeurons(R.Typ==0))<0 & pval(WhichNeurons(R.Typ==0))<0.05)/sum(WhichNeurons(R.Typ==0));
fracAll = sum(rho<0 & pval<0.05)/sum(R.Typ==0);

%% change in cues from first trial to last trial
NN=0;
for session = 1:length(RAW)
    for neuron = 1:length(R.Hz{session,2}(1,:))
        NN=NN+1;

        firstwattrial = find(R.TrialInfo{session,1}.Outcome==0 & R.TrialInfo{session,1}.Choice==0,1);
        lastwattrial = find(R.TrialInfo{session,1}.Outcome==0 & R.TrialInfo{session,1}.Choice==0,1,'last');
        firstWatTrialZ(NN,1) = (R.Hz{session,2}(firstwattrial,neuron) - R.Bmean(NN))/R.Bstd(NN);
        lastWatTrialZ(NN,1) = (R.Hz{session,2}(lastwattrial,neuron) - R.Bmean(NN))/R.Bstd(NN);
        
        firstsuctrial = find(R.TrialInfo{session,1}.Outcome==1 & R.TrialInfo{session,1}.Choice==0,1);
        lastsuctrial = find(R.TrialInfo{session,1}.Outcome==1 & R.TrialInfo{session,1}.Choice==0,1,'last');
        firstSucTrialZ(NN,1) = (R.Hz{session,2}(firstsuctrial,neuron) - R.Bmean(NN))/R.Bstd(NN);
        lastSucTrialZ(NN,1) = (R.Hz{session,2}(lastsuctrial,neuron) - R.Bmean(NN))/R.Bstd(NN);
    end
end

%cue change
figure;
WhichNeurons=R.GLM{1,1}.PVal.OutByPref(:,2)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,2)>0 & R.Typ==0;

% quadrant summary for RD
axes=[0.5 2.5 -2 4.5];
a=[];


%setting up parameters
Colors={[0.2 0.1 0],[0 0 0.1],[0 0.2 0];
        [0.4 0.1 0],[0.2 0.1 0.4],[0 0.4 0.3];
        [0.8 0.3 0],[0.3 0.2 0.7],[0.1 0.6 0.4];
        [1 0.5 0],[0.5 0.3 0.9],[0.2 0.8 0.6]};
titles={'Sucrose','Water'};
  
for i=1:2 %which reward



    rewmean=[];
    rewste=[];
    if i==2
        rewmean(i,1)=nanmean(firstWatTrialZ(WhichNeurons));
        rewste(i,1)=nanste(firstWatTrialZ(WhichNeurons),1);
        rewmean(i,2)=nanmean(lastWatTrialZ(WhichNeurons));
        rewste(i,2)=nanste(lastWatTrialZ(WhichNeurons),1);
    else
        rewmean(i,1)=nanmean(firstSucTrialZ(WhichNeurons));
        rewste(i,1)=nanste(firstSucTrialZ(WhichNeurons),1);
        rewmean(i,2)=nanmean(lastSucTrialZ(WhichNeurons));
        rewste(i,2)=nanste(lastSucTrialZ(WhichNeurons),1);
    end
    
    subplot(1,2,i)
    hold on

    rewup=rewmean(i,:)+rewste(i,:);
    rewdown=rewmean(i,:)-rewste(i,:);

    a{i}=plot(1:2,rewmean(i,:),'Color',Colors{3,i},'linewidth',1);
    patch([1:2,2:-1:1],[rewup,rewdown(end:-1:1)],Colors{3,i},'EdgeColor','none');alpha(0.5);
    axis(axes);
    xlabel('Trial');
    title(titles{i});
    ylabel({'Population mean (z-score)'});
    plot([0 5],[0 0],':','color','k','linewidth',0.75);
    xticks(1:2);
    xticklabels({'First';'Last'});
    
    waterchange = signrank(firstWatTrialZ(WhichNeurons),lastWatTrialZ(WhichNeurons));
    watervsuc = ranksum(firstWatTrialZ(WhichNeurons),firstSucTrialZ(WhichNeurons));

end


%% example neuron raster -- reward
figure;
WhichNeurons=R.GLM{1,1}.PVal.OutByPref(:,6)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,6)>0 & R.Typ;
NeuronList=find(WhichNeurons);

alph=0.3;
%raster
window=[0 3];
sample_neuron_number=367;
AnalysisWindow=[0 2.75];
AnalysisWindow=[0.75 1.95];

sessrat=R.Rat(sample_neuron_number);
sessneurons=find(R.Rat==sessrat&R.Typ);
sessneuron=sessneurons==sample_neuron_number;
session=find(unique(R.Rat)==sessrat)*2;

subplot(1,2,1);
hold on;


RD=strcmp('SucroseF', RAW(session).Einfo(:,2));
RDtimes=RAW(session).Erast{RD,1};

patch([AnalysisWindow(1) AnalysisWindow(1) AnalysisWindow(2) AnalysisWindow(2)],[0 length(RDtimes)*11 length(RDtimes)*11 0],[.8 .3 0],'edgecolor','none');alpha(alph);
PlotRaster(RAW(session).Nrast{sessneuron,1},RDtimes,window);
ylabel('Trials, sorted by session progress');
yticks([]);
xticks(-1:4);
xticklabels(-1:4);
xlabel('Seconds from reward delivery');
xlim(window);

subplot(1,2,2);
hold on;

RD=strcmp('WaterF', RAW(session).Einfo(:,2));
RDtimes=RAW(session).Erast{RD,1};

patch([AnalysisWindow(1) AnalysisWindow(1) AnalysisWindow(2) AnalysisWindow(2)],[0 length(RDtimes)*11 length(RDtimes)*11 0],[.3 .2 .7],'edgecolor','none');alpha(alph);
PlotRaster(RAW(session).Nrast{sessneuron,1},RDtimes,window);
ylabel('');
yticks([]);
xticks(-1:4);
xticklabels(-1:4);
xlim(window);

%% firing over time -- RD
WhichNeurons=R.GLM{1,1}.PVal.OutByPref(:,6)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,6)>0;
%WhichNeurons=WhichNeurons==0; %to look at other neurons

AnalysisWindow = [0.75 1.95];
alph=0.3;

for task=2
    figure;

    if task==1 Sel=WhichNeurons & R.Typ==0; end
    if task==2 Sel=WhichNeurons & R.Typ; end

    %setting up parameters
    Xaxis=[-1 5];
    inttime=find(R.Tm>=Xaxis(1) & R.Tm<=Xaxis(2)); 
    PSTHtime=R.Tm(inttime);

    Colors={[0.2 0.1 0],[0 0 0.1];
            [0.4 0.1 0],[0.2 0.1 0.4];
            [0.8 0.3 0],[0.3 0.2 0.7];
            [1 0.5 0],[0.5 0.3 0.9]};
    titles={'Sucrose signal, forced','Water signal, forced'};
    for i=1:2

        subplot(2,3,i)
        hold on

        p={};
        patch([AnalysisWindow(1) AnalysisWindow(1) AnalysisWindow(2) AnalysisWindow(2)],[-3 11 11 -3],[.8 .8 .8],'edgecolor','none');alpha(alph);
        for j=1:4
            
            %plot each quadrant
            psth1=nanmean(R.Ev((i-1)*4+j).PSTHz(Sel,inttime),1); 
            sem1=nanste(R.Ev((i-1)*4+j).PSTHz(Sel,inttime),1); %calculate standard error of the mean
            upE=psth1+sem1;
            downE=psth1-sem1;


            %plotting
            p{j}=plot(PSTHtime,psth1,'Color',Colors{j,i},'linewidth',1);
            patch([PSTHtime,PSTHtime(end:-1:1)],[upE,downE(end:-1:1)],Colors{j,i},'EdgeColor','none');alpha(0.5);
        end

        plot([-4 5],[0 0],':','color','k','linewidth',0.75);
        plot([0 0],[-3 12],':','color','k','linewidth',0.75);
        axis([Xaxis(1) Xaxis(2) -3 5]);
        ylabel('Mean firing (z-score)');
        xlabel('Seconds from reward delivery');
        legend([p{:}],'Q1','Q2','Q3','Q4');
        title(titles{i})


    end


end

% quadrant summary for RD
axes=[0 5 -2 4];
a=[];
for task=2

    if task==1 Sel=WhichNeurons & R.Typ==0; end
    if task==2 Sel=WhichNeurons & R.Typ; end

    %setting up parameters
    Colors={[0.2 0.1 0],[0 0 0.1],[0 0.2 0];
            [0.4 0.1 0],[0.2 0.1 0.4],[0 0.4 0.3];
            [0.8 0.3 0],[0.3 0.2 0.7],[0.1 0.6 0.4];
            [1 0.5 0],[0.5 0.3 0.9],[0.2 0.8 0.6]};
    titles={'Binned firing','Action','Reward'};
    
    %time periods we're plotting   
    rewbin=[0.75 1.95];
    rewtoi=find(R.Tm>=rewbin(1) & R.Tm<=rewbin(2));  
    
    for i=1:2 %which reward



        rewmean=[];
        rewste=[];
        for j=1:4 %which quadrant
            %get mean for each quadrant
            rewmean(i,j)=nanmean(nanmean(R.Ev((i-1)*4+j).PSTHnosmoothz(Sel,rewtoi),2));
            rewste(i,j)=nanste(nanmean(R.Ev((i-1)*4+j).PSTHnosmoothz(Sel,rewtoi),2),1);
        end
        
        subplot(2,3,3)
        hold on
        
        rewup=rewmean(i,:)+rewste(i,:);
        rewdown=rewmean(i,:)-rewste(i,:);
        
        a{i}=plot(1:4,rewmean(i,:),'Color',Colors{3,i},'linewidth',1);
        patch([1:4,4:-1:1],[rewup,rewdown(end:-1:1)],Colors{3,i},'EdgeColor','none');alpha(0.5);
        axis(axes);
        if task==2 xlabel('Quadrant'); end
        if task==2 title(titles{1}); end
        ylabel({'Population mean (z-score)'});
        plot([0 5],[0 0],':','color','k','linewidth',0.75);
        xticks(1:4);
        if task==2 legend([a{:}],'Sucrose','Water'); end

        
    end

    

end

%% correlation with preference
figure;
%example neuron
sample_neuron_number=367; %215, 225, 367
AnalysisWindow=[0 2.75];
AnalysisWindow=[0.75 1.95];

Colors={[0.2 0.1 0],[0 0 0.1],[0 0.2 0];
        [0.4 0.1 0],[0.2 0.1 0.4],[0 0.4 0.3];
        [0.8 0.3 0],[0.3 0.2 0.7],[0.1 0.6 0.4];
        [1 0.5 0],[0.5 0.3 0.9],[0.2 0.8 0.6]};

sessrat=R.Rat(sample_neuron_number);
sessneurons=find(R.Rat==sessrat&R.Typ);
sessneuron=sessneurons==sample_neuron_number;
session=find(unique(R.Rat)==sessrat)*2;

preference = R.TrialInfo{session,1}.Preference;
sucfiring = R.Hz{session,6}(R.TrialInfo{session,1}.Outcome==1 & R.TrialInfo{session,1}.Choice==0,sessneuron);
watfiring = R.Hz{session,6}(R.TrialInfo{session,1}.Outcome==0 & R.TrialInfo{session,1}.Choice==0,sessneuron);

subplot(2,2,1);
hold on;
plot(1-(preference-min(preference))/(max(preference)-min(preference)),'Color',Colors{3,3},'linewidth',1);
scatter(find(R.TrialInfo{session,1}.Outcome==0 & R.TrialInfo{session,1}.Choice==0),(watfiring-min(watfiring))/(max(watfiring)-min(watfiring)),12,Colors{3,2});
xlabel('Trial');
ylabel({'Water preference';'Water firing (normalized)'});
title('Example neuron');

subplot(2,2,2);
hold on;
plot((preference-min(preference))/(max(preference)-min(preference)),'Color',Colors{3,3},'linewidth',1);
scatter(find(R.TrialInfo{session,1}.Outcome==1 & R.TrialInfo{session,1}.Choice==0),(sucfiring-min(sucfiring))/(max(sucfiring)-min(sucfiring)),12,Colors{3,1});
ylabel({'Sucrose preference';'Sucrose firing (normalized)'});

watrho=[];
watpval=[];
sucrho=[];
sucpval=[];
%run correlations
NN = 0;
for session = 1:length(RAW)
    wattrials = R.TrialInfo{session,1}.Outcome==0 & R.TrialInfo{session,1}.Choice==0;
    suctrials = R.TrialInfo{session,1}.Outcome==1 & R.TrialInfo{session,1}.Choice==0;
    for neuron = 1:length(R.Hz{session,6}(1,:))
        NN = NN+1;
        [watrho(NN,1),watpval(NN,1)] = corr(R.Hz{session,6}(wattrials,neuron),1-R.TrialInfo{session,1}.Preference(wattrials),'type','spearman');
        [sucrho(NN,1),sucpval(NN,1)] = corr(R.Hz{session,6}(suctrials,neuron),R.TrialInfo{session,1}.Preference(suctrials),'type','spearman');
    end
end

%summary of correlation coefficients
WhichNeurons=R.GLM{1,1}.PVal.OutByPref(:,6)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,6)>0;
b={};

colors={};
colors{2} = [0.6 0.6 0.6];
titles={'Water correlation','Sucrose correlation'};
for reward=1:2
    colors{1} = Colors{3,-reward+3};
    
    subplot(2,2,2+reward);
    hold on;

    plot([-1 1],[0.5 0.5],'color','k','linewidth',1);
    plot([0 0],[0 01],'color','k','linewidth',1);

    if reward==1 rho=watrho; end
    if reward==2 rho=sucrho; end
    for condition=1:2
        values = rho(WhichNeurons==(-condition+2)&R.Typ==1);
        [cdf,x] = ecdf(values);
        b{condition} = plot(x,cdf,'linewidth',1.5,'color',colors{condition});
        plot([mean(values) mean(values)],[0 1],'color',colors{condition},'linewidth',1.5);
    end

    legend([b{:}],'Out X Time','Others','location','northwest');
    axis([-1 1 0 1]);
    xlabel('Spearman''s \rho');
    ylabel('Cumulative fraction of population');
    
    int_pval = ranksum(rho(WhichNeurons==(-condition+2)==0&R.Typ==1),rho(WhichNeurons==(-condition+2)&R.Typ==1));
    text(-0.9,0.45,['p = ' num2str(round(int_pval,2,'significant'))]);
    title(titles{reward});
end
