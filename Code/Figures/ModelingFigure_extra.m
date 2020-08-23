load('ModelFits_extra.mat')
load('R_choice_extra.mat')
load('RAWCh_extra.mat')

mod_RD=select_AllModeledModel(os, 'RD', 'particularModels', {'mean','tv','pref','spp'},'normalizeHistogram_Flag',true,'scoreToUse', 'AIC');
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
for task=1

    if task==1 Sel=cat(1,RAW.Ses)>10; end


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

    plot(1:4,Preference(Sel,:)','Color',Colors{3,3},'linewidth',1);
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


%%  cue and reward proportions
predictable=[0.9 0.9 0.9];
unpredictable=[0.1 0.1 0.1];
figure;

props=[];
for trial=1:2 %forced or choice

    
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


%% classification summary
WhichNeurons=R.GLM{1,1}.PVal.OutByPref(:,6)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,6)>0;

figure;

masks=[];
masks(:,4) = mod_RD.mask_spp';
masks(:,3) = mod_RD.mask_pref';
masks(:,2) = mod_RD.mask_tv';
masks(:,1) = mod_RD.mask_mean';
task=0;

stacked_data = zeros(2,6);

for task=1
    for subset=1:length(masks(1,:))
        stacked_data(task+1,subset)=sum(masks(:,subset)&R.Typ==task&WhichNeurons);
    end
end

subplot(1,1,1)
hold on;
b = bar(stacked_data(task+1,:),'edgecolor','none','facecolor','flat');
axis([0 5 0 round(max(max(stacked_data))+5,-1)])
xticks(1:4);
xticklabels({'Unmodulated';'Satiety';'Preference';'Mixed'});
xtickangle(45);
yticks([0:40:round(max(max(stacked_data))+5,-1)])
ylabel('Number of neurons')
title('Best fit model (AIC)');

%% traces from model fits per rat
WhichNeurons=R.GLM{1,1}.PVal.OutByPref(:,6)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,6)>0;


task = 1;

figure;
Colors={[0.2 0.1 0],[0 0 0.1],[0 0.2 0];
        [0.4 0.1 0],[0.2 0.1 0.4],[0 0.4 0.3];
        [0.8 0.3 0],[0.3 0.2 0.7],[0.1 0.6 0.4];
        [1 0.5 0],[0.5 0.3 0.9],[0.2 0.8 0.6]};
       
        
mod_preference={};
mod_midpoint={};
preference_corr={};
distance={};

rats = unique(R.Rat);
for ratind=1:length(rats)
    rat=rats(ratind);
    overallsession = ratind;

    %fit sigmoid to choice data
    choicetrials = find(R.TrialInfo{overallsession,1}.Choice==1);
    choices = R.TrialInfo{overallsession,1}.Outcome(R.TrialInfo{overallsession,1}.Choice==1);
    
    x=choicetrials;
    y=choices;
    
    f = @(F,x) 1 ./ (1 + exp(-F(2).*(x - F(1))));
    [F_fitted,~,~,CovB] = nlinfit(x,y,f,[40 .1]);
    fit_midpoint = F_fitted(1);
    
    SE = sqrt(diag(CovB));
    fit_midpoint_SE = SE(1);
    
    all_x = 1:length(R.TrialInfo{overallsession,1}.Outcome);
    fit_preference = 1 ./ (1 + exp(-F_fitted(2).*(all_x - F_fitted(1))));
    
    subplot(length(rats),4,1+4*(ratind-1));
    hold on;    

    %plot choices
    tickwidth=1;

    pl={};
    choices=find(R.TrialInfo{overallsession,1}.Choice==1&R.TrialInfo{overallsession,1}.Outcome==1);
    for i=1:length(choices)
        pl{1}=plot([choices(i) choices(i)],[0.8 1],'color',Colors{3,1},'linewidth',tickwidth);
    end

    choices=find(R.TrialInfo{overallsession,1}.Choice==1&R.TrialInfo{overallsession,1}.Outcome==0);
    for i=1:length(choices)
        pl{2}=plot([choices(i) choices(i)],[0 0.2],'color',Colors{3,2},'linewidth',tickwidth);
    end
    
    plot(find(R.TrialInfo{overallsession,1}.Outcome==1),fit_preference(R.TrialInfo{overallsession,1}.Outcome==1),'Color',Colors{3,1},'linewidth',1);
    plot(find(R.TrialInfo{overallsession,1}.Outcome==0),1-fit_preference(R.TrialInfo{overallsession,1}.Outcome==0),'Color',Colors{3,2},'linewidth',1);

    adjusted_preference = fit_preference;
    adjusted_preference(R.TrialInfo{overallsession,1}.Outcome==0) = 1-fit_preference(R.TrialInfo{overallsession,1}.Outcome==0);
    
    if ratind==1
        title('Behavioral preference');
        legend([pl{:}],'sucrose','water');
        ylabel('Preference');
    end
    axis([0 length(R.TrialInfo{overallsession,1}.Outcome) 0 1]);
    
    NN=0;
    NC=0;
    for i=1:length(os)
        if WhichNeurons(i) & R.Rat(i)==rat & R.Typ(i)==task
            NN=NN+1;
            
            model='spp';
            mod_preference{ratind,1}(NN,:)=os(i).mod_RD.(model).V;
            mod_midpoint{ratind,1}(NN,:)=os(i).mod_RD.(model).bestParams(3) * length(os(i).Outcome);
            preference_corr{ratind,1}(NN,1)=corr(os(i).mod_RD.(model).V,adjusted_preference');
        elseif R.Rat(i)==rat & R.Typ(i)==task %all non-pref neurons
            NC=NC+1;
            model='spp';
            mod_midpoint{ratind,2}(NC,:)=os(i).mod_RD.(model).bestParams(3) * length(os(i).Outcome);
            preference_corr{ratind,2}(NC,1)=corr(os(i).mod_RD.(model).V,adjusted_preference');
        end
    end
    
    ratprefmean = mean(mod_preference{ratind,1},1);
    ratprefste = nanste(mod_preference{ratind,1},1);
    ratprefup = ratprefmean + ratprefste;
    ratprefdown = ratprefmean - ratprefste;
    
    subplot(length(rats),4,2+4*(ratind-1));
    hold on;
    for reward=1:2
        Sel = R.TrialInfo{overallsession,1}.Outcome==(2-reward);
        trials = find(Sel)';
        ratprefup = ratprefmean(Sel) + ratprefste(Sel);
        ratprefdown = ratprefmean(Sel) - ratprefste(Sel);        
        plot(find(Sel),ratprefmean(Sel),'Color',Colors{3,reward},'linewidth',1);
        patch([trials,trials(end:-1:1)],[ratprefup,ratprefdown(end:-1:1)],Colors{3,reward},'EdgeColor','none');alpha(0.5);
    end
    text(0.1,0.5,['n=' num2str(NN)]);
    axis([0 length(R.TrialInfo{overallsession,1}.Outcome) 0 1]);
        
    if ratind==1
        title('Neural estimate');
    end
    
    subplot(length(rats),4,4+4*(ratind-1));
    hold on;
    h=histogram(mod_midpoint{ratind,1},-10:10:200,'normalization','probability','edgecolor','none','facecolor',Colors{3,3},'FaceAlpha',1);

    
    height=1;
    plot([fit_midpoint fit_midpoint],[0 height],'color','r','linewidth',1);
    patch([fit_midpoint-fit_midpoint_SE fit_midpoint-fit_midpoint_SE fit_midpoint+fit_midpoint_SE fit_midpoint+fit_midpoint_SE],[0 height height 0],'r','EdgeColor','none','FaceAlpha',0.5);
    axis([-10 180 0 round(max(h.Values)+0.05,1)]);
    yticks(0:0.2:round(max(h.Values)+0.05,1));
    
    if ratind==1
        title('Indifference point');
        legend('Neural estimates','Behavioral est.');
    end
    
    if ratind==4
        xlabel('Trials');
    end
    
    subplot(length(rats),4,3+4*(ratind-1));
    hold on;
    h1=histogram(preference_corr{ratind,1},-1:0.1:1,'normalization','probability','edgecolor','none','facecolor',Colors{3,3},'FaceAlpha',1);

    
    axis([-1 1 0 round(max([h1.Values])+0.05,1)]);
    yticks(0:0.2:round(max([h1.Values])+0.05,1));

    if ratind==1
        title('Correlation');
        ylabel('Fraction');
    end
    
    if ratind==5
        xlabel('Coefficient');
    end
    
    distance{ratind,1}=abs(fit_midpoint - mod_midpoint{ratind,1});
    distance{ratind,2}=abs(fit_midpoint - mod_midpoint{ratind,2});
    
    behav_midpoints(ratind,1) = fit_midpoint;
    
end

%statistical tests for modeling results
figure; 

%correlation for preference neurons and others
corr_pref=cat(1,preference_corr{:,1});
corr_cont=cat(1,preference_corr{:,2});

subplot(1,2,1);
hold on;
[cdf,x] = ecdf(corr_pref);
plot(x,cdf,'linewidth',1.5,'color',Colors{3,3});
[cdf,x] = ecdf(corr_cont);
plot(x,cdf,'linewidth',1.5,'color',[0.6 0.6 0.6]);

plot([mean(corr_pref) mean(corr_pref)],[0 1],'color',Colors{3,3},'linewidth',1.5);
plot([mean(corr_cont) mean(corr_cont)],[0 1],'color',[0.6 0.6 0.6],'linewidth',1.5);

p = ranksum(corr_pref,corr_cont);
xlabel('Coefficient');
ylabel('Fraction of neural population');
text(-0.9,0.7,['p = ' num2str(round(p,2,'significant'))]);
axis([-1 1 0 1]);
plot([-1 1],[0.5 0.5],'color','k','linewidth',1);
plot([0 0],[-1 1],'color','k','linewidth',1);

title('Correlation with behavioral preference');
legend('Out X Time','Others','location','northeast');

%midpoint estimate for preference neurons and others
distance_pref=cat(1,distance{:,1});
distance_cont=cat(1,distance{:,2});

subplot(1,2,2);
hold on;
[cdf,x] = ecdf(distance_pref);
plot(x,cdf,'linewidth',1.5,'color',Colors{3,3});
[cdf,x] = ecdf(distance_cont);
plot(x,cdf,'linewidth',1.5,'color',[0.6 0.6 0.6]);


plot([mean(distance_pref) mean(distance_pref)],[0 1],'color',Colors{3,3},'linewidth',1.5);
plot([mean(distance_cont) mean(distance_cont)],[0 1],'color',[0.6 0.6 0.6],'linewidth',1.5);

p = ranksum(distance_pref,distance_cont);
xlabel('Error (no. trials)');
text(100,0.3,['p = ' num2str(round(p,2,'significant'))]);
axis([0 150 0 1]);
title('Accuracy of indifference point estimate');
legend('Out X Time','Others','location','southeast');


