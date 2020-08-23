load('R_choice.mat')
load('RAWCh.mat')

%for selection of best model
load('ModelFits_CV.mat')
mod_RD=select_AllModeledModel(os, 'RD', 'particularModels', {'mean','tv','pref','spp'},'normalizeHistogram_Flag',true,'scoreToUse', 'CVLH');

%for using fit parameters
load('ModelFits.mat')
%% example neuron fit
figure;

%which neuron to plot?
neuron=130;

Colors={[0.2 0.1 0],[0 0 0.1],[0 0.2 0];
        [0.4 0.1 0],[0.2 0.1 0.4],[0 0.4 0.3];
        [0.8 0.3 0],[0.3 0.2 0.7],[0.1 0.6 0.4];
        [1 0.5 0],[0.5 0.3 0.9],[0.2 0.8 0.6]};

%absolute value
models={'tv','pref','spp'};
titles={'Satiety','Preference','Mixed'};
for curr_mod=1:3
    
    subplot(1,4,1+curr_mod);
    hold on;

    model=models{curr_mod};
    measure='RPEs'; %'RPEs' or 'V'

    x=find(os(neuron).Outcome==0 & os(neuron).Included)/find(os(neuron).Outcome==0 & os(neuron).Included,1,'last');
    y=os(neuron).mod_RD.(model).(measure)(os(neuron).Outcome==0 & os(neuron).Included);

    plot(x,y,'color',Colors{3,2},'linewidth',2);

    x=find(os(neuron).Outcome==1 & os(neuron).Included)/find(os(neuron).Outcome==1 & os(neuron).Included,1,'last');
    y=os(neuron).mod_RD.(model).(measure)(os(neuron).Outcome==1 & os(neuron).Included);

    plot(x,y,'color',Colors{3,1},'linewidth',2);
    
    if curr_mod==1 legend('Water','Sucrose'); end
    title(titles{curr_mod});
    xlabel('Session progress');
    if curr_mod==1 ylabel('Reward value'); end
    axis([-.05 1.05 -0.05 1.05]);
    
end

subplot(1,4,1);
hold on;

x=find(os(neuron).Outcome==0 & os(neuron).Included)/find(os(neuron).Outcome==0 & os(neuron).Included,1,'last');
y=os(neuron).spikeCount_RD(os(neuron).Outcome==0 & os(neuron).Included);

plot(x,y,'color',Colors{3,2},'linewidth',2);

x=find(os(neuron).Outcome==1 & os(neuron).Included)/find(os(neuron).Outcome==1 & os(neuron).Included,1,'last');
y=os(neuron).spikeCount_RD(os(neuron).Outcome==1 & os(neuron).Included);

plot(x,y,'color',Colors{3,1},'linewidth',2);

legend('Water','Sucrose','location','southwest');
ylabel('Spikes/trial');
xlabel('Session progress');
axis([-.05 1.05 0 17]);
%% classification summary
WhichNeurons=R.GLM{1,1}.PVal.OutByPref(:,6)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,6)>0;

figure;

masks=[];
masks(:,4) = mod_RD.mask_spp';
masks(:,3) = mod_RD.mask_pref';
masks(:,2) = mod_RD.mask_tv';
masks(:,1) = mod_RD.mask_mean';

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
title('Best fit model (cross-validated likelihood)');


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
    overallsession = ratind * 2 - 1 + task;

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
            
%             if (mod_RD.('mask_pref')(i)) model='pref'; end
%             if (mod_RD.('mask_spp')(i)) model='spp'; end
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
    %histogram(mod_midpoint{ratind,2},0:5:200,'normalization','probability');
    
%     [cdf,x] = ecdf(mod_midpoint{ratind,1});
%     plot(x,cdf,'linewidth',1); %'color',[0 0.6 0.8]
%     [cdf,x] = ecdf(mod_midpoint{ratind,2});
%     plot(x,cdf,'linewidth',1); %'color',[0 0.6 0.8]
    
    height=1;
    plot([fit_midpoint fit_midpoint],[0 height],'color','r','linewidth',1);
    patch([fit_midpoint-fit_midpoint_SE fit_midpoint-fit_midpoint_SE fit_midpoint+fit_midpoint_SE fit_midpoint+fit_midpoint_SE],[0 height height 0],'r','EdgeColor','none','FaceAlpha',0.5);
    axis([-10 180 0 round(max(h.Values)+0.05,1)]);
    yticks(0:0.2:round(max(h.Values)+0.05,1));
    
    if ratind==1
        title('Indifference point');
        legend('Neural estimates','Behavioral est.');
    end
    
    if ratind==5
        xlabel('Trials');
    end
    
    subplot(length(rats),4,3+4*(ratind-1));
    hold on;
    h1=histogram(preference_corr{ratind,1},-1:0.1:1,'normalization','probability','edgecolor','none','facecolor',Colors{3,3},'FaceAlpha',1);
%    h2=histogram(preference_corr{ratind,2},-1:0.1:1,'normalization','probability','edgecolor','none','facecolor',[0.6 0.6 0.6],'FaceAlpha',0.6);
    
%     [cdf,x] = ecdf(mod_midpoint{ratind,1});
%     plot(x,cdf,'linewidth',1); %'color',[0 0.6 0.8]
%     [cdf,x] = ecdf(mod_midpoint{ratind,2});
%     plot(x,cdf,'linewidth',1); %'color',[0 0.6 0.8]
    
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

%% classification summary with licks
load('ModelFits_CV.mat')
mod_RD=select_AllModeledModel(os, 'RD', 'particularModels', {'mean','tv','pref','spp','licks'},'normalizeHistogram_Flag',true,'scoreToUse', 'CVLH');

load('ModelFits.mat')
WhichNeurons=R.GLM{1,1}.PVal.OutByPref(:,6)<0.05 & R.GLM{1,1}.Coefficients.OutByPref(:,6)>0;

figure;

masks=[];
masks(:,5) = mod_RD.mask_licks';
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
axis([0 6 0 round(max(max(stacked_data))+5,-1)])
xticks(1:5);
xticklabels({'Unmodulated';'Satiety';'Preference';'Mixed';'Licks'});
xtickangle(45);
yticks([0:40:round(max(max(stacked_data))+5,-1)])
ylabel('Number of neurons')
title('Best fit model (cross-validated likelihood)');



