close all; clear all;

tic;
%Parameter setting
z1= 1; % Upper decision threshold (choice 1)
z2= -1; % Lower decision threshol (choice 2)
dt=1; % Time step 5ms

trials=2500; % Trial number

k=9.66*10^(-3); % the proportionality factor to form the mean of the drift rate
Mo_Strength=[0:0.01:0.03 0.032 0.04:0.01:0.25 0.256 0.26:0.01:1]; % the the motion strength
ti=1:dt:2000; % Time < 2s
Sy=7.34; % when G both: 7.34, sigle 0.734
Sx=0.0028;
d=584; % by ms
G=Sy.*exp(Sx*(ti-d))./(1+exp(Sx*(ti-d)))+(1+(1-Sy)*exp(-Sx*d))/(1+exp(-Sx*d));
sigma0=0.0188; % Size of the noise 

AveDTc_initial=zeros(1,length(Mo_Strength));
AveDTe_initial=zeros(1,length(Mo_Strength));
ac1=zeros(1,length(Mo_Strength));

AveDTc_both=zeros(1,length(Mo_Strength));
AveDTe_both=zeros(1,length(Mo_Strength));
ac2=zeros(1,length(Mo_Strength));

AveDTc_onlyNoise=zeros(1,length(Mo_Strength));
AveDTe_onlyNoise=zeros(1,length(Mo_Strength));
ac3=zeros(1,length(Mo_Strength));

AveDTc_onlyDrift=zeros(1,length(Mo_Strength));
AveDTe_onlyDrift=zeros(1,length(Mo_Strength));
ac4=zeros(1,length(Mo_Strength));

%% time-variant gain on both drift rate and noise term with noise depending on c
j=1;
for c=Mo_Strength
    DTc1=zeros(1,trials); % Correct decision time
    DTe1=zeros(1,trials); % Error decision time 
    Mu0=k*c;% Drift rate
    sigma=sigma0*sqrt(1+c); 
    for i=1:trials % Trial number
        x1=zeros(1,length(ti)); % initial x value = '0'
        for t=1:length(ti)
            x1(t+1)=x1(t) + dt*Mu0*G(t) + sqrt(dt)*sigma*G(t)*randn;% Update x;
            if x1(t) >= 1 % Record the correct decision time
                DTc1(i)=t;
                break;
            end
            if x1(t) <= -1 % Record the error decision time
                DTe1(i)=t;
                break;
            end
        end
    end
    DTc1(DTc1==0)=[];
    DTe1(DTe1==0)=[];
    if(c==0.032)
        ma=max([DTc1 DTe1]);
        mi=min([DTc1 DTe1]);
        DTc1_032=(DTc1-mi)/(ma-mi);
        DTe1_032=(DTe1-mi)/(ma-mi);
        DTc1_032(DTc1_032==0)=[];
        DTe1_032(DTe1_032==0)=[];
    end
    if(c==0.256)
        ma=max([DTc1 DTe1]);
        mi=min([DTc1 DTe1]);
        DTc1_256=(DTc1-mi)/(ma-mi);
        DTe1_256=(DTe1-mi)/(ma-mi);
        DTc1_256(DTc1_256==0)=[];
        DTe1_256(DTe1_256==0)=[];      
    end
    AveDTc_initial(j)=dt.*mean(DTc1); %Calculate the average decision time
    AveDTe_initial(j)=dt.*mean(DTe1);
    ac1(j)=length(DTc1)/(length(DTc1)+length(DTe1));
    j=j+1;
end
ma=max([AveDTc_initial AveDTe_initial]);
mi=min([AveDTc_initial AveDTe_initial]);
AveDTc_initial=(AveDTc_initial-mi)/(ma-mi);
AveDTe_initial=(AveDTe_initial-mi)/(ma-mi);

%% time-variant gain on both drift rate and noise term with noise independent of c
j=1;
for c=Mo_Strength
    DTc2=zeros(1,trials); % Correct decision time
    DTe2=zeros(1,trials); % Error decision time 
    Mu0=k*c;% Drift rate
    for i=1:trials % Trial number
        x2=zeros(1,length(ti)); % initial x value = '0'
        for t=1:length(ti)
            x2(t+1)=x2(t) + dt*Mu0*G(t) + sqrt(dt)*sigma0*G(t)*randn;% Update x;
            if x2(t) >= 1 % Record the correct decision time
                DTc2(i)=t;
                break;
            end
            if x2(t) <= -1 % Record the error decision time
                DTe2(i)=t;
                break;
            end
        end
    end
    DTc2(DTc2==0)=[];
    DTe2(DTe2==0)=[];
    if(c==0.032)
        ma=max([DTc2 DTe2]);
        mi=min([DTc2 DTe2]);
        DTc2_032=(DTc2-mi)/(ma-mi);
        DTe2_032=(DTe2-mi)/(ma-mi);
        DTc2_032(DTc2_032==0)=[];
        DTe2_032(DTe2_032==0)=[]; 
    end
    if(c==0.256)
        ma=max([DTc2 DTe2]);
        mi=min([DTc2 DTe2]);
        DTc2_256=(DTc2-mi)/(ma-mi);
        DTe2_256=(DTe2-mi)/(ma-mi);
        DTc2_256(DTc2_256==0)=[];
        DTe2_256(DTe2_256==0)=[]; 
    end
    AveDTc_both(j)=dt.*mean(DTc2); %Calculate the average decision time
    AveDTe_both(j)=dt.*mean(DTe2);
    ac2(j)=length(DTc2)/(length(DTc2)+length(DTe2));
    j=j+1;
end
ma=max([AveDTc_both AveDTe_both]);
mi=min([AveDTc_both AveDTe_both]);
AveDTc_both=(AveDTc_both-mi)/(ma-mi);
AveDTe_both=(AveDTe_both-mi)/(ma-mi);

%% time-variant gain only on noise term

j=1;
for c=Mo_Strength
    DTc3=zeros(1,trials); % Correct decision time
    DTe3=zeros(1,trials); % Error decision time 
    Mu0=k*c;% Drift rate
    for i=1:trials % Trial number
        x3=zeros(1,length(ti)); % initial x value = '0'
        for t=1:length(ti)
            x3(t+1)=x3(t) + dt*Mu0 + sqrt(dt)*G(t)*sigma0*randn;% Update x;
            if x3(t) >= 1 % Record the correct decision time
                DTc3(i)=t;
                break;
            end
            if x3(t) <= -1 % Record the error decision time
                DTe3(i)=t;
                break;
            end
        end
    end
    DTc3(DTc3==0)=[];
    DTe3(DTe3==0)=[];
    if(c==0.032)
        ma=max([DTc3 DTe3]);
        mi=min([DTc3 DTe3]);
        DTc3_032=(DTc3-mi)/(ma-mi);
        DTe3_032=(DTe3-mi)/(ma-mi);
        DTc3_032(DTc3_032==0)=[];
        DTe3_032(DTe3_032==0)=[]; 
    end
    if(c==0.256)
        ma=max([DTc3 DTe3]);
        mi=min([DTc3 DTe3]);
        DTc3_256=(DTc3-mi)/(ma-mi);
        DTe3_256=(DTe3-mi)/(ma-mi);
        DTc3_256(DTc3_256==0)=[];
        DTe3_256(DTe3_256==0)=[]; 
    end
    AveDTc_onlyNoise(j)=dt.*mean(DTc3); %Calculate the average decision time
    AveDTe_onlyNoise(j)=dt.*mean(DTe3);
    ac3(j)=length(DTc3)/(length(DTc3)+length(DTe3));
    j=j+1;
end
ma=max([AveDTc_onlyNoise AveDTe_onlyNoise]);
mi=min([AveDTc_onlyNoise AveDTe_onlyNoise]);
AveDTc_onlyNoise=(AveDTc_onlyNoise-mi)/(ma-mi);
AveDTe_onlyNoise=(AveDTe_onlyNoise-mi)/(ma-mi);

%% time-variant gain only on drift rate
j=1;
for c=Mo_Strength
    DTc4=zeros(1,trials); % Correct decision time
    DTe4=zeros(1,trials); % Error decision time 
    Mu0=k*c;% Drift rate
    for i=1:trials % Trial number
        x4=zeros(1,length(ti)); % initial x value = '0'
        for t=1:length(ti)
            x4(t+1)=x4(t) + dt*Mu0*G(t) + sqrt(dt)*sigma0*randn;% Update x;
            if x4(t) >= 1 % Record the correct decision time
                DTc4(i)=t;
                break;
            end
            if x4(t) <= -1 % Record the error decision time
                DTe4(i)=t;
                break;
            end
        end
    end
    DTc4(DTc4==0)=[];
    DTe4(DTe4==0)=[];
    if(c==0.032)
        ma=max([DTc4 DTe4]);
        mi=min([DTc4 DTe4]);
        DTc4_032=(DTc4-mi)/(ma-mi);
        DTe4_032=(DTe4-mi)/(ma-mi);
        DTc4_032(DTc4_032==0)=[];
        DTe4_032(DTe4_032==0)=[];
    end
    if(c==0.256)
        ma=max([DTc4 DTe4]);
        mi=min([DTc4 DTe4]);
        DTc4_256=(DTc4-mi)/(ma-mi);
        DTe4_256=(DTe4-mi)/(ma-mi);
        DTc4_256(DTc4_256==0)=[];
        DTe4_256(DTe4_256==0)=[];
    end
    AveDTc_onlyDrift(j)=dt.*mean(DTc4); %Calculate the average decision time
    AveDTe_onlyDrift(j)=dt.*mean(DTe4);
    ac4(j)=length(DTc4)/(length(DTc4)+length(DTe4));
    j=j+1;
end
ma=max([AveDTc_onlyDrift AveDTe_onlyDrift]);
mi=min([AveDTc_onlyDrift AveDTe_onlyDrift]);
AveDTc_onlyDrift=(AveDTc_onlyDrift-mi)/(ma-mi);
AveDTe_onlyDrift=(AveDTe_onlyDrift-mi)/(ma-mi);

%% Plot the 4x2 MotionStrength vs DT and Accuracy figure
X1=100*Mo_Strength;
YMatrix1=[AveDTc_initial; AveDTe_initial];
YMatrix2=[AveDTc_both; AveDTe_both];
YMatrix3=[AveDTc_onlyNoise; AveDTe_onlyNoise];
YMatrix4=[AveDTc_onlyDrift; AveDTe_onlyDrift];

figure1 = figure('OuterPosition',[31.4 40.2 1288.8 795.2]);

% create subplot
subplot1 = subplot(4,2,1);
hold(subplot1,'on');
plot1 = plot(X1,YMatrix1,'LineWidth',5);
set(plot1(1),'DisplayName','Correct','Color',[0 0 1]);
set(plot1(2),'DisplayName','Error','Color',[1 0 0]);
title('Case (i)');
% Set the rest of the axes properties
set(subplot1,'FontSize',16,'LineWidth',2,'TickDir','out');
% create legend
legend1 = legend(subplot1,'show');
set(legend1,...
    'Position',[0.364390839087796 0.865014403540914 0.0896107984202805 0.0765746516829394]);

% create subplot
subplot2 = subplot(4,2,2);
hold(subplot2,'on');
plot(X1,ac1,'LineWidth',3.5,'Color',[0.717647058823529 0.274509803921569 1]);
% Set the rest of the axes properties
set(subplot2,'FontSize',16,'LineWidth',2,'TickDir','out','YTick',...
    [0.5 0.75 1],'YTickLabel',{'50','75','100'});

% create subplot
subplot3 = subplot(4,2,3);
hold(subplot3,'on');
plot2 = plot(X1,YMatrix2,'LineWidth',5);
set(plot2(1),'Color',[0 0 1]);
set(plot2(2),'Color',[1 0 0]);
% create title
title('Case (ii)');
% Set the rest of the axes properties
set(subplot3,'FontSize',16,'LineWidth',2,'TickDir','out');

% create subplot
subplot4 = subplot(4,2,4);
hold(subplot4,'on');
plot(X1,ac2,'LineWidth',3.5,'Color',[0.717647058823529 0.274509803921569 1]);
% Set the rest of the axes properties
set(subplot4,'FontSize',16,'LineWidth',2,'TickDir','out','YTick',...
    [0.5 0.75 1],'YTickLabel',{'50','75','100'});

% create subplot
subplot5 = subplot(4,2,5);
hold(subplot5,'on');
plot3 = plot(X1,YMatrix3,'LineWidth',5);
set(plot3(1),'Color',[0 0 1]);
set(plot3(2),'Color',[1 0 0]);
% create ylabel
ylabel('Normalized mean reaction time');
% create title
title('Case (iii)');
% Set the rest of the axes properties
set(subplot5,'FontSize',16,'LineWidth',2,'TickDir','out');

% create subplot
subplot6 = subplot(4,2,6);
hold(subplot6,'on');
plot(X1,ac3,'LineWidth',3.5,'Color',[0.717647058823529 0.274509803921569 1]);
% create ylabel
ylabel('% correct');
% Set the rest of the axes properties
set(subplot6,'FontSize',16,'LineWidth',2,'TickDir','out','YTick',...
    [0.5 0.75 1],'YTickLabel',{'50','75','100'});

% create subplot
subplot7 = subplot(4,2,7);
hold(subplot7,'on');
plot4 = plot(X1,YMatrix4,'LineWidth',5);
set(plot4(1),'Color',[0 0 1]);
set(plot4(2),'Color',[1 0 0]);
% create xlabel
xlabel('Motion strength (% coh)');
% create title
title('Case (iv)');
% Set the rest of the axes properties
set(subplot7,'FontSize',16,'LineWidth',2,'TickDir','out');

% create subplot
subplot8 = subplot(4,2,8);
hold(subplot8,'on');
plot(X1,ac4,'LineWidth',3.5,'Color',[0.717647058823529 0.274509803921569 1]);
% create xlabel
xlabel('Motion strength (% coh)');
% Set the rest of the axes properties
set(subplot8,'FontSize',16,'LineWidth',2,'TickDir','out','YTick',...
    [0.5 0.75 1],'YTickLabel',{'50','75','100'});

% create textbox
annotation(figure1,'textbox',...
    [0.0759228187919462 0.477157577925194 0.0430006285845223 0.0515288797132291],...
    'String',{'(e)'},...
    'FontSize',16,...
    'EdgeColor','none');

% create textbox
annotation(figure1,'textbox',...
    [0.0771812080536912 0.252341033562286 0.0430006285845223 0.0515288797132291],...
    'String',{'(g)'},...
    'FontSize',16,...
    'EdgeColor','none');

% create textbox
annotation(figure1,'textbox',...
    [0.0738255033557046 0.690633228358816 0.0423728821790974 0.0515288797132292],...
    'String',{'(c)'},...
    'FontSize',16,...
    'EdgeColor','none');

% create textbox
annotation(figure1,'textbox',...
    [0.0746644295302013 0.90811154723807 0.0430006285845223 0.0515288797132292],...
    'String',{'(a)'},...
    'FontSize',16,...
    'EdgeColor','none');

% create textbox
annotation(figure1,'textbox',...
    [0.506958578308113 0.687297671320791 0.0430006285845224 0.0515288797132292],...
    'String',{'(d)'},...
    'FontSize',16,...
    'EdgeColor','none');

% create textbox
annotation(figure1,'textbox',...
    [0.507215127132366 0.904025667163992 0.0430006285845224 0.0515288797132292],...
    'String',{'(b)'},...
    'FontSize',16,...
    'EdgeColor','none');

% create textbox
annotation(figure1,'textbox',...
    [0.51188352584643 0.464232826295358 0.0386064037465481 0.0515288797132292],...
    'String',{'(f)'},...
    'FontSize',16,...
    'EdgeColor','none');

% create textbox
annotation(figure1,'textbox',...
    [0.506539115220865 0.24500280807863 0.0430006285845224 0.0515288797132291],...
    'String',{'(h)'},...
    'FontSize',16,...
    'EdgeColor','none');

%% plot the DT distribution
% create figure
figure2 = figure;

% create subplot1
subplot1 = subplot(4,4,1);
hold(subplot1,'on');
% create histogram
histogram(DTc1_032,'DisplayName','a','FaceColor',[0 0 1],'NumBins',50);
% Set the rest of the axes properties
set(subplot1,'FontSize',14,'LineWidth',2,'TickDir','out','XTick',...
    [0 0.2 0.4 0.6 0.8 1],'XTickLabel',{'0','0.2','0.4','0.6','0.8','1'});
xlim(subplot1,[0 1]);

% create subplot2
subplot2 = subplot(4,4,2);
hold(subplot2,'on');
% create histogram
histogram(DTe1_032,'DisplayName','b','FaceColor',[1 0 0],'NumBins',50);
% Set the rest of the axes properties
set(subplot2,'FontSize',14,'LineWidth',2,'TickDir','out','XTick',...
    [0 0.2 0.4 0.6 0.8 1],'YTick',[0 25 50]);
ylim(subplot2,[0 50]);
xlim(subplot2,[0 1]);

% create subplot3
subplot3 = subplot(4,4,3);
hold(subplot3,'on');
% create histogram
histogram(DTc1_256,'DisplayName','25.6%, Correct DT','FaceColor',[0 0 1],...
    'NumBins',50);
% Set the rest of the axes properties
set(subplot3,'FontSize',14,'LineWidth',2,'TickDir','out','XTick',...
    [0 0.2 0.4 0.6 0.8 1]);
xlim(subplot3,[0 1]);

% create subplot4
subplot4 = subplot(4,4,4);
hold(subplot4,'on');
% create histogram
histogram(DTe1_256,'DisplayName','d','FaceColor',[1 0 0],'NumBins',50);
xlim(subplot4,[0 1]);
% Set the rest of the axes properties
set(subplot4,'FontSize',14,'LineWidth',2,'TickDir','out','YTick',[0 1],...
     'XTick',[0 0.2 0.4 0.6 0.8 1]);

% create subplot
subplot5 = subplot(4,4,5);
hold(subplot5,'on');
% create histogram
histogram(DTc2_032,'DisplayName','e','FaceColor',[0 0 1],'NumBins',50);
% Set the rest of the axes properties
set(subplot5,'FontSize',14,'LineWidth',2,'TickDir','out','XTick',...
    [0 0.2 0.4 0.6 0.8 1]);
xlim(subplot5,[0 1]);

% create subplot
subplot6 = subplot(4,4,6);
hold(subplot6,'on');
% create histogram
histogram(DTe2_032,'DisplayName','f','FaceColor',[1 0 0],'NumBins',50);
% Set the rest of the axes properties
set(subplot6,'FontSize',14,'LineWidth',2,'TickDir','out','XTick',...
    [0 0.2 0.4 0.6 0.8 1]);
xlim(subplot6,[0 1]);

% create subplot
subplot7 = subplot(4,4,7);
hold(subplot7,'on');
% create histogram
histogram(DTc2_256,'DisplayName','g','FaceColor',[0 0 1],'NumBins',50);
% Set the rest of the axes properties
set(subplot7,'FontSize',14,'LineWidth',2,'TickDir','out','XTick',...
    [0 0.2 0.4 0.6 0.8 1]);
xlim(subplot7,[0 1]);

% create subplot
subplot8 = subplot(4,4,8);
hold(subplot8,'on');
% create histogram
histogram(DTe2_256,'DisplayName','h','FaceColor',[1 0 0],'NumBins',50);
% Set the rest of the axes properties
set(subplot8,'FontSize',14,'LineWidth',2,'TickDir','out','YTick',[0 1],...
    'XTick',[0 0.2 0.4 0.6 0.8 1]);
 ylim(subplot8,[-0.1 1]);
 xlim(subplot8,[0 1]);

% create subplot
subplot9 = subplot(4,4,9);
hold(subplot9,'on');
% create histogram
histogram(DTc3_032,'DisplayName','i','FaceColor',[0 0 1],'NumBins',50);
xlim(subplot9,[0 1]);
% Set the rest of the axes properties
set(subplot9,'FontSize',14,'LineWidth',2,'TickDir','out','XTick',...
    [0 0.2 0.4 0.6 0.8 1],'YTick',[0 25 50 75]);
 ylim(subplot9,[0 75]);

% create subplot
subplot10 = subplot(4,4,10);
hold(subplot10,'on');
% create histogram
histogram(DTe3_032,'DisplayName','j','FaceColor',[1 0 0],'NumBins',50);
% Set the rest of the axes properties
set(subplot10,'FontSize',14,'LineWidth',2,'TickDir','out','XTick',...
    [0 0.2 0.4 0.6 0.8 1]);
xlim(subplot10,[0 1]);

% create subplot
subplot11 = subplot(4,4,11);
hold(subplot11,'on');
% create histogram
histogram(DTc3_256,'DisplayName','k','FaceColor',[0 0 1],'NumBins',50);
% Set the rest of the axes properties
set(subplot11,'FontSize',14,'LineWidth',2,'TickDir','out','XTick',...
    [0 0.2 0.4 0.6 0.8 1],'YTick',[0 100 200]);
xlim(subplot11,[0 1]);

% create subplot
subplot12 = subplot(4,4,12);
hold(subplot12,'on');
% create histogram
histogram(DTe3_256,'DisplayName','l','FaceColor',[1 0 0],'NumBins',50);
% Set the rest of the axes properties
set(subplot12,'FontSize',14,'LineWidth',2,'TickDir','out','YTick',[0 5 10],...
    'XTick',[0 0.2 0.4 0.6 0.8 1]);
xlim(subplot12,[0 1]);

% create subplot
subplot13 = subplot(4,4,13);
hold(subplot13,'on');
% create histogram
histogram(DTc4_032,'DisplayName','m','FaceColor',[0 0 1],'NumBins',50);
% Set the rest of the axes properties
set(subplot13,'FontSize',14,'LineWidth',2,'TickDir','out','XTick',...
    [0 0.2 0.4 0.6 0.8 1]);
xlim(subplot13,[0 1]);

% create subplot
subplot14 = subplot(4,4,14);
hold(subplot14,'on');
% create histogram
histogram(DTe4_032,'DisplayName','n','FaceColor',[1 0 0],'NumBins',50);
% Set the rest of the axes properties
set(subplot14,'FontSize',14,'LineWidth',2,'TickDir','out','XTick',...
    [0 0.2 0.4 0.6 0.8 1]);
xlim(subplot14,[0 1]);

% create subplot
subplot15 = subplot(4,4,15);
hold(subplot15,'on');
% create histogram
histogram(DTc4_256,'DisplayName','o','FaceColor',[0 0 1],'NumBins',50);
% Set the rest of the axes properties
set(subplot15,'FontSize',14,'LineWidth',2,'TickDir','out','XTick',...
    [0 0.2 0.4 0.6 0.8 1],'YTick',[0 50 100 150]);
xlim(subplot15,[0 1]);

% create subplot
subplot16 = subplot(4,4,16);
hold(subplot16,'on');
% create histogram
histogram(DTe4_256,'DisplayName','25.6%, Error DT','FaceColor',[1 0 0],...
    'NumBins',50);
 xlim(subplot16,[0 1]);
 ylim(subplot16,[-0.1 1]);
% Set the rest of the axes properties
set(subplot16,'FontSize',14,'LineWidth',2,'TickDir','out','YTick',[0 1],...
     'XTick',[0 0.2 0.4 0.6 0.8 1]);

% create textbox
annotation(figure2,'textbox',...
    [0.553983016579053 0.949141347424042 0.147189648200565 0.0491981505944517],...
    'String','25.6%, Correct',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.767893247068335 0.94782034346103 0.116053376465832 0.0491981505944517],...
    'String','25.6%, Error',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.295996765062676 0.252632761823454 0.0412454498812008 0.0442536318357273],...
    'String','(n)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.492867026384956 0.910866420764794 0.0412454498812004 0.0442536318357275],...
    'String','(c)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.290262859718289 0.690911772238718 0.0412454498812003 0.0442536318357275],...
    'String','(f)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.494473016916026 0.695917957892194 0.0412454498812009 0.0442536318357272],...
    'String','(g)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.48917691063484 0.250341943250318 0.0412454498812008 0.0442536318357273],...
    'String','(o)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.293845405378082 0.468793642732401 0.0412454498812009 0.0442536318357273],...
    'String','(j)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.293127443051623 0.931841477454137 0.0412454498812003 0.0442536318357275],...
    'String','(b)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.0806274430516234 0.46018614865595 0.0412454498812007 0.0442536318357273],...
    'String','(i)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.0782836930516224 0.689777985390645 0.0412454498812003 0.0442536318357281],...
    'String','(e)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.0782836930516232 0.921070502397447 0.0412454498812006 0.0442536318357275],...
    'String','(a)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.0784956699016032 0.257729560560562 0.0412454498812007 0.0442536318357273],...
    'String','(m)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.709533693051621 0.915401568157085 0.0412454498812004 0.0442536318357275],...
    'String','(d)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.710712570764252 0.702434611229002 0.0412454498812008 0.0442536318357269],...
    'String','(h)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.712042138428358 0.471347284180173 0.0412454498812007 0.0442536318357273],...
    'String','(l)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.713083805095025 0.247723180201285 0.0412454498812008 0.0442536318357274],...
    'String','(p)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.355438738374443 0.948480845442536 0.116053376465832 0.0491981505944517],...
    'String','3.2%, Error',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.146324155041108 0.9479139520185 0.116053376465831 0.0491981505944517],...
    'String','3.2%, Correct',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% create textbox
annotation(figure2,'textbox',...
    [0.419605405041107 0.00743776154230969 0.116053376465831 0.0491981505944518],...
    'String','Normalised RTs',...
    'FontWeight','bold',...
    'FontSize',24,...
    'FitBoxToText','on',...
    'EdgeColor','none');

% create textbox
text('Parent',subplot16,'FontWeight','bold','FontSize',24,'Rotation',90,...
    'String','Frequency',...
    'Position',[-4.475747508306 2.5448275862069 1.4210854715202e-14]);

% create textbox
annotation(figure2,'textbox',...
    [0.490611318910903 0.472481071028246 0.0412454498812005 0.0442536318357273],...
    'String','(k)',...
    'FontWeight','bold',...
    'FontSize',16,...
    'FitBoxToText','off',...
    'EdgeColor','none');

toc;
