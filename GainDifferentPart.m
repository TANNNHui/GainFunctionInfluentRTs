close all; clear all;

tic;
%Parameter setting
z1= 1; % Upper decision threshold (choice 1)
z2= -1; % Lower decision threshol (choice 2)
dt=1; % Time step 5ms

trials=1000000; % Trial number

k=9.66*10^(-3); % the proportionality factor to form the mean of the drift rate, k=9.66*10^(-3)when both, 1*10^(-3)
Mo_Strength=[0:0.01:0.03 0.032 0.04:0.01:0.51 0.512 0.52:0.01:1]; % the the motion strength
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

AveDTc_onlyDrift=zeros(1,length(Mo_Strength));
AveDTe_onlyDrift=zeros(1,length(Mo_Strength));
ac3=zeros(1,length(Mo_Strength));

AveDTc_onlyNoise=zeros(1,length(Mo_Strength));
AveDTe_onlyNoise=zeros(1,length(Mo_Strength));
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
        DTc1_032=(DTc1-min(DTc1))/(max(DTc1)-min(DTc1));
        DTe1_032=(DTe1-min(DTe1))/(max(DTe1)-min(DTe1));
        DTc1_032(DTc1_032==0)=[];
        DTe1_032(DTe1_032==0)=[];
    end
    if(c==0.512)
        DTc1_512=(DTc1-min(DTc1))/(max(DTc1)-min(DTc1));
        DTe1_512=(DTe1-min(DTe1))/(max(DTe1)-min(DTe1));
        DTc1_512(DTc1_512==0)=[];
        DTe1_512(DTe1_512==0)=[];      
    end
    AveDTc_initial(j)=dt.*mean(DTc1); %Calculate the average decision time
    AveDTe_initial(j)=dt.*mean(DTe1);
    ac1(j)=length(DTc1)/(length(DTc1)+length(DTe1));
    j=j+1;
end
AveDTc_initial=(AveDTc_initial-min(AveDTc_initial))/(max(AveDTc_initial)-min(AveDTc_initial));
AveDTe_initial=(AveDTe_initial-min(AveDTe_initial))/(max(AveDTe_initial)-min(AveDTe_initial));

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
        DTc2_032=(DTc2-min(DTc2))/(max(DTc2)-min(DTc2));
        DTe2_032=(DTe2-min(DTe2))/(max(DTe2)-min(DTe2));
        DTc2_032(DTc2_032==0)=[];
        DTe2_032(DTe2_032==0)=[]; 
    end
    if(c==0.512)
        DTc2_512=(DTc2-min(DTc2))/(max(DTc2)-min(DTc2));
        DTe2_512=(DTe2-min(DTe2))/(max(DTe2)-min(DTe2));
        DTc2_512(DTc2_512==0)=[];
        DTe2_512(DTe2_512==0)=[]; 
    end
    AveDTc_both(j)=dt.*mean(DTc2); %Calculate the average decision time
    AveDTe_both(j)=dt.*mean(DTe2);
    ac2(j)=length(DTc2)/(length(DTc2)+length(DTe2));
    j=j+1;
end
AveDTc_both=(AveDTc_both-min(AveDTc_both))/(max(AveDTc_both)-min(AveDTc_both));
AveDTe_both=(AveDTe_both-min(AveDTe_both))/(max(AveDTe_both)-min(AveDTe_both));

%% time-variant gain only on noise term

j=1;
for c=Mo_Strength
    DTc3=zeros(1,trials); % Correct decision time
    DTe3=zeros(1,trials); % Error decision time 
    Mu0=k*c;% Drift rate
    sigma=sigma0*sqrt(1+c);
    for i=1:trials % Trial number
        x3=zeros(1,length(ti)); % initial x value = '0'
        for t=1:length(ti)
            x3(t+1)=x3(t) + dt*Mu0 + sqrt(dt)*sigma*G(t)*randn;% Update x;
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
        DTc3_032=(DTc3-min(DTc3))/(max(DTc3)-min(DTc3));
        DTe3_032=(DTe3-min(DTe3))/(max(DTe3)-min(DTe3));
        DTc3_032(DTc3_032==0)=[];
        DTe3_032(DTe3_032==0)=[]; 
    end
    if(c==0.512)
        DTc3_512=(DTc3-min(DTc3))/(max(DTc3)-min(DTc3));
        DTe3_512=(DTe3-min(DTe3))/(max(DTe3)-min(DTe3));
        DTc3_512(DTc3_512==0)=[];
        DTe3_512(DTe3_512==0)=[]; 
    end
    AveDTc_onlyNoise(j)=dt.*mean(DTc3); %Calculate the average decision time
    AveDTe_onlyNoise(j)=dt.*mean(DTe3);
    ac3(j)=length(DTc3)/(length(DTc3)+length(DTe3));
    j=j+1;
end
AveDTc_onlyNoise=(AveDTc_onlyNoise-min(AveDTc_onlyNoise))/(max(AveDTc_onlyNoise)-min(AveDTc_onlyNoise));
AveDTe_onlyNoise=(AveDTe_onlyNoise-min(AveDTe_onlyNoise))/(max(AveDTe_onlyNoise)-min(AveDTe_onlyNoise));

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
        DTc4_032=(DTc4-min(DTc4))/(max(DTc4)-min(DTc4));
        DTe4_032=(DTe4-min(DTe4))/(max(DTe4)-min(DTe4));
        DTc4_032(DTc4_032==0)=[];
        DTe4_032(DTe4_032==0)=[];
    end
    if(c==0.512)
        DTc4_512=(DTc4-min(DTc4))/(max(DTc4)-min(DTc4));
        DTe4_512=(DTe4-min(DTe4))/(max(DTe4)-min(DTe4));
        DTc4_512(DTc4_512==0)=[];
        DTe4_512(DTe4_512==0)=[];
    end
    AveDTc_onlyDrift(j)=dt.*mean(DTc4); %Calculate the average decision time
    AveDTe_onlyDrift(j)=dt.*mean(DTe4);
    ac4(j)=length(DTc4)/(length(DTc4)+length(DTe4));
    j=j+1;
end
AveDTc_onlyDrift=(AveDTc_onlyDrift-min(AveDTc_onlyDrift))/(max(AveDTc_onlyDrift)-min(AveDTc_onlyDrift));
AveDTe_onlyDrift=(AveDTe_onlyDrift-min(AveDTe_onlyDrift))/(max(AveDTe_onlyDrift)-min(AveDTe_onlyDrift));

%% Plot the 4x2 MotionStrength vs DT and Accuracy figure
figure;
% For Condition1
subplot(4,2,1)
plot(100*Mo_Strength,AveDTc_initial,'Color','b');
hold on;
plot(100*Mo_Strength,AveDTe_initial,'Color','r');
hold off;
title('Condition i');
set(gca,'FontSize',12);
ylabel('Normalized RTs');
legend('Correct Decision','Error Decision');

subplot(4,2,2)
plot(100*Mo_Strength,ac1);
title('Decision accuracy');
set(gca,'FontSize',12);

% For Condition2
subplot(4,2,3)
plot(100*Mo_Strength,AveDTc_both,'Color','b');
hold on;
plot(100*Mo_Strength,AveDTe_both,'Color','r');
hold off;
title('Condition ii');
set(gca,'FontSize',12);
ylabel('Normalized RTs');

subplot(4,2,4)
plot(100*Mo_Strength,ac2);
title('Decision accuracy');
set(gca,'FontSize',12);

%For Condition3
subplot(4,2,5)
plot(100*Mo_Strength,AveDTc_onlyDrift,'Color','b');
hold on;
plot(100*Mo_Strength,AveDTe_onlyDrift,'Color','r');
hold off;
title('Condition iii');
set(gca,'FontSize',12);
ylabel('Normalized RTs');

subplot(4,2,6)
plot(100*Mo_Strength,ac3);
title('Decision accuracy');
set(gca,'FontSize',12);

%For Condition4
subplot(4,2,7)
plot(100*Mo_Strength,AveDTc_onlyNoise,'Color','b');
hold on;
plot(100*Mo_Strength,AveDTe_onlyNoise,'Color','r');
hold off;
title('Condition iv');
set(gca,'FontSize',12);
xlabel('Motion strength [% coh]');ylabel('Normalized RTs');

subplot(4,2,8)
plot(100*Mo_Strength,ac4);
xlabel('Motion strength [% coh]');
title('Decision accuracy');
set(gca,'FontSize',12);

%% plot the DT distribution

figure; 
subplot(4,4,1);
histogram(DTc1_032,50); title('DT under condition i');
legend('3.2%, Correct DT');
subplot(4,4,2);
histogram(DTe1_032,50,'FaceColor','r'); title('DT under condition i');
legend('3.2%, Error DT');
subplot(4,4,3);
histogram(DTc1_512,50); title('DT under condition i');
legend('51.2%, Correct DT');
subplot(4,4,4);
histogram(DTe1_512,50,'FaceColor','r'); title('DT under condition i');
legend('51.2%, Error DT');

subplot(4,4,5);
histogram(DTc2_032,50); title('DT under condition ii');
legend('3.2%, Correct DT');
subplot(4,4,6);
histogram(DTe2_032,50,'FaceColor','r'); title('DT under condition ii');
legend('3.2%, Error DT');
subplot(4,4,7);
histogram(DTc2_512,50); title('DT under condition ii');
legend('51.2%, Correct DT');
subplot(4,4,8);
histogram(DTe2_512,50,'FaceColor','r'); title('DT under condition ii');
legend('51.2%, Error DT');

subplot(4,4,9);
histogram(DTc3_032,50); title('DT under condition iii');
legend('3.2%, Correct DT');
subplot(4,4,10);
histogram(DTe3_032,50,'FaceColor','r'); title('DT under condition iii');
legend('3.2%, Error DT');
subplot(4,4,11);
histogram(DTc3_512,50); title('DT under condition iii');
legend('51.2%, Correct DT');
subplot(4,4,12);
histogram(DTe3_512,50,'FaceColor','r'); title('DT under condition iii');
legend('51.2%, Error DT');

subplot(4,4,13);
histogram(DTc4_032,50); title('DT under condition iv');
legend('3.2%, Correct DT');
subplot(4,4,14);
histogram(DTe4_032,50,'FaceColor','r'); title('DT under condition iv');
legend('3.2%, Error DT');
subplot(4,4,15);
histogram(DTc4_512,50); title('DT under condition iv');
legend('51.2%, Correct DT');
subplot(4,4,16);
histogram(DTe4_512,50,'FaceColor','r'); title('DT under condition iv');
legend('51.2%, Error DT');
toc;
