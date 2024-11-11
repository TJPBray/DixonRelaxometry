%EPG_TSE_testing
clc,clear
%% PPE EPG TSE simulations

% T2_range=[0:2:150];
% T1_range=1000;
ETL=8;
% ETL=32; %CPMG
nograd=1+2*ETL;
test_seq.grad = ones(1,nograd);

test_seq.rf = [[90 0; 90 180], repmat([0;180],1, ETL-1)];
% test_seq.rf = [[90 0; 90 180], repmat([0;180],1, ETL-1)];

% Fill in echo_events variable (seq.events)
echo_events=['rf','grad','relax',repmat({,'rf','grad','relax','grad','relax'},1,ETL)];
test_seq.events = echo_events;

% test_seq.T1 = 1000; test_seq.T2 = 30; %fixed T1 after testing, 0-150ms T2, look into T1

test_seq.name='TSE';

min_T2_ref=[142	191	124	148	84	203	85	359	67	48	151	161];
[min_T2_ref_sort,pos]=sort(min_T2_ref);
min_match_T1_ref=[924,1437,741,1282,750,1199,559,1796,372,249,1750,1563];
min_match_T1_ref_sort(1:12)=min_match_T1_ref(pos);

max_T2_ref=[157	214	136	167	95	224	94	390	73	52	174	184];
max_T2_ref_sort(1:12)=max_T2_ref(pos);
max_match_T1_ref=[745,1160,596,1034,604,966,450,1448,299,200,1415,1262];
max_match_T1_ref_sort(1:12)=max_match_T1_ref(pos);
% ESP_range=[9.5, 11.3, 13.1, 15, 16.8, 18.6]./2;
% ESP = 9.5/2;
% ESP = 9.2/2; % TSE Dixon 20ms
% ESP = 13.8/2; % TSE Dixon 30ms
% ESP = 18.4/2; % TSE Dixon 40ms
% ESP = 27.6/2; % TSE Dixon 60ms
% ESP = 36.8/2; % TSE Dixon 80ms
% ESP = 46/2; % TSE Dixon 100ms

all_ESP=[9.2/2, 13.8/2, 18.4/2, 27.6/2, 36.8/2, 46/2];
% all_ESP=5.4;
% all_ESP=[10/2, 15/2, 20/2, 35/2, 40/2, 50/2];
all_effective_TE=[20,30,40,60,80,100];

% for i=1:6 %effective TE  
for i=1   
ESP=all_ESP(i);
% seq.time (ESP=4.8ms)
echo_timing=zeros(1,length(echo_events));
% echo_timing(1:3)=[0, ESP, ESP];
echo_timing(1:3)=[0, ESP*2, ESP*2];

for repecho=1:ETL
    
    if repecho==1
%         echo_timing(4:8)=[ESP, ESP*2, ESP*2, ESP*3, ESP*3];
        echo_timing(4:8)=[ESP*2, 16.2, 16.2, 16.2+ESP/2, 16.2+ESP/2];
%     elseif repecho==2
%         newrf_pos=4+5*(repecho-1);
%         echo_timing(newrf_pos:newrf_pos+4)=...
%             [ESP*(repecho+1), repmat([ESP*(repecho+2)],1,2), ...
%             repmat([ESP*(repecho+3)],1,2)];
    else
        newrf_pos=4+5*(repecho-1);
        echo_timing(newrf_pos:newrf_pos+4)=...
            [echo_timing(newrf_pos-1), repmat([echo_timing(newrf_pos-1)+ESP/2],1,2), ...
            repmat([echo_timing(newrf_pos-1)+2*ESP/2],1,2)];
%         newrf_pos=4+5*(repecho-1);
%         event1=(5+2*(repecho-3));
%         echo_timing(newrf_pos:newrf_pos+4)=...
%             [ESP*event1, repmat([ESP*(event1+1)],1,2), ...
%             repmat([ESP*(event1+2)],1,2)];
    end
    
end


test_seq.time = echo_timing;
% keyboard %mannually changed EffectiveTE_20ms-100ms
% for mini=1:length(min_T2_ref_sort)
%     test_seq.T1 = min_match_T1_ref_sort(mini);
%     test_seq.T2 = min_T2_ref_sort(mini);
%     [om_store,echoes] = EPG_custom(test_seq);
% %     Allechoes(mini,i).min_echoes.EffectiveTE_100ms=echoes;
%     Allechoes(mini,i).min_echoes=echoes;
% end
% 
% 
% for maxi=1:length(max_T2_ref_sort)
%     test_seq.T1 = max_match_T1_ref_sort(maxi);
%     test_seq.T2 = max_T2_ref_sort(maxi);
%     [om_store,echoes] = EPG_custom(test_seq);
% %     Allechoes(maxi,i).max_echoes.EffectiveTE_100ms=echoes;
%     Allechoes(maxi,i).max_echoes=echoes;
% end

test_seq.T1=1000;
test_seq.T2=[90.1636883046841];

end

% v = round((TE_lowfa(2)/2)+TE_eff(ESPno)); %gap between rf pulses + effective TE
% [~,pos] = (min(abs(TE_lowfa - v))); %finding the position where effective TE is
% Vq = echoes_lowfa(pos);
% New_vals(RFno).echo(count,ESPno)=Vq;
[om_store,echoes] = EPG_custom(test_seq);
figure,display_epg(om_store, test_seq, 1)

TE=echoes(:,1);
Signal=echoes(:,2);
index = Signal>0.1;



% TE_new=TE(1);
% TE_new(2:33)=TE(2:2:65);
% figure,plot(TE_new,Signal)

figure,plot(TE,Signal);
ylim([0 max(Signal)])

figure,plot(TE(index==1),Signal(index==1));
ylim([0 max(Signal)])



%% Test sequence

% T2_range=[0:2:150];
% T1_range=1000;
ETL=7;
% ESP_range=[9.5, 11.3, 13.1, 15, 16.8, 18.6]./2;
% ESP = 9.5/2;
% ESP = 9.2/2; % TSE Dixon 20ms
ESP = 7.5/2;
% ESP = 13.8/2; % TSE Dixon 30ms
% ESP = 18.4/2; % TSE Dixon 40ms
% ESP = 27.6/2; % TSE Dixon 60ms
% ESP = 36.8/2; % TSE Dixon 80ms
% ESP = 46/2; % TSE Dixon 100ms



nograd=1+2*ETL;
test_seq.grad = ones(1,nograd);

test_seq.rf = [[90 0; 90 180], repmat([0;180],1, ETL-1)];

% Fill in echo_events variable (seq.events)
echo_events=['rf','grad','relax',repmat({,'rf','grad','relax','grad','relax'},1,ETL)];
test_seq.events = echo_events;


% seq.time (ESP=4.8ms)
echo_timing=zeros(1,length(echo_events));
echo_timing(1:3)=[0, ESP, ESP];
for repecho=1:ETL
    
    if repecho==1
        echo_timing(4:8)=[ESP, ESP*2, ESP*2, ESP*3, ESP*3];
    elseif repecho==2
        newrf_pos=4+5*(repecho-1);
        echo_timing(newrf_pos:newrf_pos+4)=...
            [ESP*(repecho+1), repmat([ESP*(repecho+2)],1,2), ...
            repmat([ESP*(repecho+3)],1,2)];
    else
        newrf_pos=4+5*(repecho-1);
        event1=(5+2*(repecho-3));
        echo_timing(newrf_pos:newrf_pos+4)=...
            [ESP*event1, repmat([ESP*(event1+1)],1,2), ...
            repmat([ESP*(event1+2)],1,2)];
    end
    
end

test_seq.time = echo_timing;

% test_seq.T1 = 1000; test_seq.T2 = 30; %fixed T1 after testing, 0-150ms T2, look into T1

test_seq.name='TSE';

min_T2_ref=[142	191	124	148	84	203	85	359	67	48	151	161];
[min_T2_ref_sort,pos]=sort(min_T2_ref);
max_T2_ref=[157	214	136	167	95	224	94	390	73	52	174	184];
max_T2_ref_sort(1:12)=max_T2_ref(pos);
test_seq.T1 = 1000;

for mini=1:length(min_T2_ref_sort)
    test_seq.T2 = min_T2_ref_sort(mini);
    [om_store,echoes] = EPG_custom(test_seq);
    Allechoes(mini).min_echoes=echoes;
end


for maxi=1:length(max_T2_ref_sort)
    test_seq.T2 = max_T2_ref_sort(maxi);
    [om_store,echoes] = EPG_custom(test_seq);
    Allechoes(maxi).max_echoes=echoes;
end


[om_store,echoes] = EPG_custom(test_seq);
figure,display_epg(om_store, test_seq, 1)
%% Our sequence

clc,clear
tic
MyRF_start(1).RF_pulses = [90 0; 90 180];
MyRF_start(2).RF_pulses = [90 0; 90 150];
MyRF_start(3).RF_pulses = [90 0; 90 120];
MyRF_start(4).RF_pulses = [0 0; 90 180];

T2_range=[0:2:150];
% T2_range=[0, 130:2:150];

%200-1800ms T1 range from Marty paper
T1_range=[0,1000]; %T1 much shorter than TR in actual data
% T1_range=1000;

% Scaling_factor=5;
% ETL_range=[2,10*Scaling_factor];
% ESP=4.8/Scaling_factor; %From PPE

ETL=10;
ESP_range=[9.5, 11.3, 13.1, 15, 16.8, 18.6]./2;

% seq.rf (0/90 phase shift and 120/180)

% seq.events
% seq.time (ESP=4.8ms)

% seq.T1 = 1000; seq.T2 = 100; %fixed T1 after testing, 0-150ms T2, look into T1

% seq.grad =[1,1,1,1]; %Configuration state changes (integer) - check
% seq.name='TSE';

%% Test different combinations

keyboard
Myseq=[];
%ETL 20-30 (seq.events, seq.timing and seq.grad)

for ESP_no=1:length(ESP_range)

% for ETL_no=1:length(ETL_range)
    
    seq=[];
    
    %Unchanged parameters
    seq.name='TSE';
%     seq.T1 = T1_range; %Possibly test range

%     current_ETL=ETL_range(ETL_no);
    current_ETL=ETL;
    
%     keyboard
    nograd=1+2*current_ETL;
    seq.grad = ones(1,nograd); %Configuration state changes (integer) - check
%     seq.grad=[1,1,1,1]; 


    ESP=ESP_range(ESP_no);
    
    %90 RF pulse (0/90 phase shift), Refocusing pulse 180 or 120
    for RF_combo=1:4
        
        testno=0;
        seq.rf = [];
        seq.events = [];
        seq.time = [];
        
        current_RF_start=[];
        refocus=[];
        echo_events=[];
        echo_timing=[];
        
        
        % Fill in RF pulses variable (seq.rf)
%         keyboard
        current_RF_start=MyRF_start(RF_combo).RF_pulses;
        refocus=current_RF_start(2,2);
        
%         RF_train=zeros(2,1+current_ETL);
        RF_train=[current_RF_start, repmat([0;refocus],1, current_ETL-1)];
        seq.rf=RF_train;
        
        % Fill in echo_events variable (seq.events)
        echo_events=['rf','grad','relax',repmat({,'rf','grad','relax','grad','relax'},1,current_ETL)];
        seq.events = echo_events;
        
        % Fill in echo_timing variable (seq.time)
        echo_timing=zeros(1,length(echo_events));
        echo_timing(1:3)=[0, ESP, ESP];
        
        for repecho=1:current_ETL
            
            if repecho==1
                echo_timing(4:8)=[ESP, ESP*2, ESP*2, ESP*3, ESP*3];
            elseif repecho==2
                newrf_pos=4+5*(repecho-1);
                echo_timing(newrf_pos:newrf_pos+4)=...
                    [ESP*(repecho+1), repmat([ESP*(repecho+2)],1,2), ...
                    repmat([ESP*(repecho+3)],1,2)];
            else
                newrf_pos=4+5*(repecho-1);
                event1=(5+2*(repecho-3));
                echo_timing(newrf_pos:newrf_pos+4)=...
                    [ESP*event1, repmat([ESP*(event1+1)],1,2), ...
                    repmat([ESP*(event1+2)],1,2)];
            end
            
        end
        
        seq.time = echo_timing;
        
        
        for T1_no=1:length(T1_range)
            for T2_no=1:length(T2_range)
            
                %             keyboard
                seq.T2 = T2_range(T2_no);
                seq.T1 = T1_range(T1_no);
                [om_store,echoes] = EPG_custom(seq);
                testno=testno+1;
                
                Myseq(ESP_no,RF_combo).RFstart(testno).seq=seq;
                Myseq(ESP_no,RF_combo).RFstart(testno).om_store=om_store;
                Myseq(ESP_no,RF_combo).RFstart(testno).echoes=echoes;
                %             keyboard
                %             display_epg(om_store, seq, 1)
            end
            
        end
        Myseq(ESP_no,RF_combo).ESP=ESP;
        
%         keyboard
        
        
    end
%     keyboard
    
end
    
Myseq_info=[];

Myseq_info(1).('ETL')=ETL;
Myseq_info(1).('RF_pulse_start') = MyRF_start;  
Myseq_info(1).('ESP_range') = ESP_range.*2;
Myseq_info(1).('T2_range') = T2_range;
Myseq_info(1).('T1_range') = T1_range;


% [om_store,echoes] = EPG_custom(seq);
% display_epg(om_store, seq, 1)
toc
keyboard
%% Save all sequences and info

cd('C:\Users\rgoll\Google Drive\PostDoc Projects\T2 Project\EPG testing')
save(['Myseq_ESP_', num2str(Myseq_info.ESP_range(1)), '_',...
    num2str(Myseq_info.ESP_range(end)),...
    '_T2_', num2str(Myseq_info.T2_range(1)), '_',...
    num2str(Myseq_info.T2_range(end)),...
    '_T1_', num2str(Myseq_info.T1_range(1)), '_',...
    num2str(Myseq_info.T1_range(end)),...
    '_ETL_', num2str(Myseq_info.ETL),  '.mat'],...
    'Myseq')

save(['Myseq_ESP_', num2str(Myseq_info.ESP_range(1)), '_',...
    num2str(Myseq_info.ESP_range(end)),...
    '_T2_', num2str(Myseq_info.T2_range(1)), '_',...
    num2str(Myseq_info.T2_range(end)),...
    '_T1_', num2str(Myseq_info.T1_range(1)), '_',...
    num2str(Myseq_info.T1_range(end)),...
    '_ETL_', num2str(Myseq_info.ETL), '_info.mat'],...
    'Myseq_info')
    
    
keyboard

%% Load Myseq 

clc,clear
cd('C:\Users\rgoll\Google Drive\PostDoc Projects\T2 Project\EPG testing')
% uiopen
[FILENAME, PATHNAME] = uigetfile

load([PATHNAME, FILENAME])

%% Single TSE train (Plot different T2 and 180 v 150/120)
% close all
ESPno=1;
% T1_T2no=50+76;

% T1_T2no=[1:10:76,77:10:152]
figtot=length([77:10:152]);
% figure
% count=0;

All_echoes(1).T2=[];
All_echoes(1).echo_180=[];
All_echoes(1).echo_150=[];
All_echoes(1).echo_120=[];

for RFno=2:3
    figure
    set(gcf,'color','w')
    count=0;
    for T1_T2no=[87:10:152]
        count=count+1;
        subplot(2,4,count)
        
        current_seq_highfa=Myseq(ESPno,1).RFstart(T1_T2no).seq;
        current_om_store_highfa=Myseq(ESPno,1).RFstart(T1_T2no).om_store;
        current_echoes_highfa=Myseq(ESPno,1).RFstart(T1_T2no).echoes;
        plot(current_echoes_highfa([1:2:end],1),current_echoes_highfa([1:2:end],2))
        
        
        
        current_seq_lowfa=Myseq(ESPno,RFno).RFstart(T1_T2no).seq;
        current_om_store_lowfa=Myseq(ESPno,RFno).RFstart(T1_T2no).om_store;
        current_echoes_lowfa=Myseq(ESPno,RFno).RFstart(T1_T2no).echoes;
        
        hold on
        plot(current_echoes_lowfa([1:2:end],1),current_echoes_lowfa([1:2:end],2))
        
        
        title(['T2 = ', num2str(current_seq_lowfa.T2),...
            ', T1 = ', num2str(current_seq_lowfa.T1)],...
            'FontSize',14)
        
        ylim([0,1])
        
        
%         xlabel([FILENAME(end-10:end-4), ' Time (ms)'])
%         xlabel(['ESP = ', num2str(Myseq(ESPno,1).ESP), ' ms'])
        xlabel('TE', 'FontSize',14)
        
        if RFno==2
            ylabel('Echo intensity (180 v 150)', 'FontSize',14)
            All_echoes(count).echo_150 = [current_echoes_lowfa([1:2:end],1),current_echoes_lowfa([1:2:end],2)];
        else
            ylabel('Echo intensity (180 v 120)', 'FontSize',14)
            All_echoes(count).echo_120 = [current_echoes_lowfa([1:2:end],1),current_echoes_lowfa([1:2:end],2)];
        end
        
        All_echoes(count).echo_180 = current_echoes_highfa;
        All_echoes(count).T2=Myseq(ESPno,1).RFstart(T1_T2no).seq.T2;
        
    end

end

% [om_store,echoes] = EPG_custom(current_seq);
% figure,display_epg(om_store, current_seq, 1)


%% T2 fit - single TSE train
T2_range=[87:10:152];
tic
for T1_T2no=1:length(T2_range)
    
    original_Echo_180=All_echoes(T1_T2no).echo_180(:,2); 
    original_TE_180=All_echoes(T1_T2no).echo_180(:,1);
    
%     original_Echo_150=All_echoes(T1_T2no).echo_150(:,2); 
%     original_TE_150=All_echoes(T1_T2no).echo_150(:,1);
%     
%     original_Echo_120=All_echoes(T1_T2no).echo_120(:,2); 
%     original_TE_120=All_echoes(T1_T2no).echo_120(:,1);
    
    [fitresult_org_180, gof_org_180] = T2Fit_180(original_TE_180, original_Echo_180);
%     [fitresult_org_150, gof_org_150] = T2Fit_180(original_TE_150, original_Echo_150);
%     [fitresult_org_120, gof_org_120] = T2Fit_180(original_TE_120, original_Echo_120);
    
    Predicted_T2_org_180=(1/fitresult_org_180.b)*-1;
%     Predicted_T2_org_150=(1/fitresult_org_150.b)*-1;
%     Predicted_T2_org_120=(1/fitresult_org_120.b)*-1;

    
    T2_fitting_org(1,T1_T2no).T2= All_echoes(T1_T2no).T2;
    T2_fitting_org(1,T1_T2no).Fitted_T2_org_180=Predicted_T2_org_180;
%     T2_fitting_org(1,T1_T2no).Fitted_T2_org_150=Predicted_T2_org_150;
%     T2_fitting_org(1,T1_T2no).Fitted_T2_org_120=Predicted_T2_org_120;
    
end
toc

for i=1:7
    T2(i)=T2_fitting_org(i).T2;
    T2_180(i)=T2_fitting_org(i).Fitted_T2_org_180;
    T2_150(i)=T2_fitting_org(i).Fitted_T2_org_150;
    T2_120(i)=T2_fitting_org(i).Fitted_T2_org_120;
end

figure,plot(T2, T2_180)
set(gcf,'color','w')
hold on,plot(T2, T2_150)
hold on,plot(T2, T2_120)
xlim([0,160])
ylim([0,160])
title('T2 fit @ 180, 150 and 120', 'Fontsize',14)
xlabel('Actual T2', 'Fontsize',14)
ylabel('Fitted T2', 'Fontsize',14)
legend('180','150','120')

%% Load Myseq 

clc,clear
% cd('G:\My Drive\PostDoc Projects\T2 Project\EPG model data')
% uiopen
[FILENAME, PATHNAME] = uigetfile
load([PATHNAME, FILENAME])

%% Plot different effective echo times and T2
my_colours{1}=[0 0 1];
my_colours{2}=[0.9290 0.6940 0.1250];
my_colours{3}=[1 0 0];
% my_colours{3}=[0.4660 0.6740 0.1880]; % or [0 1 0]
my_colours{4}=[0.6350 0.0780 0.1840];
my_colours{5}=[0.8500 0.3250 0.0980];
my_colours{6}=[0.4940 0.1840 0.5560];

% my_colours{1}=[0.8500,0.3250, 0.0980];
TE_eff=[52.27 62.27 72.27 82.27 92.27 102.27];
New_vals(1).echo=[];
for RFno=1:3
%     hold on
    figure
    set(gcf,'color','w')
    for ESPno=1:6%ESPno=1:6 %(size(Myseq),2)
        %     keyboard
%         xlabel('TE','FontSize',14)
%         ylabel('Echo Intensity','FontSize',14)
%         xlabel('T2','FontSize',14)
%         ylabel('TE','FontSize',14)
%         zlabel('Echo Intensity','FontSize',14)
        xlabel('Effective TE','FontSize',14)
        ylabel('TE in TSE train','FontSize',14)
        zlabel('Echo Intensity','FontSize',14)
        
        count=0;
        for T1_T2no=97%T1_T2no=[87:10:152]
            count=count+1;
            %         current_seq_lowfa=Myseq(ESPno,RFno).RFstart(T1_T2no).seq;
            %         current_om_store_lowfa=Myseq(ESPno,RFno).RFstart(T1_T2no).om_store;
            
            current_echoes_lowfa=Myseq(ESPno,RFno).RFstart(T1_T2no).echoes;
            
            current_T2=Myseq(ESPno,RFno).RFstart(T1_T2no).seq.T2;
            
            current_TE_eff=TE_eff(ESPno);
            if RFno==1
                TE_lowfa=current_echoes_lowfa(:,1);
                echoes_lowfa=current_echoes_lowfa(:,2);
                T2_lowfa=current_T2.*ones([size(current_echoes_lowfa,1),1]);
                TE_eff_plot=current_TE_eff.*ones([size(current_echoes_lowfa,1),1]);
            else
                TE_lowfa=current_echoes_lowfa([1:2:end],1);
                echoes_lowfa=current_echoes_lowfa([1:2:end],2);
                T2_lowfa=current_T2.*ones([size(current_echoes_lowfa,1)/2,1]);
                TE_eff_plot=current_TE_eff.*ones([size(current_echoes_lowfa,1)/2,1]);
            end

%             plot(TE_lowfa, echoes_lowfa) , 'Color', 'b'
            plot3(TE_eff_plot,TE_lowfa, echoes_lowfa, 'Color', my_colours{RFno})
%             plot3(T2_lowfa,TE_lowfa, echoes_lowfa)
            set(gca,'Ydir','reverse')
            hold on
            plot3(TE_eff_plot(1),TE_lowfa(6), echoes_lowfa(6),'o',...
                'Color', my_colours{RFno},'MarkerSize',20)
            
%             keyboard 'Color', 'b', 
%             Vq = interp1(TE_lowfa,echoes_lowfa,TE_eff(ESPno));

            v = round((TE_lowfa(2)/2)+TE_eff(ESPno));
            [~,pos] = (min(abs(TE_lowfa - v)));
            Vq = echoes_lowfa(pos);
            New_vals(RFno).echo(count,ESPno)=Vq;
        end
    end
    
    if RFno==1
        pulse_str=' 180 ';
    elseif RFno==2
        pulse_str=' 150 ';
    elseif RFno==3
        pulse_str=' 120 ';
    end
        
        title(['TSE', pulse_str, 'refocusing pulse'], 'FontSize',14)
        xlim([50,105])
%         xlim([0,150])
    
    
end

% title('TSE signal with different refocusing pulses')
% legend('180 refocusing pulse', '150 refocusing pulse', '120 refocusing pulse')

%New_vals 7x6 with 7 T2 values and 6 effective echoes
%% Plot Effective TE only and different T2
figure
for RFno=1:2:3
%     figure
    set(gcf,'color','w')
    count=0;
    T2_range=97 %[87:10:152];
    for T1_T2no=1:length(T2_range)
        % for TE_eff_no=1:length(TE_eff)
       
        current_echo_TEeff=New_vals(RFno).echo(T1_T2no,:);
%         keyboard
        Fitting(RFno,T1_T2no).TE_eff=TE_eff;
        Fitting(RFno,T1_T2no).echo_intensity=current_echo_TEeff;
        plot(TE_eff, current_echo_TEeff, 'LineWidth',2,'Color', my_colours{RFno})
        ylabel('Echo Intensity','FontSize',14)
        xlabel('Effective TE','FontSize',14)
        %     T2_plot=T2_range(T1_T2no).*ones([size(current_echo_TEeff,2),1]);
        %     plot3(T2_plot,TE_eff,current_echo_TEeff)
        %     set(gca,'Xdir','reverse','Ydir','reverse')
        
        hold on
        %     legend
    end
% %     token = split(num2str(T2_range), '  ')
%     lgd=legend({'40'})
% %     lgd = legend({'20','40','60','80','100','120','140'});
%     title(lgd,'T2 (ms)')
    
    if RFno==1
        pulse_str=' 180 ';
    elseif RFno==2
        pulse_str=' 150 ';
    elseif RFno==3
        pulse_str=' 120 ';
    end
    
    title(['TSE', pulse_str, 'refocusing pulse'],'FontSize',14)
    ylim([0,0.7])
end
title('TSE 180 and 120 refocusing pulse')
ax = gca
ax.FontSize = 22; 


lgd=legend([{'40 (180 refocussing pulse)'},{'40 (120 refocussing pulse)'}])
title(lgd,'T2 (ms)')
%% T2 fit
tic
T2_range=[87:10:152];
for T1_T2no=1:length(T2_range)
    
    current_Echo_180=Fitting(1,T1_T2no).echo_intensity;
    current_Echo_150=Fitting(2,T1_T2no).echo_intensity;
    current_Echo_120=Fitting(3,T1_T2no).echo_intensity;
    
    [fitresult_180, gof_180] = T2Fit_180(TE_eff, current_Echo_180);
    [fitresult_150, gof_150] = T2Fit_180(TE_eff, current_Echo_150);
    [fitresult_120, gof_120] = T2Fit_180(TE_eff, current_Echo_120);
    

    Predicted_T2_180=(1/fitresult_180.b)*-1;
    Predicted_T2_150=(1/fitresult_150.b)*-1;
    Predicted_T2_120=(1/fitresult_120.b)*-1;
    
    %T2_range(T1_T2no) is wrong?????
    T2_fitting_results(1,T1_T2no).T2=All_echoes(T1_T2no).T2;
    T2_fitting_results(1,T1_T2no).Fitted_T2_180=Predicted_T2_180;
    T2_fitting_results(1,T1_T2no).Fitted_T2_150=Predicted_T2_150;
    T2_fitting_results(1,T1_T2no).Fitted_T2_120=Predicted_T2_120;
%     T2_fitting_results(1,T1_T2no).Fitted_T2=Predicted_T2_180;
    
%     T2_fitting_results(2,T1_T2no).T2=T2_range(T1_T2no);
%     T2_fitting_results(2,T1_T2no).Fitted_T2=Predicted_T2_150;
%     
%     T2_fitting_results(3,T1_T2no).T2=T2_range(T1_T2no);
%     T2_fitting_results(3,T1_T2no).Fitted_T2=Predicted_T2_120;
end
toc
for i=1:7
    Effective_T2(i)=T2_fitting_results(i).T2;
    T2_180(i)=T2_fitting_results(i).Fitted_T2_180;
    T2_150(i)=T2_fitting_results(i).Fitted_T2_150;
    T2_120(i)=T2_fitting_results(i).Fitted_T2_120;
end

figure,plot(Effective_T2, T2_180)
set(gcf,'color','w')
hold on,plot(Effective_T2, T2_150)
hold on,plot(Effective_T2, T2_120)
xlim([0,160])
ylim([0,160])
title('T2 fit @ 180, 150 and 120', 'Fontsize',14)
xlabel('Actual T2', 'Fontsize',14)
ylabel('Fitted T2', 'Fontsize',14)
legend('180','150','120')

% for i=1:7
%     corr_factor(i)=T2_range(i)/((-1/fittedmodel(i).b))
% end

% current_Echo_180=Fitting(1,1).echo_intensity;
% -1/fittedmodel1.b;


% val(x) = a*exp(b*x)
% S(TE)= kSo.exp(?TE/T2)
% 
% b*TE = TE* -1/T2
% b=-1/T2
% T2=-1/b
% (1/fitresult.b)*-1
%% Junk
keyboard
for RFno=2:3
    figure
    set(gcf,'color','w')
    count=0;
    for T1_T2no=[87:10:152]
        count=count+1;
        subplot(2,4,count)
        
        current_seq_highfa=Myseq(ESPno,1).RFstart(T1_T2no).seq;
        current_om_store_highfa=Myseq(ESPno,1).RFstart(T1_T2no).om_store;
        current_echoes_highfa=Myseq(ESPno,1).RFstart(T1_T2no).echoes;
        plot(current_echoes_highfa([1:2:end],1),current_echoes_highfa([1:2:end],2))
        
        
        current_seq_lowfa=Myseq(ESPno,RFno).RFstart(T1_T2no).seq;
        current_om_store_lowfa=Myseq(ESPno,RFno).RFstart(T1_T2no).om_store;
        current_echoes_lowfa=Myseq(ESPno,RFno).RFstart(T1_T2no).echoes;
        
        hold on
        plot(current_echoes_lowfa([1:2:end],1),current_echoes_lowfa([1:2:end],2))
        
        
        title(['T2 = ', num2str(current_seq_lowfa.T2),...
            ', T1 = ', num2str(current_seq_lowfa.T1)],'FontSize',14)
        
%         xlabel([FILENAME(end-10:end-4), ' Time (ms)'])
        xlabel(['ESP = ', num2str(Myseq(ESPno,1).ESP), ' ms'],'FontSize',14)
        
        if RFno==2
            ylabel('Echo intensity (180 v 150)')
        else
            ylabel('Echo intensity (180 v 120)')
        end
        
        
    end

end

figure,display_epg(om_store, seq, 1)