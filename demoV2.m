% demoV1

%% 1. Set up single sequence for simulation

%1.1 Specify sequence type
seq.name='TSE';

%1.2 Specify relaxation parameters
seq.T1 = 400; seq.T2 = 100; % arbitrary choice of parameters for initial demonstration

%1.4 Specify refocusing flip angle
alpha = 120;

%1.5 Set up RF scheme (assuming first 90 is 90, 0) 
use_y90 = 1; N = 10;

if use_y90 == 1
    seq.rf(:,1) = [90,90]';
    seq.rf(:,2:N) = repmat([0,alpha]',1,N-1);
else 
    seq.rf(:,1:N) = repmat([0,alpha]',1,N); 
end

%1.6 Specify timing of sequence events
esp = 10; 
dt = esp/2; % time evolves in 0.5*esp steps (dt = 0.5*esp -> dk = 1) 
seq.time = [0 dt dt];
seq.events = {'rf','grad','relax'};

for n = 1:N-1
    % Order of operators : T(rf)->S(grad)->E(relax) "TSE",easy to remember!
    seq.events{end+1} = 'rf';
    seq.events{end+1} = 'grad';
    seq.events{end+1} = 'relax';
    seq.events{end+1} = 'grad';
    seq.events{end+1} = 'relax';
    seq.time = [seq.time (2*n-1)*dt 2*n*dt 2*n*dt (2*n+1)*dt (2*n+1)*dt];
end
seq.grad = ones(1,2*N-1);

%1.7 Run EPG for specified sequence
[om_store,echoes] = EPG_custom(seq);
figure,display_epg(om_store, seq, 1)

%1.8 Plot signal for chosen sequence
%TE1 includes all 'echoes' - some of which are not really echoes as they
%coincide with RF pulses
TE1=echoes(:,1); 
Signal1=echoes(:,2);

%TE2 only includes the true echoes (even) 
TE2=echoes(1:2:size(echoes,1),1); %Plot only even values as the odd ones are not echoes (these correspond with the pulses)
Signal2=echoes(1:2:size(echoes,1),2);

figure
subplot(1,2,1)
plot(TE1,Signal1);
ylim([0 max(Signal1)])

subplot(1,2,2)
plot(TE2,Signal2);
ylim([0 max(Signal2)])


%% 2. For the same single sequence, loop through different T2 values 

%Note that sequence parameters have already been specified above -
%therefore only need to re-specify T2 in the loop

%Create figure outside the loop
figure

for k = 10:10:100

%2.1 Specify T2 value within the loop 
seq.T2 = k; 

%2.2 Run EPG for specified sequence
[om_store,echoes] = EPG_custom(seq);

%2.3 Plot signal for chosen sequence
%TE2 only includes the true echoes (even) 
TE2=echoes(1:2:size(echoes,1),1); %Plot only even values as the odd ones are not echoes (these correspond with the pulses)
Signal2=echoes(1:2:size(echoes,1),2);

plot(TE2,Signal2);
ylim([0 max(Signal2)])
hold on 

end

hold off

%% 3. For the same sequence, loop through different effective TE values to get effective decay curves 

%Note that sequence parameters have already been specified above -
%therefore only need to re-specify TE in the loop

%Create figure outside the loop
figure

for k = 1:1:10

%3.1 Specify echo spacing within the loop  
esp = 10*k; 
dt = esp/2; % time evolves in 0.5*esp steps (dt = 0.5*esp -> dk = 1) 

seq.time = [0 dt dt];
seq.events = {'rf','grad','relax'};

for n = 1:N-1
    % Order of operators : T(rf)->S(grad)->E(relax) "TSE",easy to remember!
    seq.events{end+1} = 'rf';
    seq.events{end+1} = 'grad';
    seq.events{end+1} = 'relax';
    seq.events{end+1} = 'grad';
    seq.events{end+1} = 'relax';
    seq.time = [seq.time (2*n-1)*dt 2*n*dt 2*n*dt (2*n+1)*dt (2*n+1)*dt];
end

seq.grad = ones(1,2*N-1);

%3.2 Run EPG for specified sequence
[om_store,echoes] = EPG_custom(seq);

%3.3 Plot signal for chosen sequence 
TE2=echoes(1:2:size(echoes,1),1); %Plot only even values as the odd ones are not echoes (these correspond with the pulses)
Signal2=echoes(1:2:size(echoes,1),2);
espVec(1:numel(TE2)) = esp; %Vector of values with the same echo spacing value to enable 3D plotting 

plot3(TE2,espVec,Signal2,'LineWidth',1);
xlim([-1 max(TE2)])
zlim([0 max(Signal2)])
hold on 

%3.4 Add point for effective TE 
% (ONLY FOR ILLUSTRATION - NEEDS FURTHER CONSIDERATION AS TO WHICH ECHO THE EFFECTIVE TE SHOULD BE AT)
teIndex = round(numel(Signal2)/2); 

teEffVec(k) = TE2(teIndex);
sVec(k) = Signal2(teIndex);

plot3(TE2(teIndex),espVec(teIndex),Signal2(teIndex),'-o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF', 'LineWidth',1);

end

xlabel('TE (ms)')
ylabel('Echo spacing')
zlabel('Signal (au)')

hold off

%3.5 Plot separate decay curve for effective TEs obtained from all of the
%echo trains

figure
plot(teEffVec',sVec','LineWidth',1)
xlabel('Effective TE (ms)')
ylabel('Signal at effective TE')
ylim([0 1])
hold on 
%Plot monoexponential decay for comparison
plot(teEffVec',exp(-teEffVec/seq.T2)','LineWidth',1)
hold off
legend('Effective decay curve', 'True decay curve')

%% 4. Now implement multiple effective TEs as well as multiple different T2 values to obtain the mapping between true T2 and measured T2 for a specified sequence  

% Loop over T2 values

for j = 1:1:10

    seq.T2 = 10*j; 

% Loop over effective TE values 

%Create figure outside the loop
figure

for k = 1:1:10

%4.1 Specify echo spacing within the loop  
esp = 10*k; 
dt = esp/2; % time evolves in 0.5*esp steps (dt = 0.5*esp -> dk = 1) 

seq.time = [0 dt dt];
seq.events = {'rf','grad','relax'};

for n = 1:N-1
    % Order of operators : T(rf)->S(grad)->E(relax) "TSE",easy to remember!
    seq.events{end+1} = 'rf';
    seq.events{end+1} = 'grad';
    seq.events{end+1} = 'relax';
    seq.events{end+1} = 'grad';
    seq.events{end+1} = 'relax';
    seq.time = [seq.time (2*n-1)*dt 2*n*dt 2*n*dt (2*n+1)*dt (2*n+1)*dt];
end

seq.grad = ones(1,2*N-1);

%4.2 Run EPG for specified sequence
[om_store,echoes] = EPG_custom(seq);

%4.3 Plot signal for chosen sequence 
TE2=echoes(1:2:size(echoes,1),1); %Plot only even values as the odd ones are not echoes (these correspond with the pulses)
Signal2=echoes(1:2:size(echoes,1),2);
espVec(1:numel(TE2)) = esp; %Vector of values with the same echo spacing value to enable 3D plotting 

plot3(TE2,espVec,Signal2,'LineWidth',1);
xlim([-1 max(TE2)])
zlim([0 max(Signal2)])
hold on 

%4.4 Add point for effective TE 
% (ONLY FOR ILLUSTRATION - NEEDS FURTHER CONSIDERATION AS TO WHICH ECHO THE EFFECTIVE TE SHOULD BE AT)
teIndex = round(numel(Signal2)/2); 

teEffVec(k) = TE2(teIndex);
sVec(j,k) = Signal2(teIndex); %NOW STORE SIGNAL INDEXED BY T2 AS WELL AS EFFECTIVE TE
sVecGroundTruth(j,k) = exp(-teEffVec(k)/seq.T2);

plot3(TE2(teIndex),espVec(teIndex),Signal2(teIndex),'-o','Color','b','MarkerSize',10,'MarkerFaceColor','#D9FFFF', 'LineWidth',1);

end

xlabel('TE (ms)')
ylabel('Echo spacing')
zlabel('Signal (au)')

hold off

end

%4.5 Plot separate decay curved for effective TEs obtained from all of the
%echo trains, for different T2 
figure
plot3(repmat(teEffVec,10,1)',repmat(10*([1:1:10]'),1,10)', sVec','LineWidth',1)
xlabel('Effective TE (ms)')
ylabel('T2')
zlabel('Signal')
zlim([0 1])

hold on
plot3(repmat(teEffVec,10,1)',repmat(10*([1:1:10]'),1,10)', sVecGroundTruth','LineWidth',1)
hold off

legend('Measured signal', 'Grouth truth monoexponential decay')