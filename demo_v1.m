

%% No relaxation, perfect 180 pulses
[om_store,echoes,seq] = EPGsim_TSE(180,8,10,1,'none')

figure, plot(echoes(:,1),echoes(:,2))

ylim([0 max(echoes(:,2))])

figure,display_epg(om_store, seq, 1)

%% No relaxation, imperfect pulses
[om_store,echoes,seq] = EPGsim_TSE(150,8,10,1,'none')

figure, plot(echoes(:,1),echoes(:,2))

ylim([0 max(echoes(:,2))])

figure,display_epg(om_store, seq, 1)