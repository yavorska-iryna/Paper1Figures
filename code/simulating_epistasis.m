% modeling our linearity/additivity/epistasis analysis
% to see how various models affect additivity of modulation index
%
%how to model additivity versus interactions, as a way to see what our
%analyses would look like in either case? Our data clearly look like
%additivity, which means pathway independence, but I was wondering what it
%would look like if the alternative hypothesis were true (i.e. that running
%effect was in fact mediated by VIP pathway).
%Feb-March 2020

% good sources that explain epistasis
% https://academic.oup.com/hmg/article/11/20/2463/616080
% https://en.wikipedia.org/wiki/Epistasis



%model the effect of running and laser, and interactions.

%%%%%%%%%%%%%%%%%%%%%%%%
%
%   a is evoked, b is spont
%
%%%%%%%%%%%%%%%%%%%%%%%%
clear, close all
ncells=500;
a=(10*abs(randn(1,ncells))); %random set of firing rates, sound ON
b= abs(a-(.8*a-abs(randn(1,ncells))));  %random set of firing rates, sound OFF (spont)
figure; hist(a); title('evoked'); 
figure; hist(b); title('spont');

soundMI=(a-b)./(a+b);
indx = find(soundMI<0);
soundMI(indx) = NaN;
% soundMI = normrnd(.5,.2,[1,ncells]);
fprintf('\n mean sound-MI %.1f', mean(soundMI))
figure;
hist(soundMI)
title('sound MI')

% model the effect of running.
% rather than modeling running as a direct effect on soundMI,
% the model really should be at the level of firing
% rates. We aren't imagining that running affects MI, but rather that it
% affects evoked and spontaneous, in diverse ways.

runshifta=randn(size(a)); %random effect of running on evoked
runshifta=3*runshifta - .4; %scale and shift so it's diverse and net suppressive
runshiftb=randn(size(b)); %effect of running on spont
runshiftb=2*runshiftb + .4;  %scale and shift
arun=a+runshifta;
asit=a;
brun=b+runshiftb;
bsit=b;
brun(brun<0)=0; %negative firing rates are not allowed
arun(arun<0)=0;


figure
plot(asit, arun, 'o', mean(asit), mean(arun), 'ro')
hold on
yl=ylim; line([0 yl(2)], [0 yl(2)], 'linestyle', '--')
xlabel('evoked FR sitting')
ylabel('evoked FR running')
title('Evoked Activity fig2a')

figure
plot(bsit, brun, 'o', mean(bsit), mean(brun), 'ro')
hold on
yl=ylim; line([0 yl(2)], [0 yl(2)], 'linestyle', '--')
xlabel('spont FR sitting')
ylabel('spont FR running')
title('Spontaneous Activity fig2b')

soundmi_sit=(asit-bsit)./(asit+bsit);
soundmi_run=(arun-brun)./(arun+brun);
soundmi_run(isnan(soundmi_run))=0;
soundmi_run(soundmi_run<0)=NaN;

figure
plot(soundmi_sit, soundmi_run, 'o')
hold on
xlabel('sound modulation index sit')
ylabel('sound modulation index running')
title('Running Sound Modulation Index fig4a')
yl=ylim; line([-1 yl(2)], [-1 yl(2)], 'linestyle', '--')
ylim([-1.5 1.5])
line(xlim, [0 0])
line([0 0],ylim )

%now we model the effects of laser, as independent of the running effect

lasershifta=randn(size(a));  %random effect of laser on evoked
lasershifta=3*lasershifta + .5;   %scale and shift
lasershiftb=randn(size(b));   %effect of laser on spont
lasershiftb=3*lasershiftb + .5;
alaseron=a+lasershifta;
alaseroff=a;
blaseron=b+lasershiftb;
blaseroff=b;
blaseron(blaseron<0)=0; %negative firing rates are not allowed
alaseron(alaseron<0)=0;


figure
plot(alaseroff, alaseron, 'o', mean(alaseroff), mean(alaseron), 'ro')
hold on
yl=ylim; line([0 yl(2)], [0 yl(2)], 'linestyle', '--')
xlabel('evoked FR laser off')
ylabel('evoked FR laser on')
title('Evoked Activity fig3a')

figure
plot(blaseroff, blaseron, 'o', mean(blaseroff), mean(blaseron), 'ro')
hold on
yl=ylim; line([0 yl(2)], [0 yl(2)], 'linestyle', '--')
xlabel('spont FR laser off')
ylabel('spont FR laser on')
title('Spontaneous Activity fig3b')

soundmi_laseroff=(alaseroff-blaseroff)./(alaseroff+blaseroff);
soundmi_laseron=(alaseron-blaseron)./(alaseron+blaseron);
soundmi_laseron(isnan(soundmi_laseron))=0;
soundmi_laseron(soundmi_run<0)=NaN;

figure
plot(soundmi_laseroff, soundmi_laseron, 'o')
hold on
xlabel('sound modulation index laseroff')
ylabel('sound modulation index laseron')
title('VIP activation Sound Modulation Index fig4b')
ylim([-2 2])
yl=ylim; line([-1 yl(2)], [-1 yl(2)], 'linestyle', '--')

%now we do the additivity test, as in Fig 4c.d

% first we compute the running effect and laser effect in the data
running_effect=soundmi_run-soundmi_sit;
laser_effect=soundmi_laseron-soundmi_laseroff;
figure
plot(running_effect, laser_effect, 'o')
xlim([-3 3])
ylim([-3 3])
grid on
xlabel('Running Effect (MI diff)')
ylabel('Laser Effect (MI diff)')
[rho, p]=corr(running_effect', laser_effect', 'type', 'spearman');
title(sprintf('\x03C1=%.4f, p=%.4f (fig 4c)', rho, p)) % \x03C1 is unicode for rho

% Next, we model the combined effect. So far, we have modeled running and
% laser, separately. In other words, sit and laseroff are identical. We
% simulate data recorded in the combined condition (running and laser on)
% based on a hypothesis about how the signals interact in the brain. First
% let's look tat the hypothesis that they are independent pathways, modeled
% as an additive combination of laser & running effects on firing rates
% (not sound_MI). Then we compare to the additive prediction, which is the
% sum of sound_MI effects.

figure
subplot1(1,3)
% model the combination (Sound MI of run & laser on)
% Sound MI run laser on = (response run laser on- spont run laser on) / (response run laser on + spont run laser on)
arunlaseron= a+lasershifta+runshifta;
brunlaseron= b+lasershiftb+runshiftb;
soundmi_run_laseron=(arunlaseron-brunlaseron)./(arunlaseron+brunlaseron);
soundmi_run_laseron(soundmi_run<0)=NaN;

combined_effect=soundmi_run_laseron - soundmi_laseroff;
predicted_combined_effect=running_effect + laser_effect;
subplot1(1)
plot(combined_effect, predicted_combined_effect, 'o')
hold on; line([-3 3],[-3 3], 'linestyle', '--')
lsline
xlim([-3 3])
ylim([-3 3])
grid on
xlabel('combined effect (running laser-on trials)')
ylabel('additive predicted combined effect (running effect + laser effect)')
[rho, p]=corr(combined_effect', predicted_combined_effect', 'type', 'spearman');
% [rho, p]=corrcoef(combined_effect, predicted_combined_effect); %corrcoef
% doesn't let you do spearman
title(sprintf('additive model, \x03C1=%.2f, p=%.1g (fig 4d)', rho, p))
% It's interesting that even though the effects are additive, by
% construction, there is still a lot of scatter. The rho is not that high.


%so what's the alternative model?
% -either laser and running effects are linearly dependent, instead of
% independent
% -or, the combination is not just adding them, but incorporating
% saturation or some such nonlinearity.

% Model the interaction as sub-additive, i.e. y= a + b - a*b
% (as in Cordell 2002, link at top)
arunlaseron= (a+runshifta+lasershifta)-runshifta.*lasershifta;
brunlaseron= (b+runshiftb+lasershiftb)-runshiftb.*lasershiftb;
soundmi_run_laseron=(arunlaseron-brunlaseron)./(arunlaseron+brunlaseron);

combined_effect=soundmi_run_laseron - soundmi_laseroff;
predicted_combined_effect=running_effect + laser_effect;
subplot1(2)
plot(combined_effect, predicted_combined_effect, 'o')
hold on; line([-3 3],[-3 3], 'linestyle', '--')
lsline
xlim([-3 3])
ylim([-3 3])
grid on
xlabel('sub-additive combined effect (running laser-on trials)')
ylabel('additive predicted combined effect (running effect + laser effect)')
[rho, p]=corr(combined_effect', predicted_combined_effect', 'type', 'spearman');
title(sprintf('sub-addititive model of epistasis, \x03C1=%.2f, p=%.4f (fig 4d)', rho, p))
% rho is small and not significant

% Model the interaction as logical OR (for scalars we will use MAX)
% the idea that the effect of one precludes observing any additional effect of the
%other.
% modeled actual effect:
arunlaseron= a+max(lasershifta,runshifta);
brunlaseron= b+max(lasershiftb+runshiftb);
soundmi_run_laseron=(arunlaseron-brunlaseron)./(arunlaseron+brunlaseron);

combined_effect=soundmi_run_laseron - soundmi_laseroff;
predicted_combined_effect=running_effect + laser_effect; %additive
subplot1(3)
plot(combined_effect, predicted_combined_effect, 'o')
hold on; line([-3 3],[-3 3], 'linestyle', '--')
lsline
xlim([-3 3])
ylim([-3 3])
grid on
xlabel('actual combined effect (MAX) (running laser-on trials)')
ylabel('additive predicted combined effect (running effect + laser effect)')
[rho, p]=corr(combined_effect', predicted_combined_effect', 'type', 'spearman');
title(sprintf('MAX model of epistasis laser=run, \x03C1=%.2f, p=%.1g (fig 4d)', rho, p))
% the correlation is as strong and significant as the independent model
% but the scatterplot does not look as tightly correlated

set(gcf, 'pos',[30 535 1560 420]) % adjust for your screen

% this is from Fig 4 legend
%
% Sound MI sit laser off  = (response sit laser off - spont sit laser off) / (response sit laser off + spont sit laser off)
%
% Sound MI run laser off = (response run laser off - spont run laser off) / (response run laser off + spont run laser off)
%
% Sound MI sit laser on  = (response sit laser on - spont sit laser on) / (response sit laser on + spont sit laser on)
%
% Sound MI run laser on = (response run laser on- spont run laser on) / (response run laser on + spont run laser on)
%
% Running Effect = Sound MI run laser off - Sound MI sit laser off
%
% Laser Effect = Sound MI sit laser on - Sound MI sit laser off
%
% Combined Effect  = Sound MI run laser on - Sound MI sit laser off
%
% Predicted combined Effect = Running Effect + Laser Effect

 f=findobj('type', 'figure');
 for i=1:length(f); figure(f(i)); end


