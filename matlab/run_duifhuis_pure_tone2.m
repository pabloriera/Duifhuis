clear functions
A=1;
tmax = 0.1;

n_osc = 300;

dx = 0.035/n_osc;
x = (0:n_osc-1)*dx;

f_base_exp_map                 = 22507;
kappa_exp_map                  = 65.1;
f_resonance = f_base_exp_map * 10.^( -kappa_exp_map * x );


fs = 1e5;
n_t = ceil(fs*tmax);
t = (0:n_t-1)/fs;

freqs = 200*(1:5);

%stimulus = zeros(length(t),1);
% Estimulo

%for f0 = freqs
%    stimulus = stimulus + A*sin(2*pi*t*f0)';
%end
%sweep
f1=200;
f2=5000;
fi=f1+t*(f2-f1)/tmax;
phase=2*pi*cumsum(fi)/fs;
stimulus=A*cos(phase)';

stimulus = stimulus.*tukeywin(n_t);

stimulus = (stimulus(3:n_t)+stimulus(1:n_t-2)-2*stimulus(2:n_t-1));
stimulus = [stimulus; zeros(n_t-length(stimulus),1)];

tic
[Y_t V_t] = duifhuis(n_osc, n_t, fs, stimulus);
toc

if (sum(sum(isnan(Y_t)))~=0), disp('NANs');end

figure(1245)
subplot(3,1,1)
curva=10*log10(mean(Y_t(2:end,:).^2,2));
semilogx(f_resonance(2:end),curva)
axis([120 20000 min(curva) max(curva)])
set(gca,'Xdir','reverse')
subplot(3,1,2)
imagesc(t,f_resonance, Y_t)
axis xy
subplot(3,1,3)
plot(t,stimulus)


%%
% fs2 = 44100;
% n_t2 = ceil(fs2*tmax);
% t2 = (0:n_t2-1)/fs2;
% 
% f0 = 1000;
% A = 1;
% 
% stimulus_nob = zeros(length(t2),1);
% 
% for f0 = freqs
%     stimulus_nob = stimulus_nob + A*sin(2*pi*t2*f0)';
% end
% 
% 
% stimulus_nob = stimulus_nob.*tukeywin(n_t2);
% 
% figure(1246)
% tic
% X_t = nobili(stimulus_nob,0)';
% toc
% subplot(3,1,1)
% plot(10*log10(mean(X_t.^2,2)))
% %set(gca,'Xdir','reverse')
% subplot(3,1,2)
% imagesc(X_t)
% subplot(3,1,3)
% plot(t2,stimulus_nob)


