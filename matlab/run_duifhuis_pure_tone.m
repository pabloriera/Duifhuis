tmax = 0.1;

n_osc = 300;

dx = 0.035/n_osc;
x = (0:n_osc-1)*dx;

f_base_exp_map                 = 22507;
kappa_exp_map                  = 65.1;
f_resonance = f_base_exp_map * 10.^( -kappa_exp_map * x );

A = 1;
fs = 1e5;
n_t = ceil(fs*tmax);
t = (0:n_t-1)/fs;

freqs = 500*(1:1);

stimulus = zeros(length(t),1);

for f0 = freqs
    stimulus = stimulus + A*sin(2*pi*t*f0)';
end

stimulus = stimulus.*tukeywin(n_t);

stimulus = (stimulus(3:n_t)+stimulus(1:n_t-2)-2*stimulus(2:n_t-1));
stimulus = [stimulus; zeros(n_t-length(stimulus),1)];

rho = 1000;

damping.dv1 = 1;
damping.dv3 = 1E12;
damping.dy1 = 0;
damping.dy3 = 0;

tic
[Y_t V_t] = duifhuis(stimulus,fs,n_osc,rho,damping);
toc

if (sum(sum(isnan(Y_t)))~=0), disp('NANs');end

figure(1245)
subplot(3,1,1)
semilogx(f_resonance(2:end),10*log10(mean(Y_t(2:end,:).^2,2)))
set(gca,'Xdir','reverse')
xlim([min(f_resonance) max(f_resonance)])
subplot(3,1,2)
imagesc(t,f_resonance, Y_t)
axis xy
subplot(3,1,3)  
plot(t,stimulus)


%% NOBILI


fs2 = 44100;
n_t2 = ceil(fs2*tmax);
t2 = (0:n_t2-1)/fs2;

f0 = 1000;
A = 1;

stimulus_nob = zeros(length(t2),1);

for f0 = freqs
    stimulus_nob = stimulus_nob + A*sin(2*pi*t2*f0)';
end


stimulus_nob = stimulus_nob.*tukeywin(n_t2);

figure(1246)
tic
X_t = nobili(stimulus_nob,0)';
toc
subplot(3,1,1)
plot(10*log10(mean(X_t.^2,2)))
%set(gca,'Xdir','reverse')
subplot(3,1,2)
imagesc(X_t)
subplot(3,1,3)
plot(t2,stimulus_nob)


