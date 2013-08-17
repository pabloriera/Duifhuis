clear functions

n_osc = 300;
fs    = 1e5;


[Y fs2] = wavread('vowel.wav');
Y = Y(2000:8000);

L = length(Y);
Y2 = (Y(3:L)+Y(1:L-2)-2*Y(2:L-1));



n_t2 = length(Y2);
n_t = ceil(fs*n_t2/fs2);
t = (0:n_t-1)/fs;
t2 = (0:n_t2-1)/fs2;
stimulus = interp1(t2,Y2,t);
n_t = length(stimulus);
dur = n_t/fs
stimulus = stimulus.*tukeywin(length(stimulus))';

%%
tic
X_t = duifhuis(n_osc, n_t-30, fs, stimulus);
toc

%%
n_osc = size(X_t,1)/2;

Y_t = X_t(1:n_osc,:);
V_t = X_t(n_osc+1:end,:);

figure(1245)

subplot(3,1,2)
imagesc(Y_t)

subplot(3,1,3)
plot(stimulus)
xlim([0 length(stimulus)])

subplot(3,1,1)
plot(10*log10(mean(Y_t(2:end,:).^2,2)))
xlim( [1 size(Y_t,1)-1 ] )

% figure(1454)
% 
% ph = plot(Y_t(:,1),'.-');
% axis([0 size(Y_t,1) min(min(Y_t)) max(max(Y_t))])
% 
% for i = 1:size(Y_t,2)
%     
%     set(ph,'ydata',Y_t(:,i))
%     drawnow
% end


%%
%nobili
stimulus_nob = Y';
stimulus_nob = stimulus_nob.*tukeywin(length(stimulus_nob))';

figure(1246)
tic
X_t = nobili(stimulus_nob,1)';
toc
subplot(3,1,1)
plot(10*log10(mean(X_t.^2,2)))
subplot(3,1,2)
imagesc(X_t)
subplot(3,1,3)
plot(t2,stimulus_nob)



