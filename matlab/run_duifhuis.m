n_osc = 300;
fs    = 1e5;

[Y fs2] = wavread('vowel.wav');
Y = (Y(3:L)+Y(1:L-2)-2*Y(2:L-1));
Y = Y(2000:5000);
n_t2 = length(Y);
n_t = ceil(fs*n_t2/fs2);
t = (0:n_t-1)/fs;
t2 = (0:n_t2-1)/fs2;
stimulus = interp1(t2,Y,t);
n_t = length(stimulus);
dur = n_t/fs
stimulus = stimulus.*tukeywin(length(stimulus))';

%%

rho = 1000;
tic
[Y_t V_t] = duifhuis(stimulus,fs,n_osc,rho);
toc


%%
figure(1453)

subplot(3,1,1)
imagesc(Y_t)
colorbar

subplot(3,1,2)
plot(stimulus)
xlim([0 length(stimulus)])


subplot(3,1,3)
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