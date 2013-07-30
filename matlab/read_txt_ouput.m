load X_t.txt

n_osc = size(X_t,1)/2;

Y_t = X_t(1:n_osc,:);V_t = X_t(n_osc+1:end,:);


figure(1453)

subplot(3,1,1)
imagesc(Y_t)
colorbar

load stimulus.txt

subplot(3,1,2)
plot(stimulus)
xlim([0 length(stimulus)])

subplot(3,1,3)
plot(10*log10(mean(Y_t.^2,2)))
xlim([0 size(Y_t,1)])

figure(1454)

ph = plot(Y_t(:,1),'.-');
axis([0 size(Y_t,1) min(min(Y_t)) max(max(Y_t))])

for i = 1:size(Y_t,2)
    
    set(ph,'ydata',Y_t(:,i))
    drawnow
end
