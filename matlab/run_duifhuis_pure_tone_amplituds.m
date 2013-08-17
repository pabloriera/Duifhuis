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


figure(12487)
cla
hold on

f0 = 500;
loc_peak = find(f_resonance<f0,1);

stimulus = zeros(length(t),1);

A = logspace(log10(0.001),log10(.5),10);

res1 = zeros(length(A),1);
res2 = zeros(length(A),1);

cmap = jet(length(A));

damping = [];

for i = 1:length(A);
   
    stimulus = A(i)*sin(2*pi*t*f0)';
    
    stimulus = stimulus.*tukeywin(n_t);
    
%     stimulus = (stimulus(3:n_t)+stimulus(1:n_t-2)-2*stimulus(2:n_t-1));
%     stimulus = [stimulus; zeros(n_t-length(stimulus),1)];
    
    rho = 1000;
    
    damping.dv1 = -1;
    damping.dv3 = 0;1E12;
    damping.dvy2 = 1E18;
    damping.dy3 = 0;
    
    % for each oscillator the ascelleration is
    % g = (dv1 V + dv3 V^3 + dvy2 V Y^2 )*omega/20 + (Y + dy3 Y^3)*omega^2
    
%     damping.dy3 = 0;
    
    tic
    [Y_t V_t] = duifhuis(stimulus,fs,n_osc,rho,damping);
    toc
    
    if (sum(sum(isnan(Y_t)))~=0), disp('NANs');end
    
    Ydb = 10*log10(mean(Y_t(2:end,:).^2,2));
    
    semilogx(f_resonance(2:end),Ydb,'color',cmap(i,:))
    set(gca,'Xdir','reverse')
    xlim([min(f_resonance) max(f_resonance)])
    
    res1(i) = Ydb(loc_peak);
    [m l]= max(Ydb);
    res2(i) = Ydb(l);
    
end

plot([f0 f0],[-320 -160],'k')

%%

figure(12488)
cla
semilogx(A,res1,'b.-')
hold on
semilogx(A,res2,'g.-')