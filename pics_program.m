% The example code. This version of the pics_program.m provide the ‘BG’ part of pictures in the paper.
clear;
close all;
%%
% Model coefficents
order = ["$(a)$","$(b)$","$(c)$","$(d)$","$(e)$","$(f)$"];
c0 = 343;
test_i = 6;
n = 100;
m = 100;

% Line color
C = [248 118 109;
    190	190	190]/256;

%% Figure 4
M0 = 0.1;
U0 = M0*c0;
Ustar = 0.03*U0;

% aerofoil parameters
r_lam_h = [8 4 2 1 0.4 0.2]; % ratio of lambda to h
r_h_c = [0.025 0.025 0.05 0.005 0.05 0.05]; % ratio of h to c
kc = 10.^linspace(-1,2,100); % frequencies required
omega = kc*c0;
freq = omega/2/pi;

fig_phi=figure(2);
fig_phi.Position = [720 50 720 720];
pic_phi=tiledlayout("flow",'TileSpacing','compact','Padding','compact');
disp('Phi')
disp('Figure 4...');
for ii = 1:6
    U0 = M0*c0;
    Ustar = 0.03*U0;
    c = 1;
    obsLocation = [0, 0, c];
    h = c*r_h_c(ii);
    lambda = r_lam_h(ii)*h;

    specSawSC = PredBG_s(lambda, h*2, 1, M0, obsLocation, freq,c0,n,m,1);
    specBaseSC = PredBG_s(0.04, 1e-6, 1, M0, obsLocation, freq,c0,n,m,1);

    nexttile(pic_phi);
    semilogx(kc,10*log10(specBaseSC/4e-10), 'Color',C(1,:),'LineStyle','-','LineWidth',2); hold on;
    semilogx(kc,10*log10(specSawSC/4e-10),'Color',C(2,:),'LineStyle','-','LineWidth',2); hold on;
    ylim([-60 30]);
    yticks([-60 -30 0]);
    xlim([0.1 100]);
    text(0.05,28,order(ii),'Interpreter','latex','FontName','Time New Roman');
    if ii == 1
        legend('Baseline','Serrated','Location','southeast');
    end
    if ii > 4
        xlabel('$kc$','Interpreter','latex','FontName','Time New Roman');
    end
    if ii == 3
        ylabel('PSD (dB ref 4\times 10^{-10}Pa^2/Hz)'); % 增大y轴标签字体
    end
end
fontsize(pic_phi,16,"pixels");
fontname(pic_phi,"Times New Roman");
%% Figure 5
point_num = 121; % number of observer point
kc_list = [1 10 20 50]; % frequencies required

disp('Figure 5...');

% aerofoil parameters
lambda = 0.04;
h = 0.1;

% observation location
r_seq = 10.^linspace(-1,1,point_num);
theta = pi/2;
obsLoc.x = r_seq*cos(theta);
obsLoc.y = 0;
obsLoc.z = r_seq*sin(theta);



fig = figure;
pic = tiledlayout(1,2,'TileSpacing','compact','Padding','loose');
fig.Position = [0 50 1200 360];


for ii = 1:2
    kc = kc_list(ii);
    freq2= kc*c0/2/pi;
    freq = zeros(point_num,1)+freq2;
    omega = 2 * pi * freq;
    for jj = 1:point_num
        obsLocation = [0,0,r_seq(jj)];
        specSawSC(jj,:) = PredBG_s(lambda, h*2, 1, 0.1, obsLocation, freq2,c0,100,100,1);
        specBaseSC(jj,:) = PredBG_s(0.04, 1e-6, 1, 0.1, obsLocation, freq2,c0,100,100,1);
    end

    nexttile;
    plot(r_seq,10*log10(specSawSC/(4*10^(-10))),'Color',C(1,:),'LineWidth',2); hold on;
    plot(r_seq,10*log10(specBaseSC/(4*10^(-10))),'Color',C(2,:),'LineWidth',2);hold on;

    if ii == 2
        text(0.043,40,order(ii),'Interpreter','latex','FontName','Times New Roman','FontSize',22);
    else
        text(0.03,40,order(ii),'Interpreter','latex','FontName','Times New Roman','FontSize',22);
        legend('Present Ser', 'Present Base','Location','southwest','FontSize',16);
    end

    % axis settings
    ylim([-40 40]);
    set(gca,'FontSize',22);
    set(gca,'XScale','log');
    xlabel(gca,'$r/c$','Interpreter','latex','FontName','Times New Roman','FontSize',22);
    set(gca,"FontName",'Times New Roman');
end
ylabel(pic, 'PSD (dB ref $4\times 10^{-10}$ Pa$^2$/Hz)','Interpreter','latex','FontName','Times New Roman','FontSize',16);

%% observation points settings
point_num = 30;
theta = linspace(0, point_num*pi/(point_num-1), point_num+1);
r = c;
x_1_seq = r*cos(theta);
x_3_seq = r*sin(theta);
%% Figure 6
m = 30;
c=1;
disp('directivity calculation');
disp('Figure 6...');

% aerofoil parameters
r_lam_h = 0.4;
r_h_c = 0.1;
h = c*r_h_c;
epsilon = 2*h;
lambda = r_lam_h*h;



kc = [1 3 5 10 20 50];
freq = kc/c*c0/2/pi;

M = 0.1; %0.1,0.4
for jj = 1:point_num
    obsLocation = [x_1_seq(jj),0,x_3_seq(jj)];
    PhiBaseline_BG = PredBG_s(lambda, 0.0000001*c, c, M, obsLocation, freq,c0,100,m,1);
    Phi_BG = PredBG_s(lambda, epsilon, c, M, obsLocation, freq,c0,100,m,1);
    Phi_dir_base(jj,:) = PhiBaseline_BG;
    Phi_dir(jj,:) = Phi_BG;
end
fig_dir1=figure;
fig_dir1.Position = [0 50 1200 960];
pic_dir1 = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
for ii = 1:length(kc)
    nexttile;
    polarplot(theta(1:point_num),10*log10(Phi_dir(1:point_num,ii)/(4*10^(-10))),'Color',C(1,:),'LineWidth',2,'LineStyle','-.'); hold on;
    polarplot(theta(1:point_num),10*log10(Phi_dir_base(1:point_num,ii)/(4*10^(-10))),'Color',C(2,:),'LineWidth',2);hold on;
    rlim([-50 30]);
    thetalim([0 180]);
    if ii==1
        legend('Present Ser', 'Present Base','Location',[0.191558088721173,0.718894675925926,0.131315935591338,0.080729166666667],'FontSize',12);
    end
    text(140*pi/180,80,order(ii),'Interpreter','latex','FontName','Time New Roman');
    set(gca,"FontName",'Times New Roman');
    set(gca,'FontSize',22);
end


%% Figure 7
fig = figure;
pic = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
fig.Position = [0 50 1200 960];

kc = 10;
freq2= kc*c0/2/pi;
freq = zeros(point_num+1,1)+freq2;
r_lam_h = [8 4 2 1 0.4 0.2];
r_h_c = [0.025 0.025 0.05 0.005 0.05 0.05];
disp('Figure 7...');
for ii = 1:6
    lambda = r_lam_h(ii)*r_h_c(ii);
    h = r_h_c(ii);
    for jj = 1:point_num
        obsLocation = [x_1_seq(jj),0,x_3_seq(jj)];
        specSawSC(jj,:) = PredBG_s(lambda, h*2, 1, 0.1, obsLocation, freq2,c0,n,m,1);
        specBaseSC(jj,:) = PredBG_s(0.04, 1e-6, 1, 0.1, obsLocation, freq2,c0,n,m,1);
    end

    nexttile;
    polarplot(theta(1:point_num),10*log10(specSawSC(1:point_num)/(4*10^(-10))),'Color',C(1,:),'LineWidth',2); hold on;
    polarplot(theta(1:point_num),10*log10(specBaseSC(1:point_num)/(4*10^(-10))),'Color',C(2,:),'LineWidth',2);hold on;
    rlim([-50 30]);
    thetalim([0 180]);


    text(140*pi/180,80,order(ii),'Interpreter','latex','FontName','Times New Roman','FontSize',22);
    set(gca,"FontName",'Times New Roman');
    set(gca,'FontSize',22);
    if ii == 1
        legend('Present Ser', 'Present Base','Location',[0.191558088721173,0.718894675925926,0.131315935591338,0.080729166666667],'FontSize',12);
    end
end

