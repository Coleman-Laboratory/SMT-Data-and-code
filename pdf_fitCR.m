function [par1, par2, dtau1, dtau2]=pdf_fitCR(dwet, binsize, bleachrate, parg, pfig)

% pdf_fitting.m is an updated version of dwet_expfit.m. It fits the PDF of
% the residence time with 4 different models (1~4 populations). The quality
% of fitting is assessed via F-test. 
%
% updates:
%   09/30/2015 V2.1
%   Bugs fixed. Fitting curve is plotted with higher resolution.
%   09/14/2015 V 2.0
%   Fitting PDF and output fitting parameters.
%   Fix the CCDF zero point error. Note that bin size needs to be > 1.
%   Add plot options: set pfig = 1 for creating plots.
%   06/09/2015 V 1.1
%   Create a table for the fitted parameters and the F-tests.
%   06/04/2015 V 1.0
%   Imported from hdwet_expfit.m. The filtering process is removed.
%   03/13/2019 report the number
%   04/02/2019 report f23. if it is negative, switch to 3 exponetial fits 
%   05/29/2019 report 95% confident 

%dwet = sptana(1).ndwet;
%bleachrate = sptana(1).bleachrate;
%binsize = sptana(1).acqu;

dt_bin = 0:binsize:max(dwet);
[dt_cts, dt_bin] = hist(dwet, dt_bin);
bx = 0:0.01:max(dwet);

pdf1 = dt_cts/sum(dt_cts);
cdf1 = cumsum(dt_cts)/sum(dt_cts);
ccdf1 = 1 - cdf1;           % ccdf: 1 - CDF data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single exponential fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if parg > 0
    
    pdf_exp1 = @(dwet, tau)exppdf(dwet, tau);       

    f01 = 10; lb1 = 0; ub1 = 400; 
    options = statset('MaxIter',3000, 'MaxFunEvals',3000);

    param1 = mle(dwet, 'pdf', pdf_exp1,...
        'start', f01, 'lower', lb1, 'upper', ub1, 'options', options);  % fit to pdf_exp1

    pdf1_fit = pdf_exp1(dt_bin, param1);
    pdf1_fitc = pdf_exp1(bx, param1);
    ccdf1_fit = pdf1_fitc/pdf1_fitc(1); %1 - exp(-dt_bin/param1); % ccdf_fit: 1 - CDF fitting

    chis_p1 = sum((pdf1(2:end) - pdf1_fit(2:end)).^2);
    %chis_c1 = sum((ccdf1(2:end) - ccdf1_fit(2:end)).^2);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% double exponential fitting %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if parg > 1
    pdf_exp2 = @(dwet, p1, tau1, tau2)...
        p1*exppdf(dwet, tau1) + (1-p1)*exppdf(dwet, tau2); 

    f02 = [0.8,  2,  30];    % initial fitting values [population for short time, t-short, t-long]
    lb2 = [0,  0,  0];       % lower bounds 
    ub2 = [1,  400,  400];   % upper bounds

    options = statset('MaxIter',3000, 'MaxFunEvals', 3000);   % maximum iteration

    param2 = mle(dwet, 'pdf', pdf_exp2, 'start', ...
        f02, 'lower', lb2, 'upper', ub2, 'options', options);  % fit to pdf_exp1

    pdf2_fit = pdf_exp2(dt_bin, param2(1), param2(2), param2(3));
%     cdf2_fit =  param2(1) * (1-exp(-dt_bin/param2(2))) +...
%         (1-param2(1)) * (1-exp(-dt_bin/param2(3)));
    cdf2_fitc = param2(1) * (1-exp(-bx/param2(2))) +...
        (1-param2(1)) * (1-exp(-bx/param2(3)));
    ccdf2_fit = 1 - cdf2_fitc;

    chis_p2 = sum((pdf1(2:end) - pdf2_fit(2:end)).^2);
    %chis_c2 = sum((ccdf1(2:end) - ccdf2_fit(2:end)).^2);     %/(length(dwet_bin)-1);

    %F12_c = ((chis_c1 - chis_c2)/2)/(chis_c2/(length(dt_bin)-1 - 3));
    F12_p = ((chis_p1 - chis_p2)/2)/(chis_p2/(length(dt_bin)-1 - 3));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fitting plot & parameter table %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if pfig == 1


    figure('name', 'CCDF'); clf;
    bar(dt_bin, [1, ccdf1(1:end-1)], 0.75, 'FaceColor', 'b', 'LineStyle', 'none');hold on;
    
    if parg > 1
        plot(bx, ccdf1_fit, 'Color', [0.5, 0.5, 0.5], 'LineStyle', '--', 'LineWidth', 3);
        plot(bx, ccdf2_fit, '-r', 'LineWidth', 3)

      
    hold off;
    xlim([-2, 100]); ylim([0, 1.1]);

end

    Xtx = 2; Ytx = 172; Htx = 10.8; sp = '     ';
    figure; clf; 
    if parg > 0
        avar1 = mlecov(param1, dwet, 'pdf', pdf_exp1); stderr1 = sqrt(diag(avar1));
        stau1 = (1/param1 - 1/bleachrate)^-1;

        text(Xtx, Ytx, '1\circ', 'FontSize', 14, 'FontWeight', 'bold'); hold on;
        text(Xtx, Ytx-Htx, ['\tau_1: ',num2str(stau1), ' \pm ', num2str(stderr1(1))], ...
            'FontSize', 14);
    end
    if parg > 1
        avar2 = mlecov(param2, dwet, 'pdf', pdf_exp2); stderr2 = sqrt(diag(avar2));
        dtau1 = (1/param2(2) - 1/bleachrate)^-1;
        dtau2 = (1/param2(3) - 1/bleachrate)^-1;

        text(Xtx, Ytx-3*Htx, '2\circ', 'FontSize', 14, 'FontWeight', 'bold');
        text(Xtx, Ytx-4*Htx, ['p_1: ',num2str(param2(1)), ' \pm ', num2str(stderr2(1)), sp, ...
            '\tau_1: ', num2str(dtau1), ' \pm ', num2str(stderr2(2))], 'FontSize', 14);
        text(Xtx, Ytx-5*Htx, ['p_2: ',num2str(1-param2(1)), ' \pm ', num2str(stderr2(1)), sp, ...
            '\tau_2: ', num2str(dtau2), ' \pm ', num2str(stderr2(3))], 'FontSize', 14);
        text(Xtx, Ytx-6*Htx, ['F_2_1: ', num2str(F12_p)], 'FontSize', 14, 'FontWeight', 'bold');
        par1 = param2(1);
        par2 = 1-param2(1);
    end
    
    hold off; axis off;
    xlim([0, 180]); ylim([0, 180]);

end





