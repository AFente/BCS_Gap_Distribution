% -----------------------------------------------------------------------
%                           DESCRIPTION
% -----------------------------------------------------------------------
% This script calculates the BCS density of states using a normalized gap
% distribution function. The DOS is then convoluted with the derivative of
% the Fermi distribution to introduce the effect of temperature. 
% The result is then plotted with the experimental data in order to compare
% with the calculation. Optionally, the DOS and the distribution can be
% plotted as well.
%
% Note: the gap distribution used in the calculation is a linear
%   interpolation of the imput file in the values of Delta.
%
% Note: the normalization of the conductance is done dividing by the value
%   away from the quasi-particle peaks. Whoever, this number changes upon
%   changing the different parameters. (TO FIX)
%
%
% INPUTS:
% -----------------------------------------------------------------------
% GapDistributionFile: Gap distribution functions with the relative weights
%   of the different values of the gap. Two column ASCII txt file.
%
% ExperimentalDataFile: Normalzied conductance experimental data points as
%   a two column txt file in ASCII.
%
% Temperatura: Temperature for the calculation in K
%
% Energia: Column vector with the values of the energy in meV to perfom the
%   calculation
%
% Delta: Values of the gap in mV to extrapolate the input distribution and
%   perform the calculation.
% -----------------------------------------------------------------------
%
% OUTPUTS:
% -----------------------------------------------------------------------
% SaveFile: two column ASCII file with the result of the calculation.
% -----------------------------------------------------------------------

% CUSTOM OPTIONS:
% -----------------------------------------------------------------------
    TipoFuente          = 'Arial';
    TamahoFuenteEjes    = 14;
    TamahoFuenteTitulos = 16;
    TamanhoLinea        = 2;
    TamanhoPuntos       = 10;
% -----------------------------------------------------------------------

% -----------------------------------------------------------------------
%                       INPUT PARAMETERS
% -----------------------------------------------------------------------
    Temperatura = 0.1;              % K
    Energia     = (-5:0.005:5)';    % meV
    Delta       = (0:0.02:1.5);     % meV
    GapDistributionFile  = 'GapDistribution.txt';
    ExperimentalDataFile = 'ExperimentalData.txt';
    SaveFile             = 'NormalizedBCSDOS.txt';
% -----------------------------------------------------------------------

% Constants:
% -----------------------------------------------------------------------
    kB   = 8.617e-2; % Boltzmann constant in meV/K
    Beta = 1/(kB*Temperatura);
% -----------------------------------------------------------------------
% 
% Data loading
% -----------------------------------------------------------------------
    DistribucionGaps = load(GapDistributionFile); 
        Gamma = interp1(DistribucionGaps(:,1),DistribucionGaps(:,2),Delta,'linear','extrap');
        clear DistribucionGaps;
    
    Curva = load(ExperimentalDataFile);

% Normalization of the gap distribution to the integral    
% -----------------------------------------------------------------------
    Normalizacion = trapz(Delta,Gamma);
    Gamma = Gamma/Normalizacion;
        clear Normalizacion;

% Control
% -----------------------------------------------------------------------
    PasoVoltaje = abs(Energia(2)-Energia(1));
    PasoDistribucion = abs(Delta(2)-Delta(1));

    if PasoVoltaje > PasoDistribucion
        fprintf('\n');
        fprintf('El paso en energia es mayor que el paso de la distribución \n');
    end
% -----------------------------------------------------------------------


% Fermi distribution
% -----------------------------------------------------------------------
    FermiDist	= 1./(1+exp(Energia*Beta));
    dFermiDist	= (Beta*exp(Beta*Energia))./((1+exp(Energia*Beta)).^2); % Analitic expression
%     dFermiDist	= -diff(FermiDist); % Numerical expression

% Fermi distribution representation (control plot)
% -----------------------------------------------------------------------
%     Fig2 = figure(568);
%         Fig2.Color = [1 1 1];
%         Ejes2 = axes('Parent',Fig2,'Box','on');
%         hold(Ejes2,'on');
% 
%         Ejes2_h1 = plot(Energia,FermiDist,'-','Parent',Ejes2);
%             Ejes2_h1.Color = [0 0.4470 0.7410];
%             Ejes2_h1.LineWidth = TamanhoLinea;
%         Ejes2_h2 = plot(Energia,dFermiDist,'-','Parent',Ejes2);
%             Ejes2_h2.Color = [0.8500 0.3250 0.0980];
%             Ejes2_h2.LineWidth = TamanhoLinea;
% 
%         xlabel(Ejes2,'Energia (meV)',...
%             'FontSize',TamahoFuenteTitulos,...
%             'FontName',TipoFuente);
%         title(Ejes2,'Fermi distribution',...
%             'FontSize',TamahoFuenteTitulos,...
%             'FontName',TipoFuente);
%         Ejes2.FontSize = TamahoFuenteEjes;
%         hold(Ejes2,'off');
% 
%         Ejes2.XLim = [min(Energia) max(Energia)];
%         Ejes2.YLim = [min(dFermiDist) max(dFermiDist)];
% -----------------------------------------------------------------------


% Calculation of the BCS DOS using the gap distribution used in 
% P. Martínez-Samper et al., Phys. C 385, 233-243 (2003).
% http://www.sciencedirect.com/science/article/pii/S0921453402022967
% -----------------------------------------------------------------------
    DOS = Energia;
        DOS(:) = 0;
    DOS_AUX = Delta;
        DOS_AUX(:) = 0;
   
    for j=1:length(Energia)
        for i = 1: length(Delta)
            if abs(Energia(j)) > (1+PasoVoltaje/2)*Delta(i)
                DOS_AUX(i) = Gamma(i)*((abs(Energia(j)))/sqrt(Energia(j)^2-Delta(i)^2));
            else
                DOS_AUX(i) = 0;
            end
            DOS(j) = sum(DOS_AUX);
        end
    end

    clear DOS_AUX i j;
    
    DOS = DOS/DOS(end); 
% -----------------------------------------------------------------------


% Convolution of the calculated DOS and the derivative of Fermi function
% -----------------------------------------------------------------------
    Conductancia = conv(dFermiDist,DOS,'same');
    ConductanciaNormalizada = Conductancia/200; % ¡¡MAGIC NUMBER!!    
% -----------------------------------------------------------------------

% Plotting the result
% -----------------------------------------------------------------------

Fig1 = figure(567);
       Fig1.Position = [150 300 1600 420];
       Fig1.Color = [1 1 1];
       
       Sub1 = subplot(1,3,1,'Parent',Fig1);
           Sub1_h1 = plot(Delta,Gamma,'-',...
               'LineWidth',TamanhoLinea,...
               'Color',[0 0.4470 0.7410]);
               Sub1_h1.Parent = Sub1;
           
           ylabel(Sub1,'Relative weight',...
               'FontSize',TamahoFuenteTitulos,...
               'FontName',TipoFuente)    
           xlabel(Sub1,'\Delta (meV)',...
               'FontSize',TamahoFuenteTitulos,...
               'FontName',TipoFuente);
           title(Sub1,'Gap distribution',...
               'FontSize',TamahoFuenteTitulos,...
               'FontName',TipoFuente);
                
           Sub1.FontSize = TamahoFuenteEjes;
           Sub1.FontName = TipoFuente;
           Sub1.Box = 'on';

    Sub2 = subplot(1,3,2,'Parent',Fig1);
        Sub2_h1 = plot(Energia,DOS,'-',...
            'LineWidth',TamanhoLinea,...
            'Color',[0 0.4470 0.7410]);
            Sub2_h1.Parent = Sub2;
    
        xlabel(Sub2,'Energy (meV)');
        ylabel(Sub2,'DOS',...
            'FontSize',TamahoFuenteTitulos,...
            'FontName',TipoFuente);
        title(Sub2,'Normalized density of states',...
            'FontSize',TamahoFuenteTitulos,...
            'FontName',TipoFuente);
    
        Sub2.FontSize = TamahoFuenteEjes;
        Sub2.FontName = TipoFuente;
        Sub2.Box = 'on';
       
    Sub3 = subplot(1,3,3,'Parent',Fig1);
        hold(Sub3,'on');
        Sub3_h1 = plot(Curva(:,1),Curva(:,2),...
            'o','MarkerSize',TamanhoPuntos,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',[0.8500 0.3250 0.0980]);
            Sub3_h1.Parent = Sub3;
            
        Sub3_h2 = plot(Energia,ConductanciaNormalizada,...
            '-','LineWidth',3,'Color',[0 0.4470 0.7410]);
            Sub3_h2.Parent = Sub3;

        xlabel(Sub3,'Energy (meV)',...
            'FontSize',TamahoFuenteTitulos,...
            'FontName',TipoFuente);
            Sub3.XLim = [min(Curva(:,1)) , max(Curva(:,1))];
        ylabel(Sub3,'Normalized conductance',...
            'FontSize',TamahoFuenteTitulos,...
            'FontName',TipoFuente);
            Sub3.YLim = [0 max(Curva(:,2))+0.2];
        title(Sub3,'Normalized conductance',...
            'FontSize',TamahoFuenteTitulos,...
            'FontName',TipoFuente);

        Sub3.FontSize = TamahoFuenteEjes;
        Sub3.FontName = TipoFuente;
        Sub3.Box = 'on';
        hold(Sub3,'off');
   
% Saving data
% -----------------------------------------------------------------------
    dlmwrite(SaveFile, [Energia,Conductancia],...
        'delimiter','\t','newline','pc');
