% ------------------------------------------------------
% CONDUCTANCIA BCS A PARTIR DE UNA DISTRIBUCI�N DE GAPS
% ------------------------------------------------------
%
% Es necesario que la distribuci�n est� guardada en un archivo .m al que se
% llama aqu� distribucion.m y dentro de �l un array de dos clumnas llamado
% DistribucionGaps. Del mismo modo, se usa un archivo llamado
% datosExperimento.m que contiene un ejemplo de una curva de conductancia
% normalizada en funci�n del voltaje bias para compararla con el ajuste.
%
%   distribucion.m
%       DistribucionGaps = [ValorGap, PesoDistribucion];
%
%   datosExperimento.m
%       Curva = [VoltajeBias, ConductanciaNormalizada];
%
% Se gererar� un archivo Fit.txt con la curva ajustada (variables Energia y
% Conductancia en el workspace) que se puede arrastrar directamente a
% Origin.
%
% Notas: 
%
%   Resulta de gran utilidad combinar este script con la funci�n move
% data points de Origin. Hay que copiar la distribuci�n creada en Origin al
% archivo que la contiene y modificarla seg�n las necesidades. Es
% recomendable usar pocos puntos (no m�s de 200) en la distribuci�n durante
% este proceso. Sin embargo, para el c�lculo final es bueno interpolar la
% distribuci�n definitiva usando muchos m�s puntos (de nuevo resulta
% sencillo en Origin.
%
%   La temperatura y los valores de energ�a para los que se quiere
% calcular, es necesario introducirlos a mano.
%
%   La densidad de estados se normaliza a un "magic number" que depende del
% n�mero de puntos utilizado para calcularla (vector Energia). Modificar la
% l�nea 137 si es necesario.
% -----------------------------------------------------------

% ----------------------------------------------------------------------
%                       INTRODUCIR PAR�METROS:
% ----------------------------------------------------------------------
    Temperatura = 0.1;              % K
    kB          = 8.617e-2;         % meV/K
    Energia     = (-5:0.005:5)';    % meV
    Delta = (0:0.02:1.5);
% ----------------------------------------------------------------------

% CARGA DE DATOS
% --------------------
    DistribucionGaps = load('DistribucionIsa_Actual2.txt'); % Archivo que contiene la distribuci�n
        Gamma = interp1(DistribucionGaps(:,1),DistribucionGaps(:,2),Delta,'linear','extrap');
%         Delta = DistribucionGaps(:,1);
%         Gamma = DistribucionGaps(:,2);
        clear DistribucionGaps;
        PasoDistribucion = abs(Delta(2)-Delta(1));
    
    Curva = load('GapIsaExp.txt'); % Archivo que contiene los puntos a ajustar
    
    Normalizacion = trapz(Delta,Gamma);

    Fig1 = figure(567);
        Fig1.Position = [150 300 1600 420];
        subplot(1,3,1,'Parent',Fig1);
        plot(Delta,Gamma/Normalizacion,'-','LineWidth',2,'Color',[0 0.4470 0.7410]);
        xlabel('\Delta (meV)','FontSize',14);
        title('Gap distribution','FontSize',14);
        set(gcf,'Color',[1 1 1]);
        set(gca,'FontSize',12);
% --------------------



Beta = 1/(kB*Temperatura);
PasoVoltaje = abs(Energia(2)-Energia(1));

% CONTROL
% -----------------------------------
    if PasoVoltaje > PasoDistribucion
        fprintf('\n');
        fprintf('El paso en energia es mayor que el paso de la distribuci�n \n');
    end
% -----------------------------------


% DISTRIBUCI�N DE FERMI
% ----------------------------------
    FermiDist	= 1./(1+exp(Energia*Beta));
    dFermiDist	= (Beta*exp(Beta*Energia))./((1+exp(Energia*Beta)).^2); % Anal�tica
%     dFermiDist	= -diff(FermiDist); % Num�rica

% figure
%     plot(Energia,FermiDist,'b-',Energia,dFermiDist,'r-')
%     xlabel('Energia (meV)','FontSize',14)
%     title('Fermi distribution');
%     set(gcf,'Color',[1 1 1]);
%     set(gca,'FontSize',12);
% -----------------------------------


% C�LCULO DE LA DENSIDAD DE ESTADOS BCS NORMALIZADA PARA LA DISTRIBUCI�N DE
% GAPS
% ---------------------------
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
% --------------------------


% C�LCULO DE LA CONDUCTANCIA CONVOLUCIONANDO CON LA DISTRIBUCI�N DE FERMI
% --------------------------------
    Conductancia = conv(dFermiDist,DOS,'same');
% --------------------------------


% REPRESENTACI�N DEL RESULTADO
% ---------------------------------

    subplot(1,3,2,'Parent',Fig1);
    plot(Energia,DOS,'-','LineWidth',2,'Color',[0 0.4470 0.7410]);
    
    xlabel('Energy (meV)','FontSize',14);
    ylabel('DOS','FontSize',14);
    title('Normalized density of states','FontSize',14);
    
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontSize',12);
    
    ConductanciaNormalizada = Conductancia/200; % ��MAGIC NUMBER!!
       
    subplot(1,3,3,'Parent',Fig1);
    hold on
    plot(Curva(:,1),Curva(:,2),...
        'o','MarkerSize',10,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[0.8500 0.3250 0.0980]);
    plot(Energia,ConductanciaNormalizada,...
        '-','LineWidth',3,'Color',[0 0.4470 0.7410]);
    
    xlabel('energy (meV)','FontSize',18);
        xlim([-2 , 2]);
    ylabel('Normalized conductance','FontSize',16);
        ylim([0 2.5]);
    title('Normalized conductance','FontSize',16);
    
    set(gcf,'Color',[1 1 1]);
    set(gca,'FontSize',14);
    box on
    hold off
   
% GUARDADO DE DATOS
% ---------------------------
    dlmwrite('Fit.txt', [Energia,Conductancia],'delimiter','\t','newline','pc');
