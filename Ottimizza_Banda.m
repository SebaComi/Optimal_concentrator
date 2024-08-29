clc
% clear
close all
load('LIMITI_per_il_gradiente_proiettato.mat')
% Definizione di dati
DATA_global.tipo_mesh = 7;
DATA_global.tipo_base =      'hex';  % 'corone circolari';'rettangoli';
DATA_global.raggio_ext = 0.6;  % [m]
DATA_global.raggio_osta = 0.02; % [m]
DATA_global.raggio_fuoco = 0.01; % [m]
DATA_global.area_fuoco = DATA_global.raggio_fuoco^2*pi;
DATA_global.lato_hex = 0.005; % [m]
DATA_global = Dati_cloak(DATA_global);
n_angoli = 6;
n_freq = 6;
delta_angoli = 120/n_angoli;
delta_freq = 2/n_freq;

centri_banda = [10, 10, 25, 25, 40, 40];
lambda = sqrt(2.2e9/1000) ./centri_banda/1000;
% raggi_mant = [0.075, 0.14, 0.075, 0.12, 0.055, 0.095]; % [m]
raggi_mant = lambda.*[0.5  1  1.25  2  1.5  2.5]; % [m]
raggi_mant(2) = 0.140;
tic
parfor n_simu = 1:6
% for n_simu = 1
tic
    
centro_banda = centri_banda(n_simu);
[freq, angoli] = meshgrid( ...
                    (centro_banda -1 +delta_freq/2):delta_freq:(centro_banda + 1 - delta_freq/2), ...
                    (-60+delta_angoli/2):delta_angoli:(60-delta_angoli/2) );

DATA = DATA_global;
DATA.omega = freq(:) *1000*2*pi;
DATA.k_amb = DATA.omega / DATA.c_0;
DATA.direz = angoli(:) * pi/180;
DATA.n_frq = length(DATA.omega);
DATA.hh_mesh = 2*pi / max(DATA.k_amb) / 5;
DATA.mesh = 'mesh_cloak';
DATA.raggio_mant = raggi_mant(n_simu); % [m]

% Creo MESH e FEM
[MESH, FE_SPACE] = crea_MESH_e_FEM(DATA, true);

%% Definizione basi
DATA = def_basi(DATA,MESH);

% Info
fileID = fopen(['Banda' num2str(centro_banda) '/log - ' char(date) ' - ' num2str(DATA.raggio_mant,'%.3f')   '.txt'],'w');
fprintf(fileID,'\n **** PROBLEM''S SIZE INFO ****\n');
fprintf(fileID,' * Numero di vertici       = %d \n',MESH.numVertices);
fprintf(fileID,' * Numero di nodi          = %d \n',MESH.numNodes);
fprintf(fileID,' * Numero di elementi      = %d \n',MESH.numElem);
fprintf(fileID,' * Numero di contorni      = %d \n',MESH.n_contorni);
fprintf(fileID,' * Numero di domini        = %d \n',MESH.n_domini);
fprintf(fileID,' * Numero di basi di controllo  = %d \n\n',DATA.n_basi);


%% Generazione matrici dominio
quali = {'A_0','A_j', 'B_0','B_j','C_e','D  ','E  ','B_b','F_b'};
        [ A_0,  A_j,   B_0,  B_j,  C_e,  D,    E,    B_b,  F_b] = Genera_matrici_parfor(DATA, MESH, FE_SPACE,quali);

%% Inizializziamo
fprintf('Il ciclo %d inizia a risolvere',n_simu)
% Y_inc nei vertici
Y_inc = zeros(MESH.numNodes,DATA.n_frq);
for ff = 1:DATA.n_frq
    Y_inc(:,ff) = DATA.y_inc_fun({MESH.nodes(1,:),MESH.nodes(2,:)},DATA.k_amb(ff),DATA.direz(ff)).';
end

% Inizializzo i vettori Y, P e U
Y = zeros(MESH.numNodes,DATA.n_frq);
P = Y;
U = zeros(DATA.n_basi,1);
U_old = U;
V = U;
V_old = V;

dJ_U = U;
dJ_V = dJ_U; 

%% Cose per il costo
[A_costo,N_costo] = calcola_vicini(DATA);

%% Iteriamo
% Costanti di iterazione
tau =  30;
MaxIt =  min(2000, DATA.n_basi*5);
lmb_u = 1e-6;  % costo di U -> Bulk
lmb_v = 1e-6;   % costo di V -> rho

tol = 1e-3;
uscita = false;

ii = 1; if tau == 0,  MaxIt = 2;     end
J = zeros(MaxIt,1); J_control = J;
J(1) = 1e10;
J_parz = zeros(DATA.n_frq,1);
tempo_medio = 0; contatore = 0; fattore = 1.2;
%%
while ii < MaxIt && ~uscita
tic
dJ_U(:) = 0;
dJ_V(:) = 0;
%% Definiamo Bu e Av     
    Av = ttv(A_j,exp(-V)-1,3);
    Av = (spmatrix(Av) + A_0) / DATA.rho_0;

    Bu = ttv(B_j,exp(-U)-1,3);
    Bu = (spmatrix(Bu) + B_0) / DATA.B_0;

    for ff = 1:DATA.n_frq
        omega = DATA.omega(ff);
        omega2 = omega^2;
        k_0 = DATA.k_amb(ff); 
        Atot = Av - Bu * omega2 + (1i*k_0 + 1/DATA.raggio_ext/2)/DATA.rho_0 * C_e;
        Fy = D{ff} * ( (exp(-U) -1) /DATA.B_0 * omega2 ) - E{ff} * ((exp(-V) -1) /DATA.rho_0);
        Y(:,ff) = Atot \ Fy;
        
        J_parz(ff) = (Y_inc(:,ff) + Y(:,ff))' * B_b * (Y_inc(:,ff) + Y(:,ff) ) / DATA.area_fuoco;
    end

    for ff = 1:DATA.n_frq
        % Definiamo Ap (= Atot)
        omega = DATA.omega(ff);
        omega2 = omega^2;
        k_0 = DATA.k_amb(ff);
        % Definiamo Ay             
        Atot = Av - Bu * omega2 + (1i*k_0 + 1/DATA.raggio_ext/2)/DATA.rho_0 * C_e;
        % Definiamo Fp
        coef = prod([J_parz(1:DATA.n_frq ~= ff); 1])^(1/DATA.n_frq)*J_parz(ff)^(1/DATA.n_frq-1);   % *1 nel caso in cui ci sia solo una frq 
        Fp = - coef * (B_b * Y(:,ff) + F_b{ff} );

        % Solve for P          
        P(:,ff) = conj(Atot)  \  Fp;
        
        %% Calcolo del costo e dei gradienti per la freq ff  
        dJ_loc = real( P(:,ff)' *   ( E{ff} + spmatrix(ttv(A_j, Y(:,ff),2)) ));
        dJ_loc = dJ_loc(:) ./ exp(V) /DATA.rho_0 ;
        dJ_V = dJ_V + dJ_loc;

        dJ_loc = real( P(:,ff)' *   ( D{ff} + spmatrix(ttv(B_j, Y(:,ff), 2) ) ));
        dJ_loc = -omega2/DATA.B_0 * dJ_loc(:) ./ exp(U);
        dJ_U = dJ_U + dJ_loc;
    end
J(ii+1) = -prod(J_parz)^(1/DATA.n_frq);

% Costo sul gradiente del controllo
J_control(ii+1) = lmb_u*2 *transpose(U) *(diag(N_costo) - A_costo) *U   +   lmb_v*2 *transpose(V) *(diag(N_costo)-A_costo) *V;
J(ii+1) = J(ii+1) + J_control(ii+1);
dJ_V = dJ_V + lmb_v*4*(diag(N_costo) - A_costo)*V;
dJ_U = dJ_U + lmb_u*4*(diag(N_costo) - A_costo)*U;

if ii == 1
    dJ0 = norm(dJ_U);
end

%% Stopping criteria    
if real(J(ii+1)) > real(J(ii))%  ||  ~DATA.vincolo(U)
    if tau > 1e-5%/2^DATA.n_frq
        tau = tau/1.6;%/1.1^DATA.n_frq;
        U = U_old;
        V = V_old;
        J(ii+1) = 0;
        ritento = true;
        contatore = 0;
        fattore = fattore^0.9;
    else
        fprintf(fileID,repmat('\b',1,9));        fprintf(fileID,'                Concluso in %.3f s.\n',toc);
        fprintf(fileID,'Ottimizzazione terminata dopo %i step: \n',ii);
        fprintf(fileID,'    anche con step piccoli, il costo non diminuisce,   tau = %5.4f \n',tau);
        uscita = 1;
    end
else
    U_old = U;
    V_old = V;
    U = U - tau * dJ_U;% / real(prod(J_parz)^(1/DATA.n_frq));
    V = V - tau * dJ_V;% / real(prod(J_parz)^(1/DATA.n_frq));
    [U,V] = proiettato(U,V,rho_lim,bulk_lim);    
    if norm([U,V]-[U_old,V_old]) > tol*dJ0 - eps %norm(dJ_U) > tol*dJ0 - eps
        ritento = false;
        ii = ii+1;
        if contatore == 10
            tau = tau * fattore;
            contatore = 0;
        else
            contatore = contatore + 1;
        end

    else
        U = U_old;
        V = V_old;
        fprintf(fileID,repmat('\b',1,9));        fprintf(fileID,'                Concluso in %.3f s.\n',toc);
        fprintf(fileID,'\nOttimizzazione terminata con successo dopo %i step. \n',ii);
        fprintf(fileID,'    Proseguendo, il costo diminuirebbe troppo poco. \n            |dJ|/|dJ0| = %3.2f*1e-4 < tol\n',1e4*norm(dJ_U)/dJ0);
        uscita = 1;
    end
end
%% Disegnamo          
if uscita || ii == MaxIt
    
    fig = figure('visible','off');
    
    for ff = 1:min([DATA.n_frq 3])
        subplot(2,3,ff)
        disegna_risultato(MESH,abs(Y(1:MESH.numVertices,ff)+Y_inc(1:MESH.numVertices,ff)))
%         title('Norma')
        title(['f = ' num2str(DATA.omega(ff)/2/pi/1000,3) ' kHz,' ...
                'Dir = (' num2str(cosdir(DATA.direz(ff))',3) ')'])
    end
    
    subplot(235)
    cla
    disegna_u(exp(U),DATA,round(10/DATA.raggio_mant*0.04))
    title('Controllo \kappa')
    xlabel('x'); ylabel('y'); view(2)
    colorbar

    subplot(236)
    cla
    disegna_u(exp(V),DATA,round(10/DATA.raggio_mant*0.04))
    title('Controllo \rho')
    xlabel('x'); ylabel('y'); view(2)
    colorbar
    
    colormap('jet')
    
    subplot(234)
    cla
    plot(exp(V),exp(U),'x')
    xlabel('\rho');ylabel('\kappa');
    drawnow
    name = ['Banda' num2str(centri_banda(n_simu)) '/Optimization - ' char(date) ' - ' num2str(DATA.raggio_mant,'%.3f')  '.jpg'];
    saveas(fig,name,'jpg');
    
end

%% Informazioni ciclo
tempo_medio = (tempo_medio*(ii-2) + toc)/(ii-1);
if ~uscita
%     fprintf(repmat('\b',1,testo))
    fprintf(fileID, 'Ciclo %d                Concluso in %.3f s.\n Media: %.2f  J = %.2f  tau = %.3f  T_tot: %.1f s\n',ii,toc,tempo_medio,J(ii),tau,tempo_medio*ii);
    for ff = 1:DATA.n_frq
        limitante = sum(MESH.vertices.^2) <= (DATA.raggio_mant/2)^2;
        limitante = find(limitante);
        [valore,dove] = max(abs(Y(limitante,ff)+Y_inc(limitante,ff)));
        fprintf(fileID, 'Massimo in (%.3f; %.3f) p = %.2f', ...
                             MESH.vertices(1,limitante(dove)),MESH.vertices(2,limitante(dove)), valore);
        fprintf(fileID, '\t J_parz = %3f \n', J_parz(ff));
    end
    if ii < MaxIt
        fprintf(fileID, 'Ciclo %d in corso',ii+1);
    end
end
end % Fine del while

if ii == MaxIt
    fprintf(fileID,'Ottimizzazione terminata dopo %i step\n',ii);
    fprintf(fileID,'    Il numero massimo di iterazioni è stato raggiunto\n');
end 

% Disp costo
J = J(3:min(ii+1,MaxIt));
J_control = J_control(3:min(ii+1,MaxIt));
if all(abs(imag(J)) < abs(real(J)) * 1e-6)
    J = real(J);
    J_control = real(J_control);
else
    warning('J è complesso con raggio lente uguale a %.4f',raggi_mant(n_simu));
end

fig1 = figure('visible','off');
plot(J)
title('Costo VS iterazioni')
saveas(fig1,['Banda' num2str(centri_banda(n_simu)) '/Costo - ' char(date) ' - ' num2str(DATA.raggio_mant,'%.3f') '.jpg'],'jpg')

fclose(fileID);

%% Salvataggi
tempo_tot = tempo_medio*ii;
m = matfile(['../Analisi Banda/Banda' num2str(centri_banda(n_simu)) '/Data - ' char(date) ' - ' num2str(DATA.raggio_mant,'%.3f')   '.mat'],'writable',true);
m.Y=Y;
m.U=U;
m.V=V;
m.J=J;
m.MESH=MESH;
m.DATA=DATA;
m.tempo_tot=tempo_tot;
m.J_control=J_control;

fprintf('Finito il ciclo numero %d in %.2f ore\n',n_simu,toc/3600)

end
toc


return

%% Postprocessing
% centro_banda = [10, 25, 40];
% raggi_mant = 0.04:0.005:0.15; % [m]
% legenda = {};
% figure(170); clf; figure(171); clf
% for ff = centro_banda
%     Costo_vs_raggio = zeros(length(raggi_mant),1);
%     for rr = 1:length(raggi_mant)
%         namefile = ['24.03.2023\Banda' num2str(ff) '\Data - 24-Mar-2023 - ' num2str(raggi_mant(rr)) '.mat']; 
%         load(namefile)
%         Costo_vs_raggio(rr) = J(end);
%         
%         Y_inc = DATA.y_inc_fun({MESH.vertices(1,:)',MESH.vertices(2,:)'},DATA.k_amb,DATA.direz);
%         
%         
%         limitante = sum(MESH.vertices.^2) <= (DATA.raggio_fuoco)^2;
%         limitante = find(limitante);
%         [valore,dove] = max(abs(Y(limitante)+Y_inc(limitante)));
% %         fprintf('Massimo in (%.3f; %.3f) p = %.2f \n', MESH.vertices(1,limitante(dove)),MESH.vertices(2,limitante(dove)), valore);
%         guadagni(rr,ff) = 10*log10(valore);
%         posizioni(:,rr,ff) = MESH.vertices(:,limitante(dove));
%     end
%     figure(170);
%     hold on
%     plot(raggi_mant/(2*pi/DATA.k_amb), Costo_vs_raggio);
%     title(['J Banda' num2str(ff)])
%     grid on
%     ylim([-1000 0])
%     xlabel('r/\lambda')
%     
%     figure(171);
%     hold on
%     plot(raggi_mant/(2*pi/DATA.k_amb), guadagni(:,ff));
%     title(['Gain Banda' num2str(ff)])
%     ylim([0 16])
%     grid on
%     xlabel('r/\lambda')
%     
%     legenda = [legenda, {[num2str(ff), ' kHz']}];
% end
% figure(170)
% legend(legenda{:})
% figure(171)
% legend(legenda{:})


%% Nested functions
%  disegno_mesh(MESH)
function [] = disegno_mesh(MESH)
    % Disegno mesh
    fprintf('Disegno mesh. '); tic
    % Visualizziamo la mesh, i contorni e le normali
    % figure(2)
    pdeplot(MESH.vertices,MESH.elements(1:3,:));
    hold on
%     colori = get(gca,'ColorOrder'); colori = [colori;colori];
    % colori(3,:) = [];
        % Domini e normali
%     for dd = 1:MESH.n_domini
%         % Elementi
%         h = pdeplot(MESH.vertices,MESH.elements(1:3,MESH.nodi_elements(:,dd)));
%         h.Color = colori(dd,:);
% 
%         % Normali 
%         cc_del_dd = unique(MESH.boundaries(5,MESH.nodi_boundaries(:,dd)));
%         for cc = cc_del_dd(:)'
%            quiver(MESH.vertices(1,MESH.nodi_contorni{cc}), MESH.vertices(2,MESH.nodi_contorni{cc}),...
%                 MESH.normali{dd,cc}(2,:) /32, MESH.normali{dd,cc}(3,:) /32, ...
%                 'Color',colori(dd,:), 'LineWidth',1,'AutoScale','off')
%         end
%     end
        % Contorni
    for cc = 1:MESH.n_contorni
        linea = MESH.boundaries(:,MESH.boundaries(5,:) == cc);
        [~,index] = sort(linea(4,:));   % aggiornamento 14/07/21:
                                        % guardo l'ascissa curvilinea sulla riga 4 e non 3 perché redbKit
                                        % sovrascrive la riga 3...
        linea = linea(:,index);
        plot(MESH.vertices(1,[linea(1,:) linea(2,end)]),  MESH.vertices(2,[linea(1,:) linea(2,end)]), ...
            'k','LineWidth',.8)
    end

    drawnow
    fprintf('             Fatto in %.3f s\n',toc)
end

