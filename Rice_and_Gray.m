function D12calc_P2 = Rice_and_Gray(T, Density, M, Tc, Vc, k12, B12)

% Calculates the tracer diffusivity of a given solute in a solvent.
% T - Absolute temperature in K in column format.
% Density - Density of the solvent in g/cucm  in column format.
% M - Molecular mass of the compounds by the following order solute, solvent.
% Tc - Critical temperature in K of the compounds by the following order solute, solvent.
% Vc - Critical volume in cucm/gmol of the compounds by the following order solute, solvent.
% k12 - Binary binary interaction parameter introduced in the LJ diameter combination rule.
% B12 - Binary binary interaction parameter



% Check input consistency
Tsize = size(T);
Dsensitysize = size(Density);


if Tsize(2)~=1 || Dsensitysize(2)~=1
    error('T and density should be column')
elseif Tsize(1)~=Dsensitysize(1)
    error('T and density column size not consistent')
end

numM = numel(M);
numTc = numel(Tc);
numVc = numel(Vc);


if numM~=numTc || numM~=numVc || numTc~=numVc
    error('Inconsistent size of M, Tc and Vc')
end



%% Constantes
Na = 6.0221367e23; %mol-1
kb = 1.380658e-16; % g cm2 s-2 K

%% calculos

solute.Vc = Vc(1);
solute.Tc = Tc(1);
solute.M = M(1);
solvent.Vc = Vc(2);
solvent.Tc = Tc(2);
solvent.M = M(2);


%% Secção comum
solvent.dlj = 0.7889e-8*solvent.Vc^(1/3);
solvent.elj = solvent.Tc/1.2593;
solute.dlj = 0.7889e-8*solute.Vc^(1/3);
solute.elj = solute.Tc/1.2593;

solute.m = solute.M / Na; %% confirmar
solvent.m = solvent.M / Na; %% confirmar

elj_12 = sqrt(solvent.elj*solute.elj);

T_1_est = T./solvent.elj;
T_2_est = T./solute.elj;
T_12_est = T./elj_12;

dlj_1_eff = 1.1532 * solvent.dlj * (1+ (1.8975*T_1_est).^0.5).^(-1/6);
dlj_2_eff = 1.1532 * solute.dlj * (1+(1.8975*T_2_est).^0.5).^(-1/6);

densd1 = Density.*Na/solvent.M;
density_est = densd1.*dlj_1_eff.^3; %% como está na Ana

compress_1 = pi/6*density_est;

m12 = solute.m*solvent.m/(solute.m+solvent.m);

g_sig_est_12 = 1./(1-compress_1).^3.*(1-compress_1+2*compress_1./(1+dlj_1_eff./dlj_2_eff)).*(1-compress_1+compress_1./(1+dlj_1_eff./dlj_2_eff));

a = -1.676382*density_est + 1.638561;
b = -8.516830*density_est + 8.631536;
c = -1.320347*density_est + 1.351067;
d = -5.062546*density_est + 5.409662;

F11 = 1 + 0.94605*density_est.^1.5 + 1.4022*density_est.^3 - 5.6898*density_est.^5 + 2.6626*density_est.^7;
F12 = (F11 + density_est.^1.7.* (a.*log(dlj_2_eff./dlj_1_eff) + b.*log(dlj_2_eff./dlj_1_eff).^2 + c.*log(solute.m/solvent.m) ) )./(1 + density_est.^3.*(d.*log(dlj_2_eff./dlj_1_eff)).^2);

for n=1:numel(F11)
    if F11(n) < 0
        error('Out of bounds, F11 < 0')
    end
    if F12(n) < 0
        error('Out of bounds, F12 < 0')
    end
end

dlj_12_P2 = (1-k12)*(solvent.dlj + solute.dlj)/2;
dlj_12_eff_P2 = 1.1532*dlj_12_P2.*(1+(1.8975*T_12_est).^0.5).^(-1/6);

    
D12calc1_HS = kb*T;
D12calc2_P2 = 8/3*densd1.*dlj_12_eff_P2.^2;
D12calc3_P2 = sqrt(2*pi*kb*T.*m12);
D12calc4_P2 = g_sig_est_12./F12+(B12)./(T_12_est.^1.5);
D12calc_P2 = D12calc1_HS./(D12calc2_P2.*D12calc3_P2.*D12calc4_P2);
