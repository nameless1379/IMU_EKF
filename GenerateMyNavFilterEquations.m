clear all;
reset(symengine);

syms q0 q1 q2 q3 real % quaternions defining attitude of body axes relative to local NED
syms dax day daz real % IMU delta angle measurements in body axes - rad
syms dax_b day_b daz_b real % delta angle bias - rad
syms daVar dvVar dabVar real; % IMU delta angle and delta velocity measurement variances
syms dt real % IMU time step - sec
syms gravity real % gravity  - m/sec^2

dAngMeas = [  dax;   day;   daz];
dAngBias = [dax_b; day_b; daz_b];

% define the quaternion rotation vector for the state estimate
quat = [q0;q1;q2;q3];
% derive the truth body to nav direction cosine matrix
Tbn = Quat2Tbn(quat);

dAngTruth = dAngMeas - dAngBias;
% define the attitude update equations
% use a first order expansion of rotation to calculate the quaternion increment
% acceptable for propagation of covariances
deltaQuat = [1;
    0.5*dAngTruth(1)*dt;
    0.5*dAngTruth(2)*dt;
    0.5*dAngTruth(3)*dt;
    ];
quatNew = QuatMult(quat,deltaQuat);
% define the IMU error update equations
dAngBiasNew = dAngBias;

stateVector    = [quat;dAngBias];
newStateVector = [quatNew;dAngBiasNew];
nStates=numel(stateVector);

% derive the state transition matrix
F = jacobian(newStateVector, stateVector);
% set the rotation error states to zero
%[F,SF]=OptimiseAlgebra(F,'SF');
% define a symbolic covariance matrix using strings to represent 
% '_l_' to represent '( '
% '_c_' to represent ,
% '_r_' to represent ')' 
% these can be substituted later to create executable code

%diagonal matrix
for i = 1:nStates
    eval(['syms OP_l_',num2str(i),'_c_',num2str(i), '_r_ real']);
    eval(['P(',num2str(i),',',num2str(i), ') = OP_l_',num2str(i),'_c_',num2str(i),'_r_;']);
end
    
for rowIndex = 2:nStates
    for colIndex = 1:rowIndex-1
        eval(['syms OP_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['P(',num2str(rowIndex),',',num2str(colIndex), ') = OP_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
        eval(['P(',num2str(colIndex),',',num2str(rowIndex), ') = OP_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end

save 'StatePrediction.mat';

%% derive the covariance prediction equations
% This reduces the number of floating point operations by a factor of 6 or
% more compared to using the standard matrix operations in code

% Error growth in the inertial solution is assumed to be driven by 'noise' in the delta angles and
% velocities, after bias effects have been removed. 

% derive the control(disturbance) influence matrix from IMu noise to state
% noise
G = jacobian(newStateVector, [dAngMeas;dAngBias]);
distMatrix = diag([daVar daVar daVar dabVar dabVar dabVar]);
Q = G*distMatrix*transpose(G);

PP = F*P*transpose(F) - P;

for rowIndex = 2:nStates
    for colIndex = 1:rowIndex-1
        PP(rowIndex, colIndex) = 0;
    end
end

%[PP,SPP]=OptimiseAlgebra(PP,'SPP');
save('StateAndCovariancePrediction.mat');
clear all;
reset(symengine);

%% acceleration measurement
load('StatePrediction.mat');
syms dvx dvy dvz real % IMU delta velocity measurements in body axes - m/sec
dVelMeas = [  dvx;   dvy;   dvz];

acclMeas = transpose(Tbn)*[0;0;1];

H_ACC = jacobian(acclMeas,stateVector); % measurement Jacobian
H_ACC_INV = H_ACC*P*transpose(H_ACC);

for rowIndex = 1:nStates
    for colIndex = 1:length(acclMeas)
        eval(['syms OK_l_',num2str(rowIndex),'_c_',num2str(colIndex), '_r_ real']);
        eval(['K(',num2str(rowIndex),',',num2str(colIndex), ') = OK_l_',num2str(rowIndex),'_c_',num2str(colIndex),'_r_;']);
    end
end

save('Acclerometer.mat')
clear all;
reset(symengine);