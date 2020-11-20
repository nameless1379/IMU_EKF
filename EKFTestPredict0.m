function [state, P_new] = EKFTestPredict0(state, P, rotation, dt, i)

P_psc = 65536;
persistent avg_d_angle dT
if(isempty(avg_d_angle) || i == 2)
    avg_d_angle = zeros(3,1);
    dT = 0;
end

gyro_var      = 1e-7;
gyro_bias_var = 1e-5;

q = state(1:4);
dAngTruth   = rotation - [state(5);state(6);0];

avg_d_angle = avg_d_angle + (dAngTruth*dt);
dT = dT + dt;

deltaQuat = [1;0.5*dAngTruth(1)*dt;0.5*dAngTruth(2)*dt;0.5*dAngTruth(3)*dt;];

quatNew = QuatMult(q,deltaQuat);
state(1:4) = quatNew;

if(mod(i,20) == 1)
    state(1:4) = state(1:4)/norm(state(1:4));
    F = [               1, -0.5*avg_d_angle(1), -0.5*avg_d_angle(2), -0.5*avg_d_angle(3),  (dT*q(2))/2,  (dT*q(3))/2;
       0.5*avg_d_angle(1),                   1,  0.5*avg_d_angle(3), -0.5*avg_d_angle(2), -(dT*q(1))/2,  (dT*q(4))/2;
       0.5*avg_d_angle(2), -0.5*avg_d_angle(3),                   1,  0.5*avg_d_angle(1), -(dT*q(4))/2, -(dT*q(1))/2;
       0.5*avg_d_angle(3),  0.5*avg_d_angle(2), -0.5*avg_d_angle(1),                   1,  (dT*q(3))/2, -(dT*q(2))/2;
                        0,                    0,                   0,                  0,            1,            0;
                        0,                    0,                   0,                  0,            0,            1];
    Q = [ gyro_var * dT^2*(1 - q(1)^2)/4, gyro_var * dT^2*-(q(1)*q(2))/4, gyro_var * dT^2*-(q(1)*q(3))/4, gyro_var * dT^2*-(q(1)*q(4))/4,  (gyro_bias_var*q(2)*dT^3)/2,  (gyro_bias_var*q(3)*dT^3)/2;
          gyro_var * dT^2*-(q(1)*q(2))/4, gyro_var * dT^2*(1 - q(2)^2)/4, gyro_var * dT^2*-(q(2)*q(3))/4, gyro_var * dT^2*-(q(2)*q(4))/4, -(gyro_bias_var*q(1)*dT^3)/2,  (gyro_bias_var*q(4)*dT^3)/2;
          gyro_var * dT^2*-(q(1)*q(3))/4, gyro_var * dT^2*-(q(2)*q(3))/4, gyro_var * dT^2*(1 - q(3)^2)/4, gyro_var * dT^2*-(q(3)*q(4))/4, -(gyro_bias_var*q(4)*dT^3)/2, -(gyro_bias_var*q(1)*dT^3)/2;
          gyro_var * dT^2*-(q(1)*q(4))/4, gyro_var * dT^2*-(q(2)*q(4))/4, gyro_var * dT^2*-(q(3)*q(4))/4, gyro_var * dT^2*(1 - q(4)^2)/4,  (gyro_bias_var*q(3)*dT^3)/2, -(gyro_bias_var*q(2)*dT^3)/2;
             (gyro_bias_var*q(2)*dT^3)/2,   -(gyro_bias_var*q(1)*dT^3)/2,   -(gyro_bias_var*q(4)*dT^3)/2,    (gyro_bias_var*q(3)*dT^3)/2,          gyro_bias_var* dT^2,                            0;
             (gyro_bias_var*q(3)*dT^3)/2,    (gyro_bias_var*q(4)*dT^3)/2,   -(gyro_bias_var*q(1)*dT^3)/2,   -(gyro_bias_var*q(2)*dT^3)/2,                            0,          gyro_bias_var* dT^2];

    P_new = F*P*F' + Q * P_psc;
    avg_d_angle = zeros(3,1);
    dT = 0;
else
    P_new = P;
end
