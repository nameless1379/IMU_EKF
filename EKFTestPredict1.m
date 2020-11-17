function [state, P] = EKFTestPredict1(state, P, rotation, dt)

gyro_var      = 1e-4;
gyro_bias_var = 1e-9;

q = state(1:4);
dAngTruth = rotation - [state(5);state(6);state(7)];
deltaQuat = [1;0.5*dAngTruth(1)*dt;
               0.5*dAngTruth(2)*dt;
               0.5*dAngTruth(3)*dt];

quatNew = QuatMult(q,deltaQuat);
quatNew = quatNew/norm(quatNew);
state(1:4) = quatNew;

F = [               1, -0.5*dAngTruth(1)*dt, -0.5*dAngTruth(2)*dt, -0.5*dAngTruth(3)*dt,  (dt*q(2))/2,  (dt*q(3))/2,  (dt*q(4))/2;
  0.5*dAngTruth(1)*dt,                    1,  0.5*dAngTruth(3)*dt, -0.5*dAngTruth(2)*dt, -(dt*q(1))/2,  (dt*q(4))/2, -(dt*q(3))/2;
  0.5*dAngTruth(2)*dt, -0.5*dAngTruth(3)*dt,                    1,  0.5*dAngTruth(1)*dt, -(dt*q(4))/2, -(dt*q(1))/2,  (dt*q(2))/2;
  0.5*dAngTruth(3)*dt,  0.5*dAngTruth(2)*dt, -0.5*dAngTruth(1)*dt,                    1,  (dt*q(3))/2, -(dt*q(2))/2, -(dt*q(1))/2;
                    0,                    0,                    0,                    0,            1,            0,            0;
                    0,                    0,                    0,                    0,            0,            1,            0;
                    0,                    0,                    0,                    0,            0,            0,            1];
Q = [ gyro_var * dt^2*(1 - q(1)^2)/4, gyro_var * dt^2*-(q(1)*q(2))/4, gyro_var * dt^2*-(q(1)*q(3))/4, gyro_var * dt^2*-(q(1)*q(4))/4,  (gyro_bias_var*q(2)*dt)/2,  (gyro_bias_var*q(3)*dt)/2, (gyro_bias_var*q(4)*dt)/2;   
      gyro_var * dt^2*-(q(1)*q(2))/4, gyro_var * dt^2*(1 - q(2)^2)/4, gyro_var * dt^2*-(q(2)*q(3))/4, gyro_var * dt^2*-(q(2)*q(4))/4, -(gyro_bias_var*q(1)*dt)/2,  (gyro_bias_var*q(4)*dt)/2,-(gyro_bias_var*q(3)*dt)/2;
      gyro_var * dt^2*-(q(1)*q(3))/4, gyro_var * dt^2*-(q(2)*q(3))/4, gyro_var * dt^2*(1 - q(3)^2)/4, gyro_var * dt^2*-(q(3)*q(4))/4, -(gyro_bias_var*q(4)*dt)/2, -(gyro_bias_var*q(1)*dt)/2, (gyro_bias_var*q(2)*dt)/2;
      gyro_var * dt^2*-(q(1)*q(4))/4, gyro_var * dt^2*-(q(2)*q(4))/4, gyro_var * dt^2*-(q(3)*q(4))/4, gyro_var * dt^2*(1 - q(4)^2)/4,  (gyro_bias_var*q(3)*dt)/2, -(gyro_bias_var*q(2)*dt)/2,-(gyro_bias_var*q(1)*dt)/2;
      (gyro_bias_var*q(2)*dt)/2, -(gyro_bias_var*q(1)*dt)/2, -(gyro_bias_var*q(4)*dt)/2,  (gyro_bias_var*q(3)*dt)/2, gyro_bias_var,             0,              0;
      (gyro_bias_var*q(3)*dt)/2,  (gyro_bias_var*q(4)*dt)/2, -(gyro_bias_var*q(1)*dt)/2, -(gyro_bias_var*q(2)*dt)/2,             0, gyro_bias_var,              0;
      (gyro_bias_var*q(4)*dt)/2, -(gyro_bias_var*q(3)*dt)/2,  (gyro_bias_var*q(2)*dt)/2, -(gyro_bias_var*q(1)*dt)/2,             0,             0,  gyro_bias_var];

P = F*P*F' + Q;
