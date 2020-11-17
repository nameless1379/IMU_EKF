function [state_out, accl_out, gyro_out] = EKFTestInput(state, rotation, dt)

accl_var = 0.5;
gyro_var = 0.0001;

gyro_out = rotation + state(5:7) + normrnd(0,sqrt(gyro_var),3,1);
dq = [1;
    0.5*rotation(1)*dt;
    0.5*rotation(2)*dt;
    0.5*rotation(3)*dt;
    ];
q = state(1:4);
q_new = QuatMult(q,dq);
q_new = q_new/norm(q_new);

Tbn = Quat2Tbn(q_new);

state_out = [q_new; state(5:7)];
accl_out  = transpose(Tbn) * [0;0;9.81] + normrnd(0,sqrt(accl_var),3,1);

end

