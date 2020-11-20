function [state, P, var] = EKFTestUpdateAccl0(state, P, P_new, accl, i)

persistent KHP Hi K KHi
if(isempty(KHP) || i == 2)
    KHP  = zeros(length(P(:,1)));
    Hi   = zeros(3,3);
    K    = zeros(length(P(:,1)), 3);
    KHi  = zeros(length(P(:,1)));
end

P_psc = 65536;

accl_var = 1.0/9.81;
accl = accl/9.81;

q = state(1:4);

H =[ -2*q(3),  2*q(4), -2*q(1), 2*q(2), 0, 0;
      2*q(2),  2*q(1),  2*q(4), 2*q(3), 0, 0;
           0, -4*q(2), -4*q(3),      0, 0, 0];
if(mod(i,20) == 1)
    P_new = P_new - KHP;     
    P_new = (P_new + P_new')/2;

    for i = 1:length(state)
        if(P_new(i,i) < 0)
            P_new(i,i) = 0;
        end
    end
    
    KHP = zeros(length(P(:,1)));
    P = P_new;
    
    Hi = H * P * H' / P_psc + diag([accl_var,accl_var,accl_var]);
    Hi = diag([Hi(1,1),Hi(2,2),Hi(3,3)]);
    for id = 1 : 3
        K(:,id) = P * H(id,:)'/ P_psc / Hi(id,id);
    end
    
    %KHi = K*H*P;
end

err = accl - [2*q(2)*q(4) - 2*q(1)*q(3);
              2*q(1)*q(2) + 2*q(3)*q(4);
                1 - 2*(q(2)^2 + q(3)^2)];

Var = zeros(1,3);
for id = 1 : 3
    Var(id) = err(id)^2 / (Hi(id,id) * accl_var);
end
var = Var(3);

if(Var(1) < 10 && Var(2) < 10 && Var(3) < 10)
    state = state + K(:,1) * err(1);
    state = state + K(:,2) * err(2);
    state = state + K(:,3) * err(3);
    KHP = KHP + K*H*P;
end

if state(5) > 0.005
    state(5) = 0.005;
elseif state(5) < -0.005
    state(5) = -0.005;
end

if state(6) > 0.005
    state(6) = 0.005;
elseif state(6) < -0.005
    state(6) = -0.005;
end

