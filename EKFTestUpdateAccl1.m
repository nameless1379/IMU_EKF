function [state, P] = EKFTestUpdateAccl1(state, P, accl)

accl_var = 1.0/9.81;
accl = accl/9.81;

q = state(1:4);

H =[ -2*q(3),  2*q(4), -2*q(1), 2*q(2), 0, 0, 0;
      2*q(2),  2*q(1),  2*q(4), 2*q(3), 0, 0, 0;
           0, -4*q(2), -4*q(3),      0, 0, 0, 0];
  
Hi = H * P * H' + diag([accl_var,accl_var,accl_var]);
err = accl - [2*q(2)*q(4) - 2*q(1)*q(3);
              2*q(1)*q(2) + 2*q(3)*q(4);
                1 - 2*(q(2)^2 + q(3)^2)];
var = (err'*inv(Hi)*err)/accl_var;
if(var > 25)
    state = state;
    P = P;
    return;
end
  
K = (P * H')/Hi;

state = state + K * err;
P = P - K*H*P;     
P = (P + P')/2;

for i = 1:length(state)
    if(P(i,i) < 0)
        P(i,i) = 0;
    end
end