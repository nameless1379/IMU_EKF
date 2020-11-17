function [state, P, var] = EKFTestUpdateAccl0(state, P, P_new, accl, i)

persistent KH
if(isempty(KH) || i == 2)
    KH = zeros(length(P(:,1)));
end

accl_var = 1.0/9.81;
accl = accl/9.81;

q = state(1:4);

H =[ -2*q(3),  2*q(4), -2*q(1), 2*q(2), 0, 0;
      2*q(2),  2*q(1),  2*q(4), 2*q(3), 0, 0;
           0, -4*q(2), -4*q(3),      0, 0, 0];
Hi = H * P * H' + diag([accl_var,accl_var,accl_var]);
Hi = diag([Hi(1,1),Hi(2,2),Hi(3,3)]);

err = accl - [2*q(2)*q(4) - 2*q(1)*q(3);
              2*q(1)*q(2) + 2*q(3)*q(4);
                1 - 2*(q(2)^2 + q(3)^2)];
var = (err'*inv(Hi)*err)/accl_var;

if(var > 10)
    %fprintf('outlier rejected with SD\r\n')
    %var
else
    K = P * H'/Hi;
    state = state + K * err;
    KH = KH + K*H;
end

if(mod(i,20) == 0)
    P_new = P_new - KH * P_new;     
    P_new = (P_new + P_new')/2;

    for i = 1:length(state)
        if(P_new(i,i) < 0)
            P_new(i,i) = 0;
        end
    end
    
    KH = zeros(length(P(:,1)));
    P = P_new;
else
    P = P;
end
