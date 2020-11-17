function [state, P] = EKFTestUpdateAccl0(state, P, P_new, accl, i)

persistent IKH
if(isempty(IKH) || i == 2)
    IKH = eye(length(P(:,1)));
end

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
    %fprintf('outlier rejected with SD\r\n')
    var
else
    K = P * H'/Hi;
    state = state + K * err;
    IKH = (eye(length(P(:,1))) -  K*H) * IKH;
end

if(mod(i,10) == 0)
    P_new = IKH*P_new;     
    P_new = (P_new + P_new')/2;

    for i = 1:length(state)
        if(P_new(i,i) < 0)
            P_new(i,i) = 0;
        end
    end
    
    IKH = eye(length(P(:,1)));
    P = P_new;
else
    P = P;
end
