load('StateAndCovariancePrediction.mat');
F = subs(F, dt, 2);
F = subs(F, [dax_b, day_b, daz_b], [0,0,0]);

PP = F*P*transpose(F) - P;
PP = PP(1:6,1:6);

for i = 2:length(PP(1,:))
    for j = 1 : i-1
        PP(i,j) = 0;
    end
end

PP = subs(PP,[OP_l_7_c_1_r_,OP_l_7_c_2_r_,OP_l_7_c_3_r_,OP_l_7_c_4_r_,OP_l_7_c_5_r_,OP_l_7_c_6_r_,OP_l_7_c_7_r_,OP_l_6_c_5_r_],[0,0,0,0,0,0,0,0]);
[PP,SPP]=OptimiseAlgebra(PP,'SPP');