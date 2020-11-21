clear
dt = 0.001;
t_total = 1000;

test_len = t_total/dt;

real_state = zeros(7,test_len);
real_euler = zeros(3,test_len);
real_state(:,1) = [1;0;0;0;0.005;-0.002;0.002];

my_state0 = zeros(7,test_len);
my_state0(:,1) = [1;0;0;0;0;0;0];
my_var0(:,1) = zeros(1,test_len);

my_state1 = zeros(7,test_len);
my_state1(:,1) = [1;0;0;0;0;0;0];
my_var1(:,1) = zeros(1,test_len);

my_euler0 = real_euler;
my_euler1 = real_euler;
P0     = zeros(6,6);
P0_new = zeros(6,6);
P1     = zeros(7,7);

figure(1);
clf;

input_Ax = 0;
input_fx = 1;
input_Ay = 0;
input_fy = 2;

for i = 2: length(real_state)
    rotation = zeros(3,1);
    rotation(1) = input_Ax * cos(2*pi* input_fx * i/1000);
    rotation(2) = input_Ay * sin(2*pi* input_fy * i/1000);
    
    [real_state(:,i), accl_input, gyro_input] = EKFTestInput(real_state(:,i-1),rotation, dt);
    real_euler(:,i) = QuatToEul(real_state(1:4, i)) * 180 / pi; 
    
    [my_state0(1:6,i), P0_new] = EKFTestPredict0(my_state0(1:6,i-1), P0, gyro_input, dt, i);
    [my_state0(1:6,i), P0, my_var0(i)] = EKFTestUpdateAccl0(my_state0(1:6,i), P0, P0_new, accl_input, i);
    my_euler0(:,i) = QuatToEul(my_state0(1:4, i)) * 180 / pi; 
    
    %[my_state1(:,i), P1] = EKFTestPredict1(my_state1(:,i-1), P1, gyro_input, dt);
    %[my_state1(:,i), P1] = EKFTestUpdateAccl1(my_state1(:,i), P1, accl_input);
    %my_euler1(:,i) = QuatToEul(my_state1(1:4, i)) * 180 / pi; 
    
    if(mod(i , 5000) == 0)
        P0 / 65536
        %P1
    end
end

subplot (2,3,1);
t = linspace(0,t_total,test_len);
hold on
plot(t,real_euler(1,:));
plot(t,my_euler0(1,:));
plot(t,my_euler1(1,:));
subplot (2,3,2);
hold on
plot(t,real_euler(2,:));
plot(t,my_euler0(2,:));
plot(t,my_euler1(2,:));
subplot (2,3,3);
hold on
plot(t,real_euler(3,:));
plot(t,my_euler0(3,:));
plot(t,my_euler1(3,:));
subplot (2,3,4);
hold on
plot(t,my_state0(5,:));
plot(t,my_state1(5,:));
subplot (2,3,5);
hold on
plot(t,my_state0(6,:));
plot(t,my_state1(6,:));
subplot (2,3,6);
hold on
plot(t,my_var0);
plot(t,my_var1);