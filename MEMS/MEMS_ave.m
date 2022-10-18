clc;
clear;

g0=9.78049;%标准地球重力加速度(m/s^2)
rad_deg=0.01745329;%度转换成弧度的比例因子pi/180
wie=7.27220417e-5;%地球自转角速度(rad/s)
Re=6378393.0;%地球半径(m)
e = 1/298.257;  %地球扁率

pitch_error=0.0001*rad_deg;%初始平台水平失准角（粗对准以后的结果）
roll_error=0.0001*rad_deg;%初始平台水平失准角（粗对准以后的结果）
yaw_error=0.0005*rad_deg;%初始平台方位失准角（粗对准以后的结果）

TT=[1          yaw_error     -roll_error;
    -yaw_error 1             pitch_error;
    roll_error -pitch_error  1];

data = readmatrix('D:\桌面\研究生\导航学习\MEMS\LOG_20220929105433.csv');

Gx	= data(:,7);
Gy 	= data(:,8);
Gz 	= data(:,9);

Ax 	= data(:,18);
Ay 	= data(:,19);
Az 	= data(:,20);
%%%%%%%%%%%%%%%%%%%%%%%%%

pitch 	= data(:,3)*rad_deg;  %弧度
roll 	= data(:,4)*rad_deg;
yaw 	= data(:,5)*rad_deg;

Ve	= data(:,15);
Vn	= data(:,16);

longi	= data(:,10);  
lati = data(:,11);

N=length(yaw);

yaw0=yaw/rad_deg;
pitch0=pitch/rad_deg;
roll0=roll/rad_deg;

Vn0=Vn;
Ve0=Ve;

lati0=lati;
longi0=longi;

%%%%%%%%%%%%%%%%%%%%%%%%%
fai(1)=yaw(1)+yaw_error;   %弧度//fai为矩阵，对第一个数值进行赋值
theta(1)=pitch(1)+pitch_error;
gama(1)=roll(1)+roll_error;

Pitch(1)=pitch0(1);
Roll(1)=roll0(1);
Yaw(1)=yaw0(1);

v1(1)=Ve0(1);
v2(1)=Vn0(1);

phi(1)=lati0(1);
lamda(1)=longi0(1);


Tbn(1,1) = cos(roll(1))*cos(yaw(1)) - sin(roll(1))*sin(pitch(1))*sin(yaw(1));
Tbn(1,2)=-cos(pitch(1))*sin(yaw(1));
Tbn(1,3) = sin(roll(1))*cos(yaw(1)) + cos(roll(1))*sin(pitch(1))*sin(yaw(1));

Tbn(2,1)= cos(roll(1))*sin(yaw(1)) + sin(roll(1))*sin(pitch(1))*cos(yaw(1));
Tbn(2,2) = cos(pitch(1))*cos(yaw(1)) ;
Tbn(2,3)= sin(roll(1))*sin(yaw(1)) - cos(roll(1))*sin(pitch(1))*cos(yaw(1));

Tbn(3,1)=-sin(roll(1))*cos(pitch(1));
Tbn(3,2)=sin(pitch(1));
Tbn(3,3)=cos(roll(1))*cos(pitch(1)); 
Tbn=TT*Tbn;

wiep=[0;wie*cos(phi(1));wie*sin(phi(1))];%%%%%%
wepp=[-v2(1)/Re;v1(1)/Re;v1(1)*tan(phi(1))/Re];

Ve_true(1)=Ve0(1);
Vn_true(1)=Vn0(1);
lati_true(1)=lati0(1);
longi_true(1)=longi0(1);
pitch_true(1)=pitch0(1);
roll_true(1)=roll0(1);
yaw_true(1)=yaw0(1);

k=1;
T=0.7;
f=1/T;
q=[0,0,0,0]';
fn=zeros(3,N-1);
quat=zeros(4,N-1);
Tnb=Tbn';
for i=1:N-1
q(1)=cos(fai(k)/2)*cos(theta(k)/2)*cos(gama(k)/2) - sin(fai(k)/2)*sin(theta(k)/2)*sin(gama(k)/2);
q(2)=cos(fai(k)/2)*sin(theta(k)/2)*cos(gama(k)/2) - sin(fai(k)/2)*cos(theta(k)/2)*sin(gama(k)/2);
q(3)=cos(fai(k)/2)*cos(theta(k)/2)*sin(gama(k)/2) + sin(fai(k)/2)*sin(theta(k)/2)*cos(gama(k)/2);
q(4)=cos(fai(k)/2)*sin(theta(k)/2)*sin(gama(k)/2) + sin(fai(k)/2)*cos(theta(k)/2)*cos(gama(k)/2);   

quat(:,k)=q;
wib=[Gx(i),Gy(i),Gz(i)]';
wpbb=wib-inv(Tnb)*(wepp+wiep);

    
    omiga=[0        -wpbb(1) -wpbb(2) -wpbb(3);
        wpbb(1)  0         wpbb(3) -wpbb(2);
        wpbb(2) -wpbb(3)  0         wpbb(1);
        wpbb(3)  wpbb(2) -wpbb(1)  0       ]*T;
    

    
    derta = sqrt((omiga(1,2))^2+(omiga(1,3))^2+(omiga(1,4))^2);
    
    q = [eye(4)*(1-derta^2/8+derta^4/384)+(1/2-derta^2/48)*omiga]*q;%泰勒展开
    q=q/sqrt(q(1)^2+q(2)^2+q(3)^2+q(4)^2);
    
   Tnb = [q(1)^2+q(2)^2-q(3)^2-q(4)^2   2*(q(2)*q(3)+q(1)*q(4))       2*(q(2)*q(4)-q(1)*q(3));% > 用四元素表示得姿态矩阵
        2*(q(2)*q(3)-q(1)*q(4))        q(1)^2-q(2)^2+q(3)^2-q(4)^2   2*(q(3)*q(4)+q(1)*q(2));
        2*(q(2)*q(4)+q(1)*q(3))        2*(q(3)*q(4)-q(1)*q(2))      q(1)^2-q(2)^2-q(3)^2+q(4)^2];

    roll(i+1) = asin(Tnb(2,3)); %> 横滚角计算
    pitch(i+1) = atan(-Tnb(1,3)/Tnb(3,3)); %> 俯仰角计算
    if abs(Tnb(3,3))>eps
        pitch(i+1) = atan(-Tnb(1,3)/Tnb(3,3));
        if Tnb(3,3)>0
            pitch(i+1) = pitch(i+1);
        elseif -Tnb(1,3)> 0
            pitch(i+1) = pitch(i+1)+pi;
        else pitch(i+1) =pitch(i+1)-pi;
        end
    elseif -Tnb(1,3)> 0
        pitch(i+1) = pi/2;
    else pitch(i+1) = -pi/2;
    end
    yaw(i+1) = atan(-Tnb(2,1)/Tnb(2,2));% > 航向角计算
    if abs(Tnb(2,2))>eps
        yaw(i+1) = atan(Tnb(2,1)/Tnb(2,2));
        if Tnb(2,2)>0
            yaw(i+1) = yaw(i+1);
        elseif Tnb(2,1)> 0
            yaw(i+1) = yaw(i+1)+pi;
        else yaw(i+1) = yaw(i+1)-pi;
        end
    elseif Tnb(2,1)>0
        yaw(i+1) = pi/2;
    else yaw(i+1) = -pi/2;
    end
    
    k=k+1;  
    fai(k)=yaw(i+1);
    theta(k)=pitch(i+1);
    gama(k)=roll(i+1);
    
    Pitch(k)=pitch(i+1)/rad_deg;
    Roll(k)=roll(i+1)/rad_deg;
    Yaw(k)=yaw(i+1)/rad_deg;
    
    pitch_true(k)=pitch0(i+1);
    roll_true(k)=roll0(i+1);
    yaw_true(k)=yaw0(i+1);
    

    fb=[Ax(i),Ay(i),Az(i)]';
    fp=Tnb'*fb;
    fn(1:3,k-1)=fp;
    
    fv(1)=fp(1)+(2*wiep(3) + wepp(3))*v2(k-1);     %% velocity update
    fv(2)=fp(2)-(2*wiep(3) + wepp(3))*v1(k-1);
    Kv1=fv';
    
    tmp_v(1)=v1(k-1) + fv(1)*T;
    tmp_v(2)=v2(k-1) + fv(2)*T;
    
    fv(1)=fp(1)+(2*wiep(3) + wepp(3))*tmp_v(2);
    fv(2)=fp(2)-(2*wiep(3) + wepp(3))*tmp_v(1);
    Kv2=fv';

    Kv=(Kv1+Kv2)/2.0;
    v1(k)=v1(k-1) + Kv(1)*T;
    v2(k)=v2(k-1) + Kv(2)*T;

    Ve_true(k)=Ve0(i+1);   
    Vn_true(k)=Vn0(i+1);
    
    wepp=[-v2(k)/Re;
        v1(k)/Re;
        v1(k)*tan(phi(k-1))/Re];%%%%%%%%
    
    phi(k)=phi(k-1)+v2(k)*T/Re;%%%%%%%%
    lamda(k)=lamda(k-1)+v1(k)*T/(Re*cos(phi(k)));%%%%%%
    

    lati_true(k)=lati0(i+1);
    longi_true(k)=longi0(i+1);
    
    wiep=[0; wie*cos(phi(k)); wie*sin(phi(k))];%%%%%%%%
    g1=g0+0.051799*sin(phi(k))*sin(phi(k));%%%%%%%%%
end

figure(1);
subplot(2,1,1);
 plot((1:k)/3600/f,Ve_true(1:k),'-r');hold on;
plot((1:k)/3600/f,v1(1:k),'-g');hold on;grid on;
ylabel('Ve(m/s)');
legend('测量值','INS解算值');
subplot(2,1,2);
plot((1:k)/3600/f,Vn_true(1:k),'-r');hold on;
plot((1:k)/3600/f,v2(1:k),'-g');hold on;
grid on;xlabel('Time(h)');ylabel('Vn(m/s)');
legend('测量值','INS解算值');

figure(2);
subplot(2,1,1);
plot((1:k)/3600/f,lati_true(1:k),'-r');hold on;plot((1:k)/3600/f,phi(1:k),'-g');hold on;
grid on;ylabel('latitude(^o)');legend('GPS测量值','INS解算值');
subplot(2,1,2);plot((1:k)/3600/f,longi_true(1:k),'-r');hold on;
plot((1:k)/3600/f,lamda(1:k),'-g');hold on;
grid on;xlabel('Time(h)');ylabel('longitude(^o)');legend('GPS测量值','INS解算值');

figure(3);
subplot(3,1,1);plot((1:k)/3600/f,pitch_true(1:k),'-r');hold on;
plot((1:k)/3600/f,Pitch(1:k),'-g');hold on;
grid on;ylabel('Pitch(^o)');legend('测量值','INS解算值');
subplot(3,1,2);plot((1:k)/3600/f,roll_true(1:k),'-r');hold on;
plot((1:k)/3600/f,Roll(1:k),'-g');hold on;
grid on;ylabel('Roll(^o)');legend('测量值','INS解算值');
subplot(3,1,3);plot((1:k)/3600/f,yaw_true(1:k),'-r');hold on;plot((1:k)/3600/f,Yaw(1:k),'-g');hold on;
grid on;xlabel('Time(h)');ylabel('Yaw(^o)');legend('测量值','INS解算值');

figure(4);
plot((longi_true(1:k)*60*1852-longi_true(1:1)*60*1852),lati_true(1:k)*60*1852-lati_true(1:1)*60*1852,'-r');hold on;
plot((lamda(1:k)*60*1852-lamda(1:1)*60*1852),phi(1:k)*60*1852-phi(1:1)*60*1852,'-g');hold on;
set(gca,'xaxislocation','top');
set(gca,'yaxislocation','right');
grid on;ylabel(' E (m)');xlabel(' N (m)');
 legend('GPS输出值','INS解算值');


 


%惯导相关的噪声统计数据
Q_wg  = (70/(57*3600))^2;                       						%陀螺白噪声
Q_wr  = (200/(57*3600))^2;                      							%陀螺马氏过程
Q_wa  = (1e-5)^2;                                                           %加计马氏过程
Q 		= diag([Q_wg Q_wg Q_wg,  Q_wr Q_wr Q_wr,  Q_wa Q_wa Q_wa]);
Tr 		= 35*ones(3,1);
Ta 		= 100*ones(3,1);
R_p     =[(0.000001)^2,(0.000001)^2,(0.000001)^2];
R 		= diag(R_p);                        %计算测量噪声方差R					                        


v3=zeros(1,k);
Vd_true=zeros(1,k);
% v=[v1',v2',v3'];
% V_dvl=[Ve_true',Vn_true',Vd_true'];
% Dv=v-V_dvl;

dep=zeros(1,k);
dep_true=zeros(1,k);
Position_ins=[phi',lamda',dep'];
Position_gps=[lati_true',longi_true',dep_true'];
Dp=Position_ins-Position_gps;


%卡尔曼滤波参数初始化
PP(1:18,1:18) = diag([1/(36*57) 1/(36*57) 1/57, 0.0001 0.0001 0.0001, 0 0 1, 1/(57*3600) 1/(57*3600) 1/(57*3600), 0.04/(57*3600) 0.04/(57*3600) 0.04/(57*3600), 1e-4 1e-4 1e-4].^2);   %初始误差协方差阵
PP0 					= PP;
X 						= zeros(18,1);  %初始状态
E_attitude 		= zeros(1,3);
E_position 		= zeros(1,3);
E_velocity 		= zeros(1,3);
tao=T;


Ve=zeros(1,k);
Vn=zeros(1,k);
Vd=zeros(1,k);
h=zeros(1,k);
for i=1:k-1
    %参数赋值
    Ve 		= v1(i);
    Vn 		= v2(i);
    Vd 		= v3(i);
    L 		= phi(i);
    h 		= dep(i);
    Fe 		= fn(1,i);
    Fn 		= fn(2,i);
    Fd 		= fn(3,i);
    Rm 		= Re*(1-2*e+3*e*sin(L)^2);
    Rn 		= Re*(1-e*sin(L)^2);
    %由四元数计算姿态阵
    q 		= quat(:,i);
    Cnb 	= [1-2*(q(3)^2+q(4)^2),     2*(q(2)*q(3)-q(1)*q(4)), 2*(q(2)*q(4)+q(1)*q(3));
           	 2*(q(2)*q(3)+q(1)*q(4)), 1-2*(q(2)^2+q(4)^2),     2*(q(3)*q(4)-q(1)*q(2));
             2*(q(2)*q(4)-q(1)*q(3)), 2*(q(3)*q(4)+q(1)*q(2)), 1-2*(q(2)^2+q(3)^2)];

    %连续系统状态转换阵 F 的时间更新
    F            = zeros(18,18);
    F(1,2)       = wie*sin(L)+Ve*tan(L)/(Rn+h);
    F(1,3)       = -(wie*cos(L)+Ve/(Rn+h));
    F(1,5)       = -1/(Rm+h);
    F(2,1)       = -(wie*sin(L)+Ve*tan(L)/(Rn+h));
    F(2,3)       = -Vn/(Rm+h);
    F(2,4)       = 1/(Rn+h);
    F(2,7)       = -wie*sin(L);
    F(3,1)       = wie*cos(L)+Ve/(Rn+h);
    F(3,2)       = Vn/(Rm+h);
    F(3,4)       = tan(L)/(Rn+h);
    F(3,7)       = wie*cos(L)+Ve*(sec(L)^2)/(Rn+h);
    F(4,2)       = -Fd;
    F(4,3)       = Fn;
    F(4,4)       = Vn*tan(L)/(Rm+h)-Vd/(Rm+h);
    F(4,5)       = 2*wie*sin(L)+Ve*tan(L)/(Rn+h);
    F(4,6)       = -(2*wie*cos(L)+Ve/(Rn+h));
    F(4,7)       = 2*wie*cos(L)*Vn+Ve*Vn*sec(L)^2/(Rn+h)+2*wie*sin(L)*Vd;
    F(5,1)       = Fd;
    F(5,3)       = -Fe;
    F(5,4)       = -2*(wie*sin(L)+Ve*tan(L)/(Rn+h));
    F(5,5)       = -Vd/(Rm+h);
    F(5,6)       = -Vn/(Rm+h);
    F(5,7)       = -(2*wie*cos(L)+Ve*(sec(L)^2)/(Rn+h))*Ve;
    F(6,1)       = -Fn;
    F(6,2)       = Fe;
    F(6,4)       = 2*(wie*cos(L)+Ve/(Rn+h));
    F(6,5)       = 2*Vn/(Rm+h);
    F(6,7)       = -2*Ve*wie*sin(L);
    F(7,5)       = 1/(Rm+h);
    F(8,4)       = 1/((Rn+h)*cos(L));
    F(8,7)       = Ve*tan(L)/((Rn+h)*cos(L));
    F(9,6)       = 1;
    
    F(1:3,10:12) = Cnb;
    F(1:3,13:15) = Cnb;
    F(4:6,16:18) = Cnb;
    F(13,13)     = -1/Tr(1);
    F(14,14)     = -1/Tr(2);
    F(15,15)     = -1/Tr(3);
    F(16,16)     = -1/Ta(1);
    F(17,17)     = -1/Ta(2);
    F(18,18)     = -1/Ta(3);
    
    %连续系统输入矩阵更新
    G 				   = zeros(18,9);
    G(1:3,1:3)   = Cnb;
    G(13:15,4:6) = eye(3,3);
    G(16:18,7:9) = eye(3,3);
    %连续系统量测阵更新
    H 					= zeros(3,18);
    H(1,7) 			= 1;
    H(2,8) 			= 1;
    H(3,9) 			= 1;
    %连续系统离散化
    A 				  = eye(18,18)+F*tao;
    B 			    = (eye(18,18)+tao*F/2)*G*tao;
    
    %卡尔曼滤波
    P 					= A*(PP0)*A'+B*Q*B';
    K 					= P*H'*inv(H*P*H'+R);
    PP0 				= (eye(18,18)-K*H)*P;
    PP0 				= (PP0+PP0')/2;
    PP(i,:) 		= diag(PP0);
    
    z = Dp(i+1,:)';
    XX = A*X+K*(z-H*A*X);
    X = XX;

    E_attitude(i+1,:) = XX(1:3)';
    E_velocity(i+1,:) = XX(4:6)';
    E_position(i+1,:) = XX(7:9)';
end
E_ve=E_velocity(1:k,1)';
E_vn=E_velocity(1:k,2)';
E_phi=E_position(1:k,1)';
E_lamda=E_position(1:k,2)';
E_pitch=E_attitude(1:k,1)';
E_roll=E_attitude(1:k,2)';
E_yaw=E_attitude(1:k,3)';

figure(5);
subplot(2,1,1);plot((1:k)/3600/f,Ve_true(1:k),'-r');hold on;plot((1:k)/3600/f,(v1(1:k)-E_ve(1:k)),'-g');hold on;
grid on;ylabel('Ve(m/s)');legend('测量值','INS/GPS组合值');
subplot(2,1,2);plot((1:k)/3600/f,Vn_true(1:k),'-r');hold on;plot((1:k)/3600/f,v2(1:k)-E_vn(1:k),'-g');hold on;
grid on;xlabel('Time(h)');ylabel('Vn(m/s)');legend('测量值','INS/GPS组合值');

figure(6);
subplot(2,1,1);plot((1:k)/3600/f,lati_true(1:k),'-r');hold on;plot((1:k)/3600/f,(phi(1:k)-E_phi(1:k)),'-g');hold on;
grid on;ylabel('latitude(^o)');legend('GPS测量值','INS/GPS组合值');
subplot(2,1,2);plot((1:k)/3600/f,longi_true(1:k),'-r');hold on;plot((1:k)/3600/f,(lamda(1:k)-E_lamda(1:k)),'-g');hold on;
grid on;xlabel('Time(h)');ylabel('longitude(^o)');legend('GPS测量值','INS/GPS组合值');

figure(7);
plot((longi_true(1:k)*60*1852-longi_true(1:1)*60*1852),lati_true(1:k)*60*1852-lati_true(1:1)*60*1852,'-r');hold on;
plot(((lamda(1:k)-E_lamda(1:k))*60*1852-longi_true(1:1)*60*1852),(phi(1:k)-E_phi(1:k))*60*1852-lati_true(1:1)*60*1852,'-g');
set(gca,'xaxislocation','top');
set(gca,'yaxislocation','right');
grid on;ylabel(' E (m)');xlabel(' N (m)');legend('GPS输出真实值','INS/GPS组合值');

num_ins=0;
num_insgps=0;
 for i=1:k
   num_ins= num_ins+(phi(i)-lati_true(i))*(phi(i)-lati_true(i))+(lamda(i)-longi_true(i))*(lamda(i)-longi_true(i));  
   num_insgps= num_insgps+(phi(i)-E_phi(i)-lati_true(i))*(phi(i)-E_phi(i)-lati_true(i))+(lamda(i)-E_lamda(i)-longi_true(i))*(lamda(i)-E_lamda(i)-longi_true(i));
 end
 rmse_ins=sqrt(num_ins/k);
 rmse_insgps=sqrt(num_insgps/k);

disp('INS经纬度的RMSE(°）：');
disp(rmse_ins);
disp('INS/GPS经纬度的RMSE(°）：');
disp(rmse_insgps);

disp('INS距离的RMSE(m)：');
disp(rmse_ins*60*1852);
disp('INS/GPS距离的RMSE(m)：');
disp(rmse_insgps*60*1852);

