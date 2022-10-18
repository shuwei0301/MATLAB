clear;
Cs=1500;               %%假设声纳在水中的速度为1500m/s
V=24*0.5144444;        %%假设船速为24节，又每节为0.514444m/s
rad_deg=0.01745329;    %%度转换成弧度的比例因子pi/180
q1=-180:1:180;         %%目标舷角，单位为角度
H=20;                  %%假设目标深度为20m
D_jing=[100,200,300,400,500];%%目标距离
q=q1*rad_deg;          %%目标舷角，单位为弧度
for i=1:5
deita(i)=asin(H/D_jing(i));%%目标高低角，与目标距离和深度有关
    for j=1:361        %%从-180到180共有361个数据
if (-90<=q1(j))&&(q1(j)<=90)%%当目标舷角度处于-90°到90°之间时
%%D1(i)=((Cs^2-V^2)*T(j))/2/(Cs+V*cos(q(j))*cos(deita(i)));
T(j)=2*D_jing(i)*(Cs+V*cos(q(j))*cos(deita(i)))/(Cs^2-V^2);%%T为声波发出到接收的总时间，受目标舷角与目标距离的影响
t_jian(j)=Cs*T(j)/2/(Cs+V*cos(q(j))*cos(deita(i)));
D_jian(j)=Cs^2*T(j)/2/(Cs+V*cos(q(j)));
else
%%D1(i)=((Cs^2-V^2)*T(j))/2/(Cs-V*cos(q(j))*cos(deita(i)));
T(j)=2*D_jing(i)*(Cs-V*cos(q(j))*cos(deita(i)))/(Cs^2-V^2);
t_jian(j)=Cs*T(j)/2/(Cs-V*cos(q(j))*cos(deita(i)));
D_jian(j)=Cs^2*T(j)/2/(Cs-V*cos(q(j)))
end
t_jing(i)=D_jing(i)/Cs;        %%精确模型的反射段时间
t_chuan(j)=T(j)/2;          %%传统模型的反射段时间
t_error1(j,i)=t_chuan(j)-t_jing(i);%%传统模型相比于精确模型的反射段时间误差
t_error2(j,i)=t_jian(j)-t_jing(i);%%简化工程模型相对与精确模型的反射段时间误差
D_chuan(j)=Cs*T(j)/2;       %%传统模型的目标距离
D_error1(j,i)=D_chuan(j)-D_jing(i);%%传统模型相比于精确模型的目标距离误差
D_error2(j,i)=D_jian(j)-D_jing(i);%%简化模型相对与精确模型的目标距离误差
    end
end
figure;plot(q1,D_error1);legend('100m','200m','300m','400m','500m');ylabel('目标距离计算误差（m）');xlabel('目标舷角（°）');title('传统模型相比精确模型的距离计算误差')
figure;plot(q1,t_error1);legend('100m','200m','300m','400m','500m');ylabel('反射段时间计算误差（s）');xlabel('目标舷角（°）');title('传统模型相比精确模型的反射时间计算误差')
figure;plot(q1,D_error2);legend('100m','200m','300m','400m','500m');ylabel('目标距离计算误差（m）');xlabel('目标舷角（°）');title('简化模型相比精确模型的距离计算误差')
figure;plot(q1,t_error2);legend('100m','200m','300m','400m','500m');ylabel('反射段时间计算误差（s）');xlabel('目标舷角（°）');title('简化模型相比精确模型的反射时间计算误差')