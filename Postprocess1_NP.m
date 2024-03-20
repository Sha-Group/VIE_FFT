% clc;clear
% combine data
% by Wei E.I. Sha (Email: weisha@zju.edu.cn)
% Update at Jan. 2019

number=41;         % total frequency points
wavelength=[];     % wavelength (nm)
metal_ext=[];      % metal extinction
metal_abs=[];      % metal absorption

% all frequency
for m=1:number;
    a=sprintf('./result/Ene_Metal%d.txt',m-1);
    data=load(a);
    wavelength=[wavelength,data(1)];
    metal_ext=[metal_ext,data(2)];
    metal_abs=[metal_abs,data(3)];
end

figure(1)
plot(wavelength,metal_ext,'r.-')
hold on;
plot(wavelength,metal_abs,'k--')
xlabel('wavelength (nm)')
ylabel('ECS and ACS (dB)')
legend('ECS','ACS')

    