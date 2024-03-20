% clc;clear
% combine data
% by Wei E.I. Sha (Email: weisha@zju.edu.cn)
% Update at Jan. 2019

number=41;         % total frequency points
wavelength=[];     % wavelength (nm)
active_abs=[];     % metal absorption

% all frequency
for m=1:number;
    a=sprintf('./result/Ene_Active%d.txt',m-1);
    data=load(a);
    wavelength=[wavelength,data(1)];
    active_abs=[active_abs,data(2)];
end

figure(1)
plot(wavelength,active_abs,'r.-')
xlabel('wavelength (nm)')
ylabel('Absorption of active layer per unit volume (a.u.)')


    