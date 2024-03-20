%  Mie theory for scattering, extinction, and absorption cross section
%  core-shell structure can be calculated
%  Wei E. I. Sha update Jan. 2019
%  email:weisha@zju.edu.cn

function Mie_LPR_Multilayer
clc;clear

% load metal data
Au_data=load('./data/Au_epr.txt');
Au_epr=Au_data(:,1)+j*Au_data(:,2);

% load active material
active_data=load('./data/P3HT_PCBM_epr.txt');
P3HT_PCBM_epr=active_data(:,1)+j*active_data(:,2);

%  parameter
epsc=1/(4*pi*9*10.^9);      %  permittivity in free space     
murc=4*pi*10.^(-7);         %  permeability in free space     
x=400:10:800;               %  wavelength in nm

%  relative permittivity
Air=ones(1,length(x))*1;                %  air
Si_epr=3.4^2*ones(1,length(x))*1;       %  silicon

%  wavelenth and frequency
wavelength=x*10^(-9);
omega=2*pi*3*10^8./wavelength;

MAX=4;               %  the number of Mie series (need increase for large sphere)
L=2;                 %  how many layer (including for the outmost layer)

for fre=1:length(wavelength)
    fre
    
    eps_array=[Au_epr(fre),Air(fre)]; %  relative permittivity (inner to outer)
    mu_array=[1,1];                   %  relative permeability
    r_array=[30];                     %  radius (inner to outer)
    
    k_array=sqrt(mu_array.*eps_array)*2*pi/(wavelength(fre));
    eps_array=eps_array*epsc;
    mu_array=mu_array*murc;
    r_array=r_array*10^(-9);
    
    %  inverse the order
    inv_k_arr=k_array(end:-1:1);
    inv_eps_arr=eps_array(end:-1:1);
    inv_mu_arr=mu_array(end:-1:1);
    inv_r_arr=r_array(end:-1:1);
    
    %  Outgoing wave reflection coef for TM wave (inner to outer)
    for m=1:L-1;
        for n=1:MAX
            R_TM_Out(m,n)=...
                (sqrt(eps_array(m+1)*mu_array(m))*ssb_der2(n,k_array(m),r_array(m))*ssb2(n,k_array(m+1),r_array(m))-...
                sqrt(eps_array(m)*mu_array(m+1))*ssb_der2(n,k_array(m+1),r_array(m))*ssb2(n,k_array(m),r_array(m)))/...
                (sqrt(eps_array(m)*mu_array(m+1))*ssb_der2(n,k_array(m+1),r_array(m))*ssb1(n,k_array(m),r_array(m))-...
                sqrt(eps_array(m+1)*mu_array(m))*ssb_der1(n,k_array(m),r_array(m))*ssb2(n,k_array(m+1),r_array(m)));
            T_TM_Out(m,n)=...
                -i*eps_array(m+1)*sqrt(mu_array(m+1)/eps_array(m))/...
                (sqrt(eps_array(m)*mu_array(m+1))*ssb_der2(n,k_array(m+1),r_array(m))*ssb1(n,k_array(m),r_array(m))-...
                sqrt(eps_array(m+1)*mu_array(m))*ssb_der1(n,k_array(m),r_array(m))*ssb2(n,k_array(m+1),r_array(m)));
        end
    end
    
    %  Standing wave reflection coef for TM wave (outer to inner)
    for m=1:L-1;
        for n=1:MAX
            R_TM_In(m,n)=...
                (sqrt(inv_eps_arr(m)*inv_mu_arr(m+1))*ssb_der1(n,inv_k_arr(m+1),inv_r_arr(m))*ssb1(n,inv_k_arr(m),inv_r_arr(m))-...
                sqrt(inv_eps_arr(m+1)*inv_mu_arr(m))*ssb_der1(n,inv_k_arr(m),inv_r_arr(m))*ssb1(n,inv_k_arr(m+1),inv_r_arr(m)))/...
                (sqrt(inv_eps_arr(m+1)*inv_mu_arr(m))*ssb_der2(n,inv_k_arr(m),inv_r_arr(m))*ssb1(n,inv_k_arr(m+1),inv_r_arr(m))-...
                sqrt(inv_eps_arr(m)*inv_mu_arr(m+1))*ssb_der1(n,inv_k_arr(m+1),inv_r_arr(m))*ssb2(n,inv_k_arr(m),inv_r_arr(m)));
            T_TM_In(m,n)=...
                -i*inv_eps_arr(m+1)*sqrt(inv_mu_arr(m+1)/inv_eps_arr(m))/...
                (sqrt(inv_eps_arr(m+1)*inv_mu_arr(m))*ssb_der2(n,inv_k_arr(m),inv_r_arr(m))*ssb1(n,inv_k_arr(m+1),inv_r_arr(m))-...
                sqrt(inv_eps_arr(m)*inv_mu_arr(m+1))*ssb_der1(n,inv_k_arr(m+1),inv_r_arr(m))*ssb2(n,inv_k_arr(m),inv_r_arr(m)));
        end
    end
    
    %  Standing wave reflection coef for TM wave£¨inner to outer£©
    R_TM_In=R_TM_In(end:-1:1,:);
    T_TM_In=T_TM_In(end:-1:1,:);
    
    %  Generalized reflection coef for TM wave
    TM_Ref(1,:)=R_TM_In(1,:);
    for m=2:L-1
        for n=1:MAX
            TM_Ref(m,n)=R_TM_In(m,n)+(T_TM_Out(m,n)*TM_Ref(m-1,n)*T_TM_In(m,n))/...
                (1-R_TM_Out(m,n)*TM_Ref(m-1,n));
        end
    end
    
    %  Outgoing wave reflection coef for TE wave (inner to outer)
    for m=1:L-1;
        for n=1:MAX
            R_TE_Out(m,n)=...
                (sqrt(eps_array(m)*mu_array(m+1))*ssb_der2(n,k_array(m),r_array(m))*ssb2(n,k_array(m+1),r_array(m))-...
                sqrt(eps_array(m+1)*mu_array(m))*ssb_der2(n,k_array(m+1),r_array(m))*ssb2(n,k_array(m),r_array(m)))/...
                (sqrt(eps_array(m+1)*mu_array(m))*ssb_der2(n,k_array(m+1),r_array(m))*ssb1(n,k_array(m),r_array(m))-...
                sqrt(eps_array(m)*mu_array(m+1))*ssb_der1(n,k_array(m),r_array(m))*ssb2(n,k_array(m+1),r_array(m)));
            T_TE_Out(m,n)=...
                -i*mu_array(m+1)*sqrt(eps_array(m+1)/mu_array(m))/...
                (sqrt(eps_array(m+1)*mu_array(m))*ssb_der2(n,k_array(m+1),r_array(m))*ssb1(n,k_array(m),r_array(m))-...
                sqrt(eps_array(m)*mu_array(m+1))*ssb_der1(n,k_array(m),r_array(m))*ssb2(n,k_array(m+1),r_array(m)));
        end
    end
    
    %  Standing wave reflection coef for TE wave (outer to inner)
    for m=1:L-1;
        for n=1:MAX
            R_TE_In(m,n)=...
                (sqrt(inv_eps_arr(m+1)*inv_mu_arr(m))*ssb_der1(n,inv_k_arr(m+1),inv_r_arr(m))*ssb1(n,inv_k_arr(m),inv_r_arr(m))-...
                sqrt(inv_eps_arr(m)*inv_mu_arr(m+1))*ssb_der1(n,inv_k_arr(m),inv_r_arr(m))*ssb1(n,inv_k_arr(m+1),inv_r_arr(m)))/...
                (sqrt(inv_eps_arr(m)*inv_mu_arr(m+1))*ssb_der2(n,inv_k_arr(m),inv_r_arr(m))*ssb1(n,inv_k_arr(m+1),inv_r_arr(m))-...
                sqrt(inv_eps_arr(m+1)*inv_mu_arr(m))*ssb_der1(n,inv_k_arr(m+1),inv_r_arr(m))*ssb2(n,inv_k_arr(m),inv_r_arr(m)));
            T_TE_In(m,n)=...
                -i*inv_mu_arr(m+1)*sqrt(inv_eps_arr(m+1)/inv_mu_arr(m))/...
                (sqrt(inv_eps_arr(m)*inv_mu_arr(m+1))*ssb_der2(n,inv_k_arr(m),inv_r_arr(m))*ssb1(n,inv_k_arr(m+1),inv_r_arr(m))-...
                sqrt(inv_eps_arr(m+1)*inv_mu_arr(m))*ssb_der1(n,inv_k_arr(m+1),inv_r_arr(m))*ssb2(n,inv_k_arr(m),inv_r_arr(m)));
        end
    end
    
    %  Standing wave reflection coef for TE wave (inner to outer)
    R_TE_In=R_TE_In(end:-1:1,:);
    T_TE_In=T_TE_In(end:-1:1,:);
    
    %  Generalized reflection coef for TE wave
    TE_Ref(1,:)=R_TE_In(1,:);
    for m=2:L-1
        for n=1:MAX
            TE_Ref(m,n)=R_TE_In(m,n)+(T_TE_Out(m,n)*TE_Ref(m-1,n)*T_TE_In(m,n))/...
                (1-R_TE_Out(m,n)*TE_Ref(m-1,n));
        end
    end
    
    %  total contribution
    sum1=0;
    sum2=0;
    for n=1:MAX
        sum1=sum1+(2*n+1)*(abs(TM_Ref(end,n))^2+abs(TE_Ref(end,n))^2);
        sum2=sum2-(2*n+1)*real(TM_Ref(end,n)+TE_Ref(end,n));
    end
    
    %  tm contribution
    sum1_tm=0;
    sum2_tm=0;
    for n=1:MAX
        sum1_tm=sum1_tm+(2*n+1)*(abs(TM_Ref(end,n))^2);
        sum2_tm=sum2_tm-(2*n+1)*real(TM_Ref(end,n));
    end
    
    %  te contribution
    sum1_te=0;
    sum2_te=0;
    for n=1:MAX
        sum1_te=sum1_te+(2*n+1)*(abs(TE_Ref(end,n))^2);
        sum2_te=sum2_te-(2*n+1)*real(TE_Ref(end,n));
    end
    
    % scattering, extinction, and absorption cross section
    sca=sum1*2*pi/abs(k_array(end))^2;
    ext=sum2*2*pi/abs(k_array(end))^2;
    absc=ext-sca;
        
    sca_tm=sum1_tm*2*pi/abs(k_array(end))^2;
    ext_tm=sum2_tm*2*pi/abs(k_array(end))^2;
    absc_tm=ext_tm-sca_tm;
    
    sca_te=sum1_te*2*pi/abs(k_array(end))^2;
    ext_te=sum2_te*2*pi/abs(k_array(end))^2;
    absc_te=ext_te-sca_te;    
    
    % normalized to geometric size
    Absc(fre)=absc/pi/r_array(end)^2;
    Sca(fre)=sca/pi/r_array(end)^2;
    Ext(fre)=ext/pi/r_array(end)^2;
       
    Absc_tm(fre)=absc_tm/pi/r_array(end)^2;
    Sca_tm(fre)=sca_tm/pi/r_array(end)^2;
    Ext_tm(fre)=ext_tm/pi/r_array(end)^2;
    
    Absc_te(fre)=absc_te/pi/r_array(end)^2;
    Sca_te(fre)=sca_te/pi/r_array(end)^2;
    Ext_te(fre)=ext_te/pi/r_array(end)^2;
    
end

figure(1)

hold on;
plot(x,(Ext),'bo')
plot(x,(Absc),'kd')
%plot(x,Ext_tm,'r.-')
%plot(x,Ext_te,'g--')

xlabel('Wavelength (nm)')
ylabel('ECS and ACS')

% Spherical bessel 1
function output=sb1(degree,k,r)
output=sqrt(pi/(2*k*r))*besselj(degree+0.5,k*r);

% Modified spherical bessel 1
function output=ssb1(degree,k,r)
output=k*r*sb1(degree,k,r);

% Derivative of modified spherical bessel 1
function output=ssb_der1(degree,k,r)
output=sb1(degree,k,r)+k*r*(degree*sb1(degree-1,k,r)-(degree+1)*sb1(degree+1,k,r))/(2*degree+1);

% Spherical bessel 2
function output=sb2(degree,k,r)
output=sqrt(pi/(2*k*r))*besselh(degree+0.5,2,k*r);

% Modified spherical bessel 2
function output=ssb2(degree,k,r)
output=k*r*sb2(degree,k,r);

% Derivative of modified spherical bessel 2
function output=ssb_der2(degree,k,r)
output=sb2(degree,k,r)+k*r*(degree*sb2(degree-1,k,r)-(degree+1)*sb2(degree+1,k,r))/(2*degree+1);

% Associated legendre
function output=legen(degree,order,input)
var=legendre(degree,input);
output=var(order+1);

% Derivative of associated legendre
function output=legen_der(degree,order,input)
output=((degree+1)*input*legen(degree,order,input)-(degree-order+1)*legen(degree+1,order,input))/(1-input*input);
