clc; clear all; close all;
cputime=0;
tic;
data_transmission = getTransmission('attenuation.xlsx');
beam_size = 1e-3;
beam_divergence = 1e-3;
energy = 1e-6;
pulse_width = 0.2e-12;
Po=energy / pulse_width / pi *(beam_size^2 /4); %input pwr in watts;
lambda_array = 300e-9 : 0.5e-9 : 800e-9;
rho_medium = 1000; %1000 kg/m^3
T_medium = 293.15;
% transmission = 0.95; %%透过率
wavelength = 355;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%wavelength
[n_lambda, b2, b3] = getNGVDTODfromLorentzLorenz(wavelength * 1e-9, rho_medium, T_medium);
index = find(data_transmission(:,1) == wavelength);
loss = data_transmission(index, 2);
ln=1;
i=sqrt(-1);
% loss = 0.05
% Po=1; %input pwr in watts
% transmission = 0.05; %%透过率
alpha=0.00; % Fiber loss value in dB/km
alph=alpha/(4.343); %Ref page#55 eqn 2.5.3 Fiber optic Comm by GP Agrawal
gamma=0.00; %fiber non linearity in /W/m
to=pulse_width/(2 * sqrt(2 * log(2))); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%initial pulse width in second
C=0; %Input chirp parameter for first calculation
% b2=51.671685e-27; %2nd order disp. (s2/m)
% Ld=(to^2)/(abs(b2)); %dispersion length in meter
Ld = 1/n_lambda;
pi=3.1415926535;
Ao=sqrt(Po); %Amplitude

%----------------------------------------------------------
dt = 0.001 * to;
tau = -50*to : dt : 50 * to;
rel_error=1e-5;
h=0.01;% step size
num_h = 1 / h + 1;
h = 1 / num_h;
% iii_step = 0.1;
% iii = 0.001:iii_step:1.001;

index_op_pulse = 0;
% for ii_index=1 : length(iii) %the various fiber lengths can be varied and this vector can be changed
% ii = iii(ii_index);
z=1.0*Ld;
h = h * z;
u=Ao*exp(-((1+i*(-C))/2)*(tau/to).^2);%page#47 G.P.AGrawal
u_initial = u;
%         figure(1)
%         plot(tau, abs(u),'r');
%         title('Input Pulse'); xlabel('Time'); ylabel('Amplitude');
%         grid on;
%         hold on;

l=max(size(u));
%%%%%%%%%%%%%%%%%%%%%%%
fwhm1=find(abs(u)>abs(max(u)/2));
fwhm1=length(fwhm1);
dw=1/l/dt*2*pi;
w=(-1*l/2:1:l/2-1)*dw;
u=fftshift(u);
w=fftshift(w);
spectrum=fft(fftshift(u)); %spectrum
spectrum_initial = spectrum;
phase_initial = angle(spectrum_initial);
array_phase = [];
xx = h : h : z;
xx = xx - h;
for jj=h:h:z
        spectrum=spectrum.*exp(-alph/2*(h/2)+i*b2/2*w.^2*(h/2) + i/6 * b3 * w.^3 * (h/2)) ;
    f=ifft(spectrum);
%     f = f .* exp(-alph/2*(h/2)+i*b2/2*w.^2*(h/2) + i/6 * b3 * w.^3 * (h/2)) ;
    f=f.*exp(i*gamma*((abs(f)).^2)*(h));
    f = f *exp(-loss*h);
    f = f * beam_size^2/(beam_size + h * beam_divergence)^2;
    spectrum=fft(f);
    spectrum=spectrum.*exp(-alph*(h/2)+i*b2/2*w.^2*(h/2) + i/6 * b3 * w.^3 * (h/2)) ;
    
    index_op_pulse = index_op_pulse + 1;

    temp_pulse_profile=ifft(spectrum);
    op_pulse(index_op_pulse,:)=abs(temp_pulse_profile);%saving output pulse at all intervals
    % op_pulse(ln,:)=abs(f);%saving output pulse at all intervals%%把透过率算上

    fwhm=find(abs(f)>abs(max(f)/2));
    fwhm=length(fwhm);
    %% PBR ratio
    ratio=fwhm/fwhm1; %PBR at every value
    pbratio(index_op_pulse)=ratio;%saving PBR at every step size
    %%phase
    array_phase(index_op_pulse,:) = angle(spectrum);
    dd=atand((abs(imag(f)))/(abs(real(f))));
    phadisp(index_op_pulse)=dd;%saving pulse phas


end
%     f=ifft(spectrum);
%     op_pulse(index_op_pulse,:)=abs(f);%saving output pulse at all intervals
% op_pulse(ln,:)=abs(f);%saving output pulse at all intervals%%把透过率算上

%     fwhm=find(abs(f)>abs(max(f)/2));
%     fwhm=length(fwhm);
%     ratio=fwhm/fwhm1; %PBR at every value
%     pbratio(ii_index)=ratio;%saving PBR at every step size
%     dd=atand((abs(imag(f)))/(abs(real(f))));
%     phadisp(ii_index)=dd;%saving pulse phase
%     ln=ln+1;
% end
toc;
cputime=toc;
figure
x = tau*1e12;
y = h:h:z;
y = y - h;
[X, Y] = meshgrid(x, y);
mesh(X, Y, op_pulse(:,:));



% title('Pulse Evolution');
set(gca,'FontName','Times New Roman','FontSize',28);

xlabel('Time (ps)','FontName','Times New Roman','FontSize',28);
ylabel('Distance (m)','FontName','Times New Roman','FontSize',28);
zlabel('amplitude','FontName','Times New Roman','FontSize',28);
colormap jet
% colormap(slanCM('thermal'))
% figure
% [X, Y] = meshgrid(x, y);
% mesh(X, Y, array_phase(:,:));
% % title('Pulse Evolution');
% set(gca,'FontName','Times New Roman','FontSize',30);
%
% xlabel('Time (ps)','FontName','Times New Roman','FontSize',30);
% ylabel('Distance (m)','FontName','Times New Roman','FontSize',30);
% zlabel('amplitude','FontName','Times New Roman','FontSize',30);
% colormap jet













view([0 90])
xlim([-2 2])
ylim([0 0.7411])
yticks([0 0.2 0.4 0.6])
yticklabels({'0','0.2','0.4','0.6'})
% colorbar
% figure
% waterfall(X, Y, op_pulse(:,:));



figure
plot(xx, pbratio(1:index_op_pulse),'k');
xlabel('Number of steps');
ylabel('Pulse broadening ratio');
grid on;
hold on;
figure
plot(xx ,phadisp(1:index_op_pulse),'k');
xlabel('distance travelled');
ylabel('phase change');
grid on;
hold on;
disp('CPU time:'), disp(cputime);
xx_pbratio = h:h:z;
xx_pbratio = xx_pbratio - h;
% xx_pbratio = linspace(0, Ld, length(iii))';
pbratio = pbratio';
xy_pbratio =[xx_pbratio, pbratio'];


figure
plot(tau, abs(u_initial),'b');
hold on
plot(tau, abs(op_pulse(index_op_pulse, :)),'r')
title('Input Pulse Temporal Width'); xlabel('Time'); ylabel('Amplitude');
grid on;
hold on;



% 时空光场
beam_size_initial = beam_size;
beam_size_end = beam_size + Ld * beam_divergence;
sigma_initial = beam_size / (2 * sqrt(2 * log(2)));
sigma_end = beam_size_end / (2 * sqrt(2 * log(2)));

x_domain = linspace(-beam_size_end * 3, beam_size_end * 3, 1001);
pulse_initial = getGaussian(x_domain, 0, sigma_initial);
pulse_end = getGaussian(x_domain, 0, sigma_end);
figure
plot(x_domain, pulse_initial);
title('initial pulse beam size');
hold on
plot(x_domain, pulse_end);
colormap jet


profile_spatial_temporal_initial =  pulse_initial' * op_pulse(1,:);  %%看最下面解释
profile_spatial_temporal_end =  pulse_end' * op_pulse(end,:);  %%看最下面解释

figure
imagesc(tau*1e12,x_domain * 1e3,profile_spatial_temporal_initial);
xlim([-4 4])
xlabel('Time / ps','FontName','Times New Roman','FontSize',28);
ylabel('y / mm','FontName','Times New Roman','FontSize',28);
title('x=0时空光场','FontName','Times New Roman','FontSize',28);
colormap jet
set(gca,'FontName','Times New Roman','FontSize',28);

figure
imagesc(tau*1e12,x_domain * 1e3,profile_spatial_temporal_end);
xlim([-4 4])

xlabel('Time / ps','FontName','Times New Roman','FontSize',28);
ylabel('y / mm','FontName','Times New Roman','FontSize',28);
title('x=Ld时空光场','FontName','Times New Roman','FontSize',28);
colormap jet
set(gca,'FontName','Times New Roman','FontSize',28);


x = tau*1e12;
y = xx_pbratio;
z = x_domain * 1e3;
% 
% xxx = x(1 : 1000 : end);   %%%脉冲宽度
% yyy = y(1 : 1 : end);   %%%传输距离
% zzz = z(1 : 10 : end);  %%%光斑尺寸
% length_xxx = length(xxx);
% length_yyy = length(yyy);
% length_zzz = length(zzz);
% 
% [X,Y,Z] = meshgrid(xxx,yyy,zzz);
% 
% V = X * 0 + Y * 0 + Z * 0;  % 101  1001  101
% temp_initial = profile_spatial_temporal_initial(1:10:end, 1:1000:end);
% temp_end = profile_spatial_temporal_end(1:10:end, 1:1000:end);
% % array_3d_profile_spatial_temporal_initial = zeros(length_xxx, 1, length_zzz);
% % array_3d_profile_spatial_temporal_end = zeros(length_xxx, 1, length_zzz);
% for i = 1 : size(temp_initial, 1)
%     for j = 1 : size(temp_initial, 2)
%             V( 1, i, j)  = temp_initial(j,i);
%             V( end, i, j) = temp_end(j,i);
%     end
% end
% % V(:, 1, :) = temp_initial;
% % V(:, end, :) = temp_end;
% 
% % V(:, 1, :) = profile_spatial_temporal_initial(1:10:end, 1:100:end);
% % V(:, end, :) = profile_spatial_temporal_end(1:10:end, 1:100:end);
% 
% pulse_distance = op_pulse;
% temp_pulse_distance = pulse_distance(:, 1 : 1000 : end);
% pulse_distance = zeros(1, length_yyy, length_zzz);
% for i = size(temp_pulse_distance, 1)
%     for j = size(temp_pulse_distance, 2)
%         V( i, j, 1 + floor(length_xxx/2)) = temp_pulse_distance(i, j);
%     end
% end
% 
% % V(1 + floor(length(zzz)/2), :, :) = pulse_distance(1 : 100 : end, :);
% % V(1 + floor(length(zzz)/2), :, :) = pulse_distance;
% xslice = [];
% yslice = [yyy(1), yyy(end)];
% 
% zslice = 0;
% figure
% slice(X,Y,Z,V,xslice,yslice,zslice)
% shading interp