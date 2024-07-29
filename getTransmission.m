function out = getTransmission(Filename)
data = xlsread(Filename);
data_wavelength = data(:,1);
data_absorption = data(:,2); %%传输1m吸收
data_scattering = data(:,3); %%传输1m散射


% interp_wavelength = min(data_wavelength) : 0.5 : max(data_wavelength);
% interp_transmission = interp1(data_wavelength,data_transmission,interp_wavelength,'spline');
% out = [interp_wavelength', interp_transmission'];


interp_wavelength = min(data_wavelength) : 0.5 : max(data_wavelength);
interp_absorption = interp1(data_wavelength,data_absorption,interp_wavelength,'spline');
interp_scattering = interp1(data_wavelength,data_scattering,interp_wavelength,'spline');

interp_loss = interp_absorption + interp_scattering;
out = [interp_wavelength', interp_loss'];





