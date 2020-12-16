% Auto_phase_correction.m
% Maria Yanez Lopez, KCL 2019.
% 
% USE:
% [out,outw,degree, N_flat]=Final_auto_phase_correction(in,inw);
% 
% DESCRIPTION:
% Perform automatic phase correction
% 
% INPUTS:
% in     = water suppressed input data in matlab structure format.
% inw    = water unsuppressed input data in matlab structure format.
%
% OUTPUTS:
% out    = Water suppressed output following phase correction  
% outw   = Water unsuppressed output following phase correction
% degree   = Phase correction in degrees following zero order phase correction

function [out,outw,degree, N_flat]=Auto_phase_correction(in,inw)

%First step: making baseline flat
tmp=in;
shifting_points=50;%default is 50
Residual = zeros(length(in.fids),shifting_points);%preallocate

for n=1:shifting_points  
    tmp.fids = circshift(in.fids,-n);
    tmp.specs = fftshift(ifft(tmp.fids));
    x_to_fit = (tmp.ppm)';
    y_to_fit = real(tmp.specs);
    %Exclude datapoints not from baseline (metabolite + lipids + water peaks)
    a=0; %default 0ppm
    b=6; %default 6ppm
    [ ~, ix_a ] = min( abs( x_to_fit-a ) );
    [ ~, ix_b ] = min( abs( x_to_fit-b ) );
    x_to_fit(ix_b:ix_a,1) = NaN; 
    y_to_fit(ix_b:ix_a,1) = NaN;
    idx = isnan(x_to_fit);
    idy = isnan(y_to_fit);
    %Fit remaining points to flat line (polynomial zero order)
    P = polyfit(x_to_fit(~idx),y_to_fit(~idy),0);
    yfit = P(1);
    %Calculate residuals
    Residual(:,n)= (y_to_fit - yfit').^2;
end
%Find n that makes baseline flattest
columns= nansum(Residual, 1);
N_flat = find((columns)==min((columns))); 

%Shifting first N points in the FID and recalculate spectrum
in.fids = circshift(in.fids,-N_flat);
in.specs=fftshift(ifft(in.fids));
inw.fids = circshift(inw.fids,-N_flat);
inw.specs=fftshift(ifft(inw.fids));

%Second step: calculating zero order phase automatically
phis = 0:0.01:2*pi;%range in radians
degrees = phis*180/pi; %range in degrees
%Automatic phasing of water suppressed spectrum
out_ph=in;
out_ph1=in;
Residuals2=zeros(length(in.fids),1,size(phis,2));

for phindex=1:size(phis,2)
    phi=phis(phindex);
    out_ph.specs=((real(in.specs)).*cos(phi)-(imag(in.specs)).*sin(phi))+1i*((real(in.specs)).*sin(phi)+(imag(in.specs)).*cos(phi));
    Residuals2(:,1,phindex) = (real(out_ph.specs) - abs(out_ph.specs)).^2;
end
Residuals2=squeeze(Residuals2);
%Take into account only the metabolic peaks
c=1.7; %default 1.7ppm
d=4; %default 4ppm
[ ~, ix_d ] = min( abs( in.ppm-d ) );
[ ~, ix_c ] = min( abs( in.ppm-c ) );
Residuals2(1:ix_d,:)=NaN; 
Residuals2(ix_c:end,:)=NaN;
columns= nansum(Residuals2, 1);
N = find((columns)==min((columns)));
degree=degrees(1,N); %Zero order phase in degrees
phi=phis(1,N); %Zero order phase in radians

%Phasing spectrum and recalculating FID
out_ph1.specs=((real(in.specs)).*cos(phi)-(imag(in.specs)).*sin(phi))+1i*((real(in.specs)).*sin(phi)+(imag(in.specs)).*cos(phi));
out_ph1.fids=fftshift(ifft(out_ph1.specs));
out=out_ph1;

%Automatic phasing of water spectrum
outw_ph=inw;
Residuals3=zeros(length(inw.fids),1,size(phis,2));
for phindex=1:size(phis,2)
    phi=phis(phindex);
    outw_ph.specs=((real(inw.specs)).*cos(phi)-(imag(inw.specs)).*sin(phi))+1i*((real(inw.specs)).*sin(phi)+(imag(inw.specs)).*cos(phi));
    Residuals3(:,1,phindex) = (real(outw_ph.specs) - abs(outw_ph.specs)).^2;
end
Residuals3=squeeze(Residuals3);
columns= nansum(Residuals3, 1);
N = find((columns)==min((columns)));
phi=phis(1,N); %Zero order phase in radians

%Phasing water spectrum and recalculating FID
outw_ph.specs=((real(inw.specs)).*cos(phi)-(imag(inw.specs)).*sin(phi))+1i*((real(inw.specs)).*sin(phi)+(imag(inw.specs)).*cos(phi));
outw_ph.fids=fftshift(ifft(outw_ph.specs));
outw=outw_ph;


