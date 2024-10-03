clear all;close all;clc
cd 'C:\MATLAB\data_fra_ecmwf\data'

info = ncinfo('cldcov_2017_12_anal.nc');
lon = ncread('cldcov_2017_12_anal.nc', 'longitude'); %længdegrad
lat = ncread('cldcov_2017_12_anal.nc', 'latitude'); %breddegrad
tcc = ncread('cldcov_2017_12_anal.nc', 'tcc'); %total cloud cover
time = ncread('cldcov_2017_12_anal.nc', 'time'); %time

info2 = ncinfo('cldcov_2017_11_anal.nc');
lon2 = ncread('cldcov_2017_11_anal.nc', 'longitude'); %længdegrad
lat2 = ncread('cldcov_2017_11_anal.nc', 'latitude'); %breddegrad
tcc2 = ncread('cldcov_2017_11_anal.nc', 'tcc'); %total cloud cover
time2 = ncread('cldcov_2017_11_anal.nc', 'time'); %time

%tccm=zeros(120,61,12)
%for i=1:9
%    tcc = ncread(['cldcov_2017_0i_anal.nc', 'tcc']);

%tccm=zeros(120,61,12)

%tccm(:,:,i)= mean(tcc,3)
%end

