clear all
close all
clc
global line voltage
line(:,10:13)=xlsread('lines-901207-line48outage.xls');
voltage(:,9:11)=xlsread('voltages-901207-line48outage.xls');
for i=1:length(line)
    line(i,1)=line(i,10);
    line(i,2)=line(i,11);
    line(i,3)=(line(i,12)+j*line(i,13))/182.25;
end
for i=1:length (voltage)
    voltage(i,1)=voltage(i,10)*cos(voltage(i,11)*pi/180)+j*voltage(i,10)*sin(voltage(i,11)*pi/180);
end