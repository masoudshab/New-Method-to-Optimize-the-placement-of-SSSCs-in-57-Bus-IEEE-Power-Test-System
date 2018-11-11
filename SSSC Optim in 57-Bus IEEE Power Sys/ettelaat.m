clear all
close all
clc
global line voltage
line(:,10:13)=xlsread('lines-9011-01.xls');
voltage(:,10:11)=xlsread('voltage9011-01.xls');
for i=1:length(line)
    line(i,1)=line(i,10);
    line(i,2)=line(i,11);
    line(i,3)=line(i,12)+j*line(i,13);
end
for i=1:length (voltage)
voltage(i,1)=voltage(i,10)*cos(voltage(i,11))+j*voltage(i,10)*sin(voltage(i,11));
end