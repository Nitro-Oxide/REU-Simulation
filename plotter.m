figure (1)
rawTable = readtable('Data_REU.xlsx','Sheet','Sheet1');
plot(rawTable,1,2)
hold on
plot(rawTable,1,3)

legend('Service Rate','Arrival Rate')

figure(2)
rawTable2 = readtable('Data_REU.xlsx','Sheet','Sheet2');
plot(rawTable2,1,2)
hold on
plot(rawTable2,1,3)
legend('Service Rate','Arrival Rate')
