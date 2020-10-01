fgs = readtable('FGS_001.txt', 'HeaderLines', 1);

x1=fgs{:,1}; y1=fgs{:,2}; 
x2=weno{:,1}; y2=weno{:,2}; 
x3=ocrwenold{:,1}; y3=ocrwenold{:,2}; 

plot(x1,y1,'black','LineWidth',2)
hold on
plot(x2,y2,'-o','MarkerIndices',1:4:length(y2),'LineWidth',1.5)
plot(x3,y3,'-.*','MarkerIndices',1:3:length(y3),'LineWidth',1.5);

xlabel('X','FontSize',18)
ylabel('\rho','FontSize',20)
legend('Fine Grid Solution','WENO5, N=400','OCRWENO-LD, N=400','FontSize',13)

