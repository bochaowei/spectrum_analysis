
%%
open('throughputLog_2022.01.20.fig');
ax = gca; 
h = findobj(gca,'Type','line');
x2 = h.XData; 
y2 = h.YData;
semilogy(x,y,x2,y2,'.')
%%
open('throughputLog_2022.01.24c.fig');
ax = gca; 
h = findobj(gca,'Type','line');
x = h.XData; 
y = h.YData;
%y(1:2500)=y2(1:2500);
%semilogy(x,y,'.','color','blue')
semilogy(x(32:4096)-x(32),y(32:4096),'.','color','blue')
ax = gca;
ax.FontSize = 16;
ax.FontWeight='bold';
ax.LineWidth = 1;
box on
ylabel('Throughput (atoms/s)')
xlabel('Time (hour)')
ylim([2*10^10,1*10^14])
xlim([0,66.6])
xticks([0,10,20,30,40,50,60])



%%calculate the amount of Rb
sum=0;
x(isnan(y))=[];
y(isnan(y))=[];
x(isinf(y))=[];
y(isinf(y))=[];
for i=1:length(x)-1
    sum=sum+(x(i+1)-x(i))*(y(i+1)+y(i))/2;
end
sum/10^12;
sum*3600/6.02/10^23*85.4678*1000