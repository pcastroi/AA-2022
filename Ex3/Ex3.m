%% Ar.A - Ex 3
clear variables
clear mex
close all
clc

f_bands = [125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000];

% r1noabs, r2noabs, r3noabs, r1abs, r2abs, r3abs
RT30_empty = [7.01 7.51 7.2; 8.19 7.81 8.2; 6.54 6.88 6.64; 6.69 7 6.78; 7.05 7.11 7.26; 6.83 6.78 6.86; 6.19 6.25 6.2; 6.22 6.28 6.24; 5.79 5.76 5.69; 5.12 5.24 5.12; 4.75 4.84 4.79; 4.28 4.27 4.34; 3.77 3.77 3.73; 3.4 3.16 3.33; 2.56 2.62 2.62; 1.95 1.95 1.9];
RT30_abs = [3.86 2.99 4.05; 3.84 3.34 3.24; 3.01 3.32 3.53; 2.68 2.89 2.79; 2.3 2.4 2.49; 2.3 2.4 2.16; 2.12 2.06 2.11; 2.32 2.18 2.06; 2 2.05 2.14; 1.86 1.94 1.95; 1.84 1.91 1.86; 1.77 1.91 1.83; 1.57 1.66 1.63; 1.52 1.61 1.6; 1.3 1.36 1.28; 1.1 1.13 1.1];
for i=1:length(f_bands)
    avgRT30_empty(i,:) = mean(RT30_empty(i,:));
    stdRT30_empty(i,:) = std(RT30_empty(i,:));
    avgRT30_abs(i,:) = mean(RT30_abs(i,:));
    stdRT30_abs(i,:) = std(RT30_abs(i,:));
end

% Plots
figure
for i=1:3
    hold on
    plot(f_bands,RT30_empty(:,i))
    plot(f_bands,RT30_abs(:,i))
end
errorbar(f_bands,avgRT30_empty,stdRT30_empty)
errorbar(f_bands,avgRT30_abs,stdRT30_abs)
set(gca, 'XScale', 'log','XMinorTick','off','XMinorGrid','off')
h = get(gca, 'Children');
colorplot=['#c30f0e'; '#181894'];
for i=1:length(h)
    if i>2
        set(h(i), 'Color', '#595959','LineStyle','--');
    else
        set(h(i), 'Color', colorplot(i,:),'LineWidth',1.5);
    end
end
xticks(f_bands)
xlim([100 5000])
title(['Average RT_{30,Empty} = ' num2str(mean(RT30_empty(:))) ' s | Average RT_{30,Abs} =' num2str(mean(RT30_abs(:))) ' s'])
legend('','','','','','Single-Meas','Avg-RT_{30,Empty}','Avg-RT_{30,Abs}',Location='best')

xlabel('Frequency (Hz)')
ylabel('Time (s)')
grid on