function h=LTplotIvBars(Iv,plot_data,color,alpha)
max_data = max(plot_data); 
min_data = min(plot_data); 
offset = (max_data - min_data);%I use these lines to calculate the y limits (+- 10% of these value) 
range_final= [min_data-offset*0.1 max_data+offset*0.1];
hold on
for ii = 1:length(Iv.tstart)
    
    h=fill([Iv.tstart(ii) Iv.tend(ii) Iv.tend(ii) Iv.tstart(ii)], [range_final(1) range_final(1) range_final(2) range_final(2)],color, 'FaceAlpha', alpha, 'LineStyle',"none");
   
end

end
% The function fill creates rentangles using the arguments in the following fashion:fill([xmin xmax xmax xmin], [ymin ymin ymax ymax],[color], 'FaceAlpha', 0.3, 'LineStyle',"none")

