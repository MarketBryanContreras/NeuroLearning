function h=LTplotIvBarsScalogram(Iv,color,alpha)
range_final= [1 128];
hold on
for ii = 1:length(Iv.tstart)
    
    h=fill([Iv.tstart(ii) Iv.tend(ii) Iv.tend(ii) Iv.tstart(ii)], [range_final(1) range_final(1) range_final(2) range_final(2)],color, 'FaceAlpha', alpha, 'LineStyle',"none");
   
end

end
% The function fill creates rentangles using the arguments in the following fashion:fill([xmin xmax xmax xmin], [ymin ymin ymax ymax],[color], 'FaceAlpha', 0.3, 'LineStyle',"none")

