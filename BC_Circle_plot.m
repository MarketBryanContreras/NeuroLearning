function [Xcord,Ycord]=BC_Circle_plot(rad, xCent, yCent, color, plot_flag)
%% The function BC_Circle_plot:
%            [Plot a circle with the radious, x center and  y center specified by the user]
%   Inputs:
%            -rad[Int]: radious
%            -xCent[Int]: X coordinate of the center
%            -yCent[Int]: y coordinate of the center

%
% 
% 
%   Outputs:
%             -output1[type]: description 
% 
% 
%  First version BC 22-Jan-2024 
%% 
theta =0:0.1:2*pi;
Xcord=xCent + rad*cos(theta);
Ycord=yCent + rad*sin(theta);

%Plot the coordinates
if plot_flag==1
    if nargin >=3;
        plot(Xcord,Ycord, 'Color',color, 'LineWidth', 2);
    else
        plot(Xcord,Ycord, 'LineWidth', 2);
    end
end