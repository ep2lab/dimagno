function do_exit = out_of_box_2d(data,front,i)
% returns 1 if all points in front are out of the box; 0 otherwise
do_exit = 1;
for j = 1:F.n_points
    if front.z_(j) >= data.exit.parameters.Xlim(1) && ...
       front.z_(j) <= data.exit.parameters.Xlim(2) && ...
       front.r_(j) >= data.exit.parameters.Ylim(1) && ...
       front.r_(j) <= data.exit.parameters.Ylim(2)
       do_exit = 0; % there are still points in the box
       return
    end
end
