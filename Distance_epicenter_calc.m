% Calculation of source-to-site epicentral distances

% This function has been vectorized, allowing it to calculate distances
% to multiple points in a single function call.

% Variable xy contains the coordinates of a point: a vector of length 2.
% Variables x and y are vectors [x1, x2, ... , xn] and [y1, y2, ... , yn]
% containing coordinates for points for which distances are required:
% [x1, y1], [x2, y2], ... , [xn, yn].

function distance = Distance_epicenter_calc(xy, x, y)

R=6371; % Earth Radius

    for i=1:length(x)

        t1=deg2rad(xy(2));
        t2=deg2rad(y(i));
        dt=deg2rad(y(i)-xy(2));
        dl=deg2rad(x(i)-xy(1));

        a=((sin(dt/2))^2) + cos(t1)*cos(t2)*((sin(dl/2))^2);
        c=2*atan2(sqrt(a),sqrt(1-a));
        distance(i)=int16(R*c);

    end 

end