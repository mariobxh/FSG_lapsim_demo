function Track = Generate_LinCurv_turn(step_length,straight,turn_length,min_radius)

max_curv = 1/min_radius;

curv_straight = zeros(1,ceil(straight/step_length));
curv_ascent = linspace(0,max_curv,ceil(turn_length/(2*step_length)));
curv_descent = linspace(max_curv,0,ceil(turn_length/(2*step_length)));

curv = [curv_straight, curv_ascent, curv_descent(2:end), curv_straight];

dS = step_length * ones(size(curv));

psi = cumtrapz(dS.*curv);

Track.x = cumtrapz(dS.*cos(psi));
Track.y = cumtrapz(dS.*sin(psi));
Track.psi = psi;
Track.curv = curv;
Track.S = cumsum(dS);
Track.N = length(Track.x);

end