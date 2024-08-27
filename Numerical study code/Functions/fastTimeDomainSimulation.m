function results = fastTimeDomainSimulation(rpm, ap, ae, fd, fz, tool, material, dynamics, param)    
    [kt, kn, kte, kne, ct, cn] = material.MaterialMatrix;
    
    [kx, mx, cx, wnx, x_modes] = dynamics.XParameters;
    [ky, my, cy, wny, y_modes] = dynamics.YParameters;
    
    % Load tool
    Nt = tool.FluteCount;
    d = tool.Diameter;                      % teeth diameter, m
    gamma = tool.HelixAngle;                     % helix angle, deg

    % Load process
    if fd == 1                          % up milling
        phistart = 0;                   % starting angle, deg
        phiexit = acos((d/2 - ae)/(d/2)); % exit angle, deg
    else                                % down milling
        phistart = pi - acos((d/2 - ae)/(d/2));
        phiexit = pi;
    end

    omega = rpm;                  % rpm
    b = ap;                   % m
    ft = fz;                   % m
	
	% Set simulation length
    if isfield(param, 'tdsRevolutions')
		rev = param.tdsRevolutions;
	else
		rev = 100;
    end
	
    if isfield(param, 'tdsTimeStepFactor')
        timeStepFactor = param.tdsTimeStepFactor;
    else
        timeStepFactor = 20;
    end
    fs = max([wnx;wny]/2/pi, [], 'all') * timeStepFactor; %  * 20;
    dt = 1/fs;
    steps_rev = round(60/(dt*omega)) - mod(round(60/(dt*omega)),Nt);
    % Round the sample frequency to get an even number of samples per rev
    dt = 60 / omega / steps_rev;
    fs = 1 / dt;
    dphi = 2*pi/steps_rev;           % deg
    V = pi*d*omega/60; % Radial velocity
    

    % Determine axial depth step, m
    if gamma == 0
        db = b;
    else
        db = d*(dphi)/2/tan(gamma);
		minDb = (max(param.Aps) - min(param.Aps))/(length(param.Aps)-1);
		if db > minDb % Enforce that the axial step has to be less than the simulation resolution
			steps_rev = ceil((steps_rev * db / minDb) / Nt) * Nt;
																											 
            dt = 60 / omega / steps_rev;
			dphi = 2*pi/steps_rev;           % deg
			db = d*(dphi)/2/tan(gamma);
		end
    end
    
	% If zero axial depth, then skip simulation and return 0
    if b == 0
        results.xDisplacement = zeros(rev*steps_rev,1);
        results.yDisplacement = zeros(rev*steps_rev,1);
        results.xForce = zeros(rev*steps_rev,1);
        results.yForce = zeros(rev*steps_rev,1);
        results.time = reshape(linspace(0,(rev*steps_rev)*dt-dt,rev*steps_rev), size(zeros(rev*steps_rev,1)));
        results.sampleFrequency = fs;
        return
    end

    % number of steps along tool axis
    steps_axial = ceil(b/db);
    lastStepDb = b - db * steps_axial;
    if lastStepDb == 0, lastStepDb = db; end
    steps = rev*steps_rev;
    
    time = linspace(0,steps*dt-dt,steps);
    
    % Initialize vectors
    if ~isfield(param, 'FluteSpacing')
        teeth = 1: steps_rev/Nt :steps_rev;
    else
        toothAbsoluteAngles = cumsum(param.FluteSpacing)-param.FluteSpacing(1);
        teeth = round(toothAbsoluteAngles / (2*pi) * steps_rev) + 1;
    end
    
    phi = (0:dphi:steps_rev*dphi-dphi)';
    phi_counter = zeros(steps_axial, Nt);
    phi_map = zeros(steps_axial, Nt, steps_rev);
    
    for i = 1:steps_rev
        phi_counter(1:end,:) = teeth - (1:steps_axial)';
        phi_counter(phi_counter<1) = phi_counter(phi_counter<1) + steps_rev;
    
        phi_map(:,:,i) = phi_counter;
        teeth = teeth+1;
        teeth(teeth>steps_rev) = 1;
    end
    
    phi_map_full_sin = sin(phi);
    phi_map_full_cos = cos(phi);
    
    phi_is_cutting_index = zeros(size(phi));
    phi_is_cutting_index(phi>=phistart & phi<=phiexit) = true;
        
    xpos = zeros(steps,1);
    ypos = zeros(steps,1);
    Forcex = zeros(steps,1);
    Forcey = zeros(steps,1);
    torque = zeros(steps,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Test Phi Map
        
    surf = zeros(steps_axial, steps_rev);

    % Euler integration initial conditions
    x = 0;
    dx = 0;
    y = 0;
    dy = 0;
    dp = zeros(1, x_modes);
    p = zeros(1, x_modes);          % x-direction modal displacements, m
    dq = zeros(1, y_modes);
    q = zeros(1, y_modes);          % y-direction modal displacements, m
    
    engaged = zeros(size(xpos));
    for cnt0 = 1:rev
        for cnt1 = 1:steps_rev

            Fx = 0;
            Fy = 0;

            for cnt3 = 1:Nt
                for cnt4 = 1:steps_axial

                    phi_counter = phi_map(cnt4, cnt3, cnt1);

                    if phi_is_cutting_index(phi_counter)

                        phi_sin = phi_map_full_sin(phi_counter);
                        phi_cos = phi_map_full_cos(phi_counter);

                        n = x*phi_sin - y*phi_cos;            % m
                        h = ft*phi_sin + surf(cnt4, phi_counter) - n;  % m

                        if h < 0
                            Ft = 0;
                            Fn = 0;
                            surf(cnt4, phi_counter) = surf(cnt4, phi_counter) + ft*phi_sin;
                        else
                            r_dot = (dx * phi_sin) - (dy*phi_cos);
                            if cnt4 ~= steps_axial
                                Ft = kt*db*h + kte*db - ct*db*r_dot/V;
                                Fn = kn*db*h + kne*db - cn*db*r_dot/V;
                                engaged(cnt1 + (steps_rev*(cnt0-1))) = engaged(cnt1 + (steps_rev*(cnt0-1))) + 1;
                            else
                                Ft = kt*lastStepDb*h + kte*lastStepDb - ct*lastStepDb*r_dot/V;
                                Fn = kn*lastStepDb*h + kne*lastStepDb - cn*lastStepDb*r_dot/V;
                                engaged(cnt1 + (steps_rev*(cnt0-1))) =  engaged(cnt1 + (steps_rev*(cnt0-1))) + 1;
                            end
                            surf(cnt4, phi_counter) = n;

                            Fx = Fx + Ft*phi_cos + Fn*phi_sin;
                            Fy = Fy + Ft*phi_sin - Fn*phi_cos;
                            torque(cnt1 + (steps_rev*(cnt0-1))) = torque(cnt1 + (steps_rev*(cnt0-1))) + Ft*d/2;
            
                        end
                    end
                end
            end

            % Numerical integration for position
            x = 0;
            dx = 0;
            y = 0;
            dy = 0;

            % x direction
            for cnt5 = 1:x_modes
                ddp = (Fx - cx(cnt5)*dp(cnt5) - kx(cnt5)*p(cnt5))/mx(cnt5);
                dp(cnt5) = dp(cnt5) + ddp*dt;
                dx = dx + dp(cnt5);
                p(cnt5) = p(cnt5) + dp(cnt5)*dt;
                x = x + p(cnt5);        % m
            end

            % y direction
            for cnt5 = 1:y_modes
                ddq = (Fy - cy(cnt5)*dq(cnt5) - ky(cnt5)*q(cnt5))/my(cnt5);
                dq(cnt5) = dq(cnt5) + ddq*dt;
                dy = dy + dq(cnt5);
                q(cnt5) = q(cnt5) + dq(cnt5)*dt;
                y = y + q(cnt5);        % m
            end

            xpos(cnt1 + (steps_rev*(cnt0-1))) = x;
            ypos(cnt1 + (steps_rev*(cnt0-1))) = y;
            Forcex(cnt1 + (steps_rev*(cnt0-1))) = Fx;
            Forcey(cnt1 + (steps_rev*(cnt0-1))) = Fy;

        end
    end
	
	if isfield(param, 'tdsTrimAmount')
		trimFactor = param.tdsTrimAmount;
	else
		trimFactor = 0.5;
    end

	startInd = round(length(xpos) * trimFactor);

    results.xDisplacement = xpos(startInd:end);
    results.yDisplacement = ypos(startInd:end);
    results.xForce = Forcex(startInd:end);
    results.yForce = Forcey(startInd:end);
    results.time = reshape(time(startInd:end), size(results.xForce));
    results.sampleFrequency = fs;
    results.torque = torque(startInd:end);
    results.engagedFrac = mean(engaged(startInd:end));
end
