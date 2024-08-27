% ForceSimSLE.m
% 2023-8-31

% Fast time-domain simulation which assumes rigid tool/workpiece. Used to
% precalculate forces for estimating cut power, bending moment, SLE, etc

function [Forcex, Forcey, Forcez, torque, times, fluteCuttingWall] = ForceSimSLE(n, b, a, fz, cutDir, cuttingForce, param, options)
    arguments
        n double
        b double
        a double
        fz double
        cutDir int32
        cuttingForce double
        param struct
        options.StepsPerRev double = 2000;
    end

    
    % Load material coefficients from input
    Ktc = cuttingForce(1);
    Krc = cuttingForce(2);
    Kac = cuttingForce(3);
    Kte = cuttingForce(4);
    Kre = cuttingForce(5);
    Kae = cuttingForce(6);
    Cr = 0;
    Ct = 0;

    % Load tool
    m = param.FluteCount;                            % number of teeth, integer
    d = param.ToolDiameter;                              % diameter, m
    if isfield(param, 'ToolCornerRadius')
        cornerRad = param.ToolCornerRadius; % Corner radius of the tool, m
    else
        cornerRad = 0;
    end
    beta = param.HelixAngle*ones(m,1);                        % helix angle vector, deg
    tooth_angle = 0;            % angles of m cutter teeth starting from zero, deg
    tooth_angle = (0:1:m-1)*(360/m);
    % Runout is NOT currently enabled in the tool model
    RO = zeros(size(beta));                         % flute-to-flute runout relative to largest flute, m

    % Load params
    omega = n;                            % spindle speed, rpm
    b = b;                                 % axial depth, m
    ft_mean = fz;                        % mean feed/tooth, m
    [phistart, phiexit] = EntryAngles(param.ToolDiameter, a, cutDir);


    % Account for feed/tooth variation due to non-uniform teeth spacing
    theta = diff([tooth_angle 360]);
    ft = ones(1, m)*ft_mean;
    for cnt = 1:m
        ft(cnt) = (ft_mean*theta(cnt)*m)/360;
    end


    % Simulation specifications
    steps_rev = options.StepsPerRev;              % steps per revolution
    timePerRev = 60/n;                      % Time per revolution (s)
    dt = timePerRev / steps_rev;
    times = dt*(0:steps_rev-1);
    dphi = 360/steps_rev;                   % angular step size between time steps, deg
    if beta == 0
        db = b;
    else
        db = d*(dphi*pi/180)/2/tan(beta(1)*pi/180); % discretized axial depth, m
    end
    steps_axial = round(b/db);              % number of steps along tool axis

    % Initialize vectors
    for cnt = 1:m
        teeth(cnt) = round(tooth_angle(cnt)/dphi) + 1;
    end
    phi = zeros(1, steps_rev);
    for cnt = 1:steps_rev
        phi(cnt) = (cnt - 1)*dphi;
    end

    if cutDir == 1 % Conventional milling
        [~, ind] = min(abs(phi-0));
        phiCutWall = phi(ind);
    else % Climb milling
        [~, ind] = min(abs(phi-180));
        phiCutWall = phi(ind);
    end


    Forcex = zeros(1, steps_rev);
    Forcey = zeros(1, steps_rev);
    Forcez = zeros(1, steps_rev);
    torque = zeros(1, steps_rev);
    fluteCuttingWall = zeros(1, steps_rev);

    % If zero axial depth, then just return zero for everything
    if b == 0
        return
    end

    %************************** MAIN PROGRAM ******************************
    for cnt1 = 1:steps_rev                % time steps, s
        for cnt2 = 1:m
            teeth(cnt2) = teeth(cnt2) + 1;      % index teeth pointer one position (rotate cutter by dphi)
            if teeth(cnt2) > steps_rev
                teeth(cnt2) = 1;
            end
        end

        Fx = 0;
        Fy = 0;
        Fz = 0;

        for cnt3 = 1:m                          % sum forces over all teeth, N
            for cnt4 = 1:steps_axial            % sum forces along axial depth of helical endmill, N
                phi_counter = teeth(cnt3) - (cnt4-1);
                while phi_counter < 1              % helix has wrapped through phi = 0 deg
                    phi_counter = phi_counter + steps_rev;
                end
                phia = phi(phi_counter);        % angle for given axial disk using max helix angle, deg
                phiactual = phi(teeth(cnt3)) - (2*(cnt4-1)*db*tan(beta(m)*pi/180)/d)*180/pi;    % actual angle for selected tooth including local helix lag, deg
                phi_counter_new = round((phiactual-phia)/dphi) + phi_counter;   % counter to select discretized actual phi for selected tooth with local helix, integer
                while phi_counter_new < 1          % helix has wrapped through phi = 0 deg
                    phi_counter_new = phi_counter_new + steps_rev;
                end
                phib = phi(phi_counter_new);    % angle for given axial disk using current helix angle, deg



                if (phib >= phistart) && (phib <= phiexit)                               % verify that tooth angle is in specified range for current disk, deg
                    % Calculate corner radius effects
                    z = cnt4*db;
                    if z < cornerRad % In the corner radius segment
                        kappa_p = acos(1-z/cornerRad);
                        theta = kappa_p - acos(1-(z-db)/cornerRad);
                        dbEff = cornerRad*theta;
                    else % Above the corner radius
                        kappa_p = pi/2;
                        dbEff = db;
                    end
                    

                    h = (ft(cnt3)*sin(phib*pi/180) + RO(cnt3))*sin(kappa_p);   % chip thickness including runout effect, m
                    if h < 0 % tooth jumped out of cut
                        ftan = 0;
                        frad = 0;
                        fax = 0;
                    else    % tooth is engaged in cut
                        if abs(phib - phiCutWall) <= dphi*2 % Flag times where part of the tool is cutting the finished wall
                            fluteCuttingWall(cnt1) = 1;
                        end

                        ftan = Ktc*dbEff*h + Kte*dbEff;
                        frad = Krc*dbEff*h + Kre*dbEff;
                        fax  = Kac*dbEff*h + Kae*dbEff;
                    end
                else    % tooth angle is outside range bounded by radial immersion
                    ftan = 0;
                    frad = 0;
                    fax  = 0;
                    kappa_p = 0;
                end

                cosphib = cosd(phib);
                sinphib = sind(phib);

                Fx = Fx + ftan*cosphib + frad*sinphib*sin(kappa_p) - fax*sinphib*cos(kappa_p);
                Fy = Fy + ftan*sinphib - frad*cosphib*sin(kappa_p) - fax*sinphib*cos(kappa_p);
                Fz = Fz - frad*cos(kappa_p) - fax*sin(kappa_p);

                torque(cnt1) = torque(cnt1) + ftan*d/2;
            end     % cnt4 loop
        end         % cnt3 loop

        Forcex(cnt1) = Fx;
        Forcey(cnt1) = Fy;
        Forcez(cnt1) = Fz;

        % Euler integration for position
    end % cnt1 loop
end


function [phistart, phiexit] = EntryAngles(toolDiameter, radialWidth, feedDirection)
    if feedDirection == 1
        phistart = 0;
        phiexit = acosd((toolDiameter/2-radialWidth)/(toolDiameter/2));
    else
        phistart = 180 - acosd((toolDiameter/2-radialWidth)/(toolDiameter/2));
        phiexit = 180;
    end
end
