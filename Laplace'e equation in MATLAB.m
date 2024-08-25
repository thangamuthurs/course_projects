% FINITE VOLUME METHOD TO SOLVE THE LAPLACE EQUATION ON A POLAR GRID
% THE DISCRETIZED EQUATION IS SOLVED USING SOR (Successive Over-Relaxation)
% Grid Size: 51 X 51

clc;
close all;

% --------------------------
% Parameters Initialization
% --------------------------
imax = 51;
jmax = 51;
nmax = 150;
rwz = 0.1;     % Inner radius
Rx = 0.1;      % Outer radius
rxy = 1.0;     % Maximum radius
Theta_a = 0.0; % Start angle (degrees)
Theta_b = 90.0;% End angle (degrees)
EPS = 1e-5;    % Convergence criterion
OM = 1.5;      % Over-relaxation factor

% --------------------------
% Grid Creation
% --------------------------
imap = imax - 1;
jmap = jmax - 1;

for i = 1:imax
    Theta(i) = deg2rad((i-1)*(Theta_b - Theta_a)/jmap);
    for j = 1:jmax
        R(j) = rwz + (j-1)*(rxy - rwz)/jmap;
        x(i,j) = R(j) * cos(Theta(i));
        y(i,j) = R(j) * sin(Theta(i));
        phix(i,j) = sin(Theta(i)) / R(j);
        phi = phix;
    end
end

% --------------------------
% Boundary Conditions
% --------------------------
for i = 1:imax
    phix(i,1) = 0;                  % Bottom boundary
    phix(i,jmax) = 1 / (rxy - rwz); % Left boundary
end

for j = 1:jmax
    phix(1,j) = sin(Theta(j)) / rwz; % Inner circle boundary
    phix(imax,j) = sin(Theta(j)) / rwz; % Outer circle boundary
end

% --------------------------
% Grid Parameters Calculation
% --------------------------
for i = 2:imap
    for j = 2:jmap
        XA = 0.25 * (x(i,j) + x(i-1,j) + x(i-1,j-1) + x(i,j-1));
        XB = 0.25 * (x(i+1,j) + x(i,j) + x(i,j-1) + x(i+1,j-1));
        XC = 0.25 * (x(i+1,j+1) + x(i,j+1) + x(i,j) + x(i+1,j));
        XD = 0.25 * (x(i,j+1) + x(i-1,j+1) + x(i-1,j) + x(i,j));
        YA = 0.25 * (y(i,j) + y(i-1,j) + y(i-1,j-1) + y(i,j-1));
        YB = 0.25 * (y(i+1,j) + y(i,j) + y(i,j-1) + y(i+1,j-1));
        YC = 0.25 * (y(i+1,j+1) + y(i,j+1) + y(i,j) + y(i+1,j));
        YD = 0.25 * (y(i,j+1) + y(i-1,j+1) + y(i-1,j) + y(i,j));

        % Calculate parameters for each side of the cell
        % Side AB
        DXA = XB - XA;
        DYA = YB - YA;
        DXj = x(i,j) - x(i,j-1);
        DYj = y(i,j) - y(i,j-1);
        SAB = abs(DXA * DYj - DXj * DYA);
        QAB(i,j) = (DXA^2 + DYA^2) / SAB;
        PAB(i,j) = (DXA * DXj + DYA * DYj) / SAB;

        % Side BC
        DXB = XC - XB;
        DYB = YC - YB;
        DXi = x(i,j) - x(i+1,j);
        DYi = y(i,j) - y(i+1,j);
        SBC = abs(DYi * DXB - DXi * DYB);
        QBC(i,j) = (DXB^2 + DYB^2) / SBC;
        PBC(i,j) = (DXB * DXi + DYB * DYi) / SBC;

        % Side CD
        DXC = XD - XC;
        DYC = YD - YC;
        DXj = x(i,j) - x(i,j+1);
        DYj = y(i,j) - y(i,j+1);
        SCD = abs(DXC * DYj - DYC * DXj);
        QCD(i,j) = (DXC^2 + DYC^2) / SCD;
        PCD(i,j) = (DXC * DXj + DYC * DYj) / SCD;

        % Side DA
        DXD = XA - XD;
        DYD = YA - YD;
        DXi = x(i,j) - x(i-1,j);
        DYi = y(i,j) - y(i-1,j);
        SDA = abs(DXi * DYD - DYi * DXD);
        QDA(i,j) = (DXD^2 + DYD^2) / SDA;
        PDA(i,j) = (DXD * DXi + DYD * DYi) / SDA;
    end
end

% --------------------------
% Iteration Using SOR
% --------------------------
for N = 1:nmax
    SUM = 0.0;
    for J = 2:jmap
        JM = J - 1;
        JP = J + 1;
        for I = 2:imap
            IM = I - 1;
            IP = I + 1;
            PHD = 0.25 * (PCD(I,J) - PDA(I,J)) * phi(IM,JP);
            PHD = PHD + (QCD(I,J) + 0.25 * (PBC(I,J) - PDA(I,J))) * phi(I,JP);
            PHD = PHD + 0.25 * (PBC(I,J) - PCD(I,J)) * phi(I,JP);
            PHD = PHD + (QDA(I,J) + 0.25 * (PCD(I,J) - PAB(I,J))) * phi(IM,J);
            PHD = PHD + (QBC(I,J) + 0.25 * (PAB(I,J) - PCD(I,J))) * phi(IP,J);
            PHD = PHD + 0.25 * (PDA(I,J) - PAB(I,J)) * phi(IM,JM);
            PHD = PHD + (QAB(I,J) + 0.25 * (PDA(I,J) - PBC(I,J))) * phi(I,JM);
            PHD = PHD + 0.25 * (PAB(I,J) - PBC(I,J)) * phi(IP,JM);
            PHD = PHD / (QAB(I,J) + QBC(I,J) + QCD(I,J) + QDA(I,J));
            DIF = PHD - phi(I,J);
            SUM = SUM + DIF^2;
            phi(I,J) = phi(I,J) + OM * DIF;
        end
    end
    RMS = sqrt(SUM / (imap-1) / (jmap-1));
    if RMS <= EPS
        fprintf('Converged after %d steps, RMS error: %e\n', N, RMS);
        break;
    end
end

if RMS > EPS
    fprintf('Convergence not achieved in %d steps', N);
end

% --------------------------
% Plotting Results
% --------------------------
contourf(x, y, phi);
colorbar;
xlabel('X', "FontSize", 10);
ylabel('Y', "FontSize", 10);
title('Isolines of phi for 51 X 51 grid', "FontSize", 10);
hold on;

for i = 1:imax
    plot(x(i,:), y(i,:), 'k-', 'LineWidth', 0.01);
end

for j = 1:jmax
    plot(x(:,j), y(:,j), 'k-', 'LineWidth', 0.01);
end
