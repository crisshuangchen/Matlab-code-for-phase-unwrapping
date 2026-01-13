% clear;
% close;
% clc;


%% Mueller matrices in RGB channels, the data comes from real measurements of HWP.
M_R = [1.0000   -0.0132   -0.0019    0.0010
       0.0263   -0.8244   -0.2969   -0.4610
       0.0054   -0.3046    0.9509   -0.0851
       -0.0063    0.4664    0.0715   -0.8786];
M_G = [1.0000   -0.0055    0.0022    0.0029
       0.0195   -0.8762   -0.2956    0.0608
       0.0193   -0.3167    0.9856   -0.0056
       -0.0013   -0.0543   -0.0203   -0.9864];
M_B = [1.0000   -0.0102   -0.0061    0.0042
      -0.0007   -0.5646   -0.2487    0.5656
       0.0135   -0.2716    0.9851    0.0832
       0.0043   -0.6361   -0.1129   -0.7065];

%% Get wrapped_theta and wrapped_delta
[wrapped_theta_R, wrapped_delta_R] = cal_theta_and_delta(M_R);
[wrapped_theta_G, wrapped_delta_G] = cal_theta_and_delta(M_G);
[wrapped_theta_B, wrapped_delta_B] = cal_theta_and_delta(M_B);

%% Effective wavelngths
lamda_R = 629.9;
lamda_G = 543.94;
lamda_B = 455.18;

%% Get approximation value
angdiff = @(a,b) abs(mod(a - b + 90, 180) - 90); %this takes into account that 90 and -90 are the same angle
theta_difference_BG = angdiff(wrapped_theta_B, wrapped_theta_G);
theta_difference_RG = angdiff(wrapped_theta_R, wrapped_theta_G);

if (wrapped_delta_G + wrapped_delta_B + wrapped_delta_R) < 0.2
    approximation_delta_R = wrapped_delta_R;
    approximation_delta_G = wrapped_delta_G;
    approximation_delta_B = wrapped_delta_B;
elseif   abs(theta_difference_BG + theta_difference_RG) < 55
    %% RGB
    Delta_n_x_d_1 = abs((wrapped_delta_B-wrapped_delta_G)*(lamda_B*lamda_G/(lamda_B-lamda_G))/(2*pi));
    Delta_n_x_d_2 = abs((wrapped_delta_R-wrapped_delta_G)*(lamda_R*lamda_G/(lamda_R-lamda_G))/(2*pi));
    Delta_n_x_d_3 = abs((wrapped_delta_R-wrapped_delta_B)*(lamda_R*lamda_B/(lamda_R-lamda_B))/(2*pi));
    Delta_n_x_d = (Delta_n_x_d_1 + Delta_n_x_d_2 + Delta_n_x_d_3)/3;

    approximation_delta_R = 2*pi*Delta_n_x_d/lamda_R;
    approximation_delta_G = 2*pi*Delta_n_x_d/lamda_G;
    approximation_delta_B = 2*pi*Delta_n_x_d/lamda_B;
elseif theta_difference_BG < 45
     %% BG
     Delta_n_x_d = abs((wrapped_delta_B-wrapped_delta_G)*(lamda_B*lamda_G/(lamda_B-lamda_G))/(2*pi));

     approximation_delta_R = 2*pi*Delta_n_x_d/lamda_R;
     approximation_delta_G = 2*pi*Delta_n_x_d/lamda_G;
     approximation_delta_B = 2*pi*Delta_n_x_d/lamda_B;

elseif  theta_difference_RG < 45
     %% RG
     Delta_n_x_d = abs((wrapped_delta_R-wrapped_delta_G)*(lamda_R*lamda_G/(lamda_R-lamda_G))/(2*pi));

     approximation_delta_R = 2*pi*Delta_n_x_d/lamda_R;
     approximation_delta_G = 2*pi*Delta_n_x_d/lamda_G;
     approximation_delta_B = 2*pi*Delta_n_x_d/lamda_B;

else %it means there is jump of order between colors, we then we used G\R
     Delta_n_x_d = abs((wrapped_delta_R-(2*pi-wrapped_delta_G))*(lamda_R*lamda_G/(lamda_R-lamda_G))/(2*pi));

     approximation_delta_R = 2*pi*Delta_n_x_d/lamda_R;
     approximation_delta_G = 2*pi*Delta_n_x_d/lamda_G;
     approximation_delta_B = 2*pi*Delta_n_x_d/lamda_B;
end


%% Get true value
k_B = sign(approximation_delta_B/(2*pi) - round(approximation_delta_B/(2*pi)));
true_delta_B = round(approximation_delta_B/(2*pi))*2*pi + k_B*wrapped_delta_B;

k_G = sign(approximation_delta_G/(2*pi) - round(approximation_delta_G/(2*pi)));
true_delta_G = round(approximation_delta_G/(2*pi))*2*pi + k_G*wrapped_delta_G;

k_R = sign(approximation_delta_R/(2*pi) - round(approximation_delta_R/(2*pi)));
true_delta_R = round(approximation_delta_R/(2*pi))*2*pi + k_R*wrapped_delta_R;

%% Print
fprintf('Wrapped delta (radians): R=%.4f, G=%.4f, B=%.4f\n', wrapped_delta_R, wrapped_delta_G, wrapped_delta_B);
fprintf('Wrapped theta (degree): R=%.4f, G=%.4f, B=%.4f\n', wrapped_theta_R, wrapped_theta_G, wrapped_theta_B);

fprintf('Approximation delta (radians): R=%.4f, G=%.4f, B=%.4f\n', approximation_delta_R, approximation_delta_G, approximation_delta_B);
fprintf('True delta (radians): R=%.4f, G=%.4f, B=%.4f\n', true_delta_R, true_delta_G, true_delta_B);


%% Function 
% Get wrapped theta and wrapped delta using MM
function [theta, delta] = cal_theta_and_delta(M)

m00 = M(1,1); m01 = M(1,2); m02 = M(1,3); m03 = M(1,4);
m10 = M(2,1); m11 = M(2,2); m12 = M(2,3); m13 = M(2,4);
m20 = M(3,1); m21 = M(3,2); m22 = M(3,3); m23 = M(3,4);
m30 = M(4,1); m31 = M(4,2); m32 = M(4,3); m33 = M(4,4);

C = 0.25*[m00+m11+m22+m33      m01+m10-1i*(m23-m32)  m02+m20+1i*(m13-m31)  m03+m30-1i*(m12-m21);
         m01+m10+1i*(m23-m32)  m00+m11-m22-m33       m12+m21+1i*(m03-m30)  m13+m31-1i*(m02-m20);
         m02+m20-1i*(m13-m31)  m12+m21-1i*(m03-m30)  m00-m11+m22-m33       m23+m32+1i*(m01-m10);
         m03+m30+1i*(m12-m21)  m13+m31+1i*(m02-m20)  m23+m32-1i*(m01-m10)  m00-m11-m22+m33];

[eig_vector, eig_value] = eig(C);
[~, idx] = max(real(diag(eig_value)));
c = eig_vector(:, idx);

J = [c(1)+c(2)    c(3)-1i*c(4);
     c(3)+1i*c(4) c(1)-c(2)]; 

phi0 = angle(J(1,1));
J = J * exp(-1i*phi0);
if real(J(1,1)) < 0
    J = -J;
end

detJ = det(J);
K = 1 / sqrt(detJ);     
s = real(K *(J(1,1) + J(2,2)) / 2);   
s = max(min(s, 1), -1);              
T = 2 * acos(s);  

den = 2*sin(T / 2);
if abs(den) < 0.0001
    womiga = K;
else
    womiga = K * (T) / den;
end


L = 1i*womiga*(J(1,1)-J(2,2));
L_ = 1i*womiga*(J(1,2)+J(2,1));

LB  = real(L);
LB_ = real(L_);


delta = sqrt(LB^2 + LB_^2);   %0~pi

if delta>pi
   delta = 2*pi -delta;
end

if delta<0.01
    delta = 0;
    theta = -90 + (90 - (-90)) * rand(1);  
else
    theta = rad2deg(atan2(LB_, LB) / 2);  %-90~90
end

end




