
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Authors: Khanh Dao Duc, Jalal Khouhak
%         Department of Mathematics 
%         The University of British Columbia (UBC)
% 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = IVT_model_fit_complete(scale_limits, time_points, protein_level, concentrations, ...
  preincubation_times, n, param_selection, Nb_step, neld_param, TIME, deltat, file_name, loss_scale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inference Procedure Function %%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% Parsing and formatting input %%%%%%%%%
parameter_lists = ["\alpha_P", "\alpha_M", "\alpha_R", ...
  "\delta_P", "\delta_M", "\delta_R","\beta", ...
  "\gamma", ];
parameter_lists_string = ["a_P", "a_M", "a_R", ...
  "d_P", "d_M","d_R", ...
  "betta","gamma"];
active_parameters = parameter_lists(param_selection);
active_parameters_string = parameter_lists_string(param_selection);

preincubation_times_set = unique(preincubation_times);

disp("Starting IVT model fit with following parameters turned on:");
disp(active_parameters_string);

%%%%%%%%% Initial parameter scale finding %%%%%%%%%
scale = zeros(1,n);
scale_perm_number = zeros(1,n);
for i = 1:n
  M{i} = scale_limits(i,1):scale_limits(i,2);
  scale_perm_number(i) = length(M{i});
end

bwd_prod = cumprod(flip(scale_perm_number));
bwd_prod = flip([1, bwd_prod(1:n-1)]);
scale_perm_tot = prod(scale_perm_number);

for i = 1:n
  scale(i) = M{i}(length(M{i}));
end
ref = loss(IVT_output(ones(1,n), 10 .^ scale),ones(1,n));
scale_temp = zeros(1,n);

for ns = 1:scale_perm_tot-1
  scale_carry_over = ns;
  for i = 1:n-1
    scale_idx_temp = floor(scale_carry_over / bwd_prod(i));
    scale_temp(i) = M{i}(scale_idx_temp + 1);
    scale_carry_over = scale_carry_over - scale_idx_temp * bwd_prod(i);
  end
  if scale_carry_over == 0
      scale_carry_over = length(M{n});
  end
  scale_temp(n) = M{n}(scale_carry_over);
  
  ref_min = loss(IVT_output(ones(1,n), 10 .^ scale_temp),ones(1,n));
  if ref_min < ref
    ref = ref_min;
    scale = scale_temp;
  end
end 

set_scale = 10 .^ scale;
disp("Scale found:");
% set_scale=50*set_scale;
disp(set_scale);


fprintf('Starting Nelder-Mead method with %d steps. \n', Nb_step);
%%%%%%%%% Nelder-Mead Method %%%%%%%%%
a = neld_param(1);
b = neld_param(2);
c = neld_param(3);
d = neld_param(4);

% set initial points x1,...,xn+1.
x(1:n,1:n+1) = 0.0;
  for i = 1 : n
    x(i,i) = sqrt ( 1.0 - sum ( x(1:i-1,i).^2 ) ); %  Set X(I,I) so that sum ( X(1:I,I)^2 ) = 1.
    for j = i + 1 : n + 1
      x(i,j) = ( - 1.0 / n - ( x(1:i-1,i)' * x(1:i-1,j) ) ) / x(i,i); %  Set X(I,J) for J = I+1 to N+1 by using the fact that XI dot XJ = - 1 / N
    end
  end
X = ones(n+1,n)+transpose(x);
X_out = zeros(1,n+1);
for k=1:n+1
    X_out(1,k)= loss(IVT_output(X(k,:),set_scale),X(k,:));
end

error=zeros(1,Nb_step);
X_star=[];

for step=1:Nb_step
  if mod(step, floor(Nb_step/10)) == 0
    fprintf('Progress: %d%% (step %d) \n', floor(step * 100/Nb_step), step);
    disp(set_scale .* X(1,:));
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [~, B]= sort(X_out);
  X = X(B,:); % order points such that loss(x1)?loss(x2)???loss(xn+1).
  X_out = X_out(B);
  x_o = zeros(1,n);

  for i=1:n
      x_o(i)=mean(X(1:n,i));%compute xo as the centroid of x1,?,xn.
  end  

  x_r = x_o + a*(x_o - X(n+1,:));  %Reflection
  x_r_out=loss(IVT_output(x_r,set_scale),x_r);

  if (X_out(1)< x_r_out) && (X_out(n)> x_r_out)
     X(n+1,:)=x_r;
      X_out(n+1) = x_r_out;
  else      
     if X_out(1)> x_r_out %expansion
        x_e = x_o + c*(x_o - X(n+1,:));
        x_e_out=loss(IVT_output(x_e,set_scale),x_e);
        if x_e_out< x_r_out 
          X(n+1,:)=x_e;
          X_out(n+1) = x_e_out;
        else
          X(n+1,:)=x_r;
          X_out(n+1) = x_r_out;  
        end
     else
       x_c = x_o + b*(x_o - X(n+1,:)); %contraction
        x_c_out=loss(IVT_output(x_c,set_scale),x_c);
        if x_c_out< X_out(n+1)
          X(n+1,:)=x_c;
          X_out(n+1) = x_c_out;
        else
          for j=1:n+1  
          X(j,:)=X(1,:)+d*(X(j,:) - X(1,:));
          X_out(j) = loss(IVT_output(X(j,:),set_scale),X(j,:)) ;
          end
        end
     end       
  end

  error(step) = min(X_out);
  X_star = [X_star; X(1,:)];
end

%%%%%%%%% Figure Generation %%%%%%%%%
fig = figure(1);
fig.Units = 'pixels';
fig.Position = [0, 0, 1000, 300*n];
hold on;
subplot(n+1,1,1)
plot(error)
title_string = [];
for it = 1:n
    title_string = [title_string, active_parameters(it), " = ", num2str(X(1,it)*set_scale(it)), ", "];
end
title(strjoin(title_string));
ylabel('Error') 
for i = 1:n
    subplot(n+1,1,i+1)
    plot(set_scale(i) * X_star(:,i))
    ylabel(active_parameters(i))

end
saveas(fig, sprintf('%s-parameter_convergence.jpg',file_name));

fig = figure(2);
prec_times_len = length(preincubation_times_set);
fig.Units = 'pixels';
fig.Position = [0, 0, 400*prec_times_len, 600];
hold on;
for inc = 1:prec_times_len
  preincubation_idx = find(preincubation_times == preincubation_times_set(inc));
%   subplot(7,prec_times_len,inc)
    subplot(4,prec_times_len,inc)

  title_string = ["Preincubation ", preincubation_times_set(inc), "min"];
  title(strjoin(title_string));
  hold on;
  xlabel('time (min)') 
  ylabel('Protein (a.u.)') 
  for s = preincubation_idx
    plot(time_points + preincubation_times(s), protein_level(s,:),'ko--')
  end

%   subplot(7,prec_times_len,inc+prec_times_len)
  subplot(4,prec_times_len,inc+prec_times_len)

  xlabel('time (min)') 
  ylabel('M_0')
  %%%%%%%%%%%%%%
%   subplot(7,prec_times_len,inc+2*prec_times_len)
  subplot(4,prec_times_len,inc+2*prec_times_len)

  xlabel('time (min)') 
  ylabel('M_1')
  %%%%%%%%%%%%%%%%%%%%%%%%
%   subplot(7,prec_times_len,inc+3*prec_times_len)
%   xlabel('time (min)') 
%   ylabel('M_0 + M_1 ')
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   %%%%%%%%%%%%%%%%%%%
%   subplot(7,prec_times_len,inc+4*prec_times_len)
%   xlabel('time (min)') 
%   ylabel('M_1/M_0 ')
%   %%%%%%%%%%%%%%%%%%%
%   subplot(7,prec_times_len,inc+5*prec_times_len)
%   xlabel('time (min)') 
%   ylabel('M_1/(M_0 + M1) ')
  %%%%%%%%%%%%%%%%%%%%%%%%
%   subplot(7,prec_times_len,inc+6*prec_times_len)
  subplot(4,prec_times_len,inc+3*prec_times_len)

  xlabel('time (min)') 
  ylabel('R_{e}')
  for s = preincubation_idx
      IVT_plot(X(1,:),concentrations(s),set_scale, preincubation_times(s), inc, prec_times_len);
  end
end
saveas(fig, sprintf('%s-IVT-fits-Plasmo.jpg',file_name));

params = set_scale .* X(1,:);
disp("Final Parameter Values:");
disp(active_parameters_string);
disp(params);

disp("Final Error:");
disp(error(end))

%%%%%%%%% Function Definitions %%%%%%%%%

%%%%%%%%% IVT Model %%%%%%%%%
    function dy = model_IVT(y, x, incubation_start)
        % Variable key
        % y(1) = P - protein level
        % y(2) = M_i - mRNA level
        % y(3) = M0 - 
        % y(4) = E - resource level
        
        % x(1): alpha_P - production rate of protein
        % x(2): alpha_M - production cost of mRNA
        % x(3): alpha_E - production cost of resources
        % x(4): delta_P - degradation rate of protein
        % x(5): delta_M - degradation rate of mRNA
        % x(6): delta_E - degradation rate of resources
        % x(7): r_P - recycle of protein into resource
        % x(8): r_M - recycle of mRNA into resource
        
        dy = zeros(4,1);            

        dy(1) = 0;
        dy(2) = 0; 
        dy(3) = 0;
        dy(4) = 0;
%         
        if incubation_start
            
            dy(1) = dy(1) + max([0,x(1)*y(3)*y(4)]) - max([0,x(4)*y(1)]); 
            
            dy(2) = dy(2) - max([0,x(7)*y(2)]);
         
            dy(3) = dy(3) + max([0,x(7)*y(2)]) - max([0,x(8)*y(3)]);
       
            dy(4) = dy(4) - max([0,(x(3))*y(3)*y(4)]);
        else
            y(2)=0;
            
            dy(1) = dy(1) + max([0,x(1)*y(3)*y(4)]) - max([0,x(4)*y(1)]); 
            
            dy(2) = 0;

            dy(3) = dy(3) + max([0,x(7)*y(2)]) - max([0,x(8)*y(3)]);

            dy(4) = dy(4) - max([0,((x(3)))*y(3)*y(4)]) ;
            
        end
    end

%%%%%%%%% Objective Function %%%%%%%%%
    function ans2 = loss(z,Z)
%         ans2 = sqrt(sum(sum((z-protein_level).^2))) + sum(-min(0,Z)) * 10^6;
        ans2 = sqrt(dot(sum((z-protein_level).^2), loss_scale)) + sum(-min(0,Z)) * 10^6; % time reverse importance
%         ans2 = sqrt(sum(sum(((z-protein_level)./ protein_level).^2 ))) + sum(-min(0,Z)) * 10^6; % time reverse importance
    end

%%%%%%%%% Simulated Data %%%%%%%%%
    function y = IVT_output(x, scale_factor)
        y = zeros(length(concentrations),length(time_points));
        for ik=1:length(concentrations)
            y(ik,:)= model_IVT_sim(concentrations(ik), x, scale_factor, preincubation_times(ik));
        end
    end

%%%%%%%%% Model simulation %%%%%%%%%
    function Y = runge_kutta(concentration, param, scale_factor, preincubation_time)
      % Runge-Kutta 4-point method
      % Initial Conditions
%       Z(1,1) = 0;
%       Z(2,1) = concentration;
%       Z(3,1) = 1;
      Z(1,1) = 0;
      Z(2,1) = concentration;
      Z(3,1) = 0;
      Z(4,1) = 1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Y = Z;
      full_time = TIME + preincubation_time;
      incubation_start = false;
      if preincubation_time == 0
        incubation_start = true;
      end

      % Simulation of the model (Runge Kutta 4 )
      N = floor(full_time/deltat) + 1;
      for i=1:N
          Time(i,1) = deltat*(i-1);
      end

      param_scaled = scale_factor .* param; % multiply scaling to all active parameters
      param_full = zeros(1, length(param_selection)); % multiply zero to all the inactive parameters
      param_full(param_selection) = param_scaled;

      for i=1:N-1 
        if (i * deltat) > preincubation_time
          incubation_start = true;
        end

        y = Z;
        k1 = model_IVT(y, param_full, incubation_start);
        
        yprime = y + deltat*k1/2;
        k2 = model_IVT(yprime, param_full, incubation_start);
        
        yprime = y + deltat*k2/2;
        k3 = model_IVT(yprime, param_full, incubation_start);
        
        yprime = y + deltat*k3/2;  
        k4 = model_IVT(yprime, param_full, incubation_start);

        Z = Z + (1/6) * deltat * (k1 + 2*k2 + 2*k3 + k4);

        Y(1,i+1) = Z(1);
        Y(2,i+1) = Z(2);
        Y(3,i+1) = Z(3);
        Y(4,i+1) = Z(4);
      end
    end

%%%%%%%%% Model output %%%%%%%%%
    function ANS = model_IVT_sim(concentration, param, scale_factor, preincubation_time)
      Y = runge_kutta(concentration, param, scale_factor, preincubation_time);
      ANS = Y(1,floor(time_points/deltat));
    end



%%%%%%%%% Model plot %%%%%%%%%
    function IVT_plot(param, c, scale_factor, preincubation_time, plot_idx, plot_cols)
      Y = runge_kutta(c, param, scale_factor, preincubation_time);
      full_time = TIME + preincubation_time;
      N = floor(full_time/deltat)+1;
      for i=1:N
        Time(i,1) = deltat*(i-1);
      end

      %%%%%%%%%%%%%%%%%%% plots/figures/save %%%%%%%%%%%%%%%%%%%%% 
%       subplot(7,plot_cols,plot_idx)
      subplot(4,plot_cols,plot_idx)

      hold on
      plot(Time,Y(1,:),'LineWidth',1);

%       subplot(7,plot_cols,plot_idx+plot_cols)
      subplot(4,plot_cols,plot_idx+plot_cols)

      hold on
      plot(Time,Y(2,:),'LineWidth',1);

      subplot(4,plot_cols,plot_idx+2*plot_cols);
%       subplot(7,plot_cols,plot_idx+2*plot_cols);

      hold on
      plot(Time,Y(3,:),'LineWidth',1);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       subplot(7,plot_cols,plot_idx+3*plot_cols);
%       hold on
%       M0_M1 = Y(3,:)+ Y(2,:);
%       plot(Time,M0_M1,'LineWidth',2);
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       subplot(7,plot_cols,plot_idx+4*plot_cols);
%       hold on
%       plot(Time,Y(3,:)./Y(2,:),'LineWidth',2);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       subplot(7,plot_cols,plot_idx+5*plot_cols);
%       hold on
%       plot(Time,(Y(3,:)./(Y(2,:)+Y(3,:))),'LineWidth',2);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
      subplot(4,plot_cols,plot_idx+3*plot_cols);
      hold on
      plot(Time,Y(4,:),'LineWidth',1);
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
end