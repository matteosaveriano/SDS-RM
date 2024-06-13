clear;
close all;

% Add unit quaternion functions, GMR lib, and LASA dataset
addpath(genpath('Lib'))
addpath(genpath('Data'))


%% Note
% I experimentally found out that the EM has poor performance with inverse
% mapping if the direction of source and target trajectories are almost
% orthogonal. This can be seen by using also the first 150 points in the
% demonstration.
%%
%% Learning and reproduction loop
demos_num_loop = 1:26; %[1, 2, 4, 7, 9];
sampling_rate = 10;
demo_iter = 1;
gmm_time = zeros(1,length(demos_num_loop));
stability_thres = 0.1;
for modIt=15%demos_num_loop
    %% Load training data
    % Load demonstrations
    demos = load_LASA_UQ_models('LASA_HandWriting_UQ/', modIt);
    q_goal = array2quat(demos{1}.quat(:,end));
    
%     figure
%     plot(demos{1}.tsPos(1,:), demos{1}.tsPos(2,:))
    
    tic;
    
    % Downsample training data     
    N = size(demos{1}.quat, 2);
    sam_ind = [1:sampling_rate:N N];
    demo_len = length(sam_ind);
    q_ts_demo = []; 
    q_ts_vel = []; 
    q_init_all = [];
    q_ts_base = [];
    q_ts_demo_init_all = [];
    for i=1:length(demos)
        dt_ = sampling_rate*demos{i}.dt;
        
        q_init_all = [q_init_all, demos{i}.quat(:,1)];
        q_demo_array{i} = demos{i}.quat(:, sam_ind);
        q_ts_demo_init_all = [q_ts_demo_init_all, demos{i}.tsPos(:, 1)];
        q_ts_demo = [q_ts_demo, demos{i}.tsPos(:, sam_ind)];
        q_ts_vel = [q_ts_vel [diff(q_ts_demo, [], 2)./dt_, zeros(3,1)]];
        
        % Base stable dynamics
        q_curr = array2quat(demos{i}.quat(:,1));

        % Compute control to have zero error after 'demo_len' steps
        kq = demo_len * dt_ / 5;
        for d=1:demo_len    
            log_q = quat_log(q_curr, q_goal);

            d_log_q = - kq * log_q;

            q_ts_base = [q_ts_base, log_q];

            log_q = log_q + d_log_q * dt_;

            % Convert to Riemmanian Manifold
            tmp = quat_exp(log_q);
            q_curr = quat_mult(tmp, q_goal);
        end
    end
    %% Learn a diffeomorphism between tangent spaces
    D = 3;
    Data = [q_ts_base; q_ts_demo];
    in = 1:D;
    out = D+1:2*D; 

    nbStates = 10;

    [Priors, Mu, Sigma] = EM_init_kmeans(Data, nbStates);
    Mu(:,1) = [mean(q_ts_demo_init_all, 2); mean(q_ts_demo_init_all, 2)];
    Mu(:,end) = Data(:,end);
    [Priors, Mu, Sigma] = EM(Data, Priors, Mu, Sigma);   
    Mu(:,end) = Data(:,end);

    fun = @(pt) GMR(Priors, Mu, Sigma, pt, in, out);
    fun_reverse = @(pt) GMR(Priors, Mu, Sigma, pt, out, in);
    fun_jac = @(pt) GMR_with_Jacobian(Priors, Mu, Sigma, pt, in, out);
    
    gmm_time(demo_iter) = toc;
%     
%     plot(q_ts_base','r')
%     hold on
%     plot(q_ts_demo','b')
%     plot(fun(q_ts_base)','k--')
%     plot(fun_reverse(q_ts_demo)','m--')
% 
% for i=1:demo_len
%     [tmp_, Jac] = fun_jac(q_ts_base(:,i));
%     
%     omega_diff = - kq  * Jac * q_ts_base(:,i);
%     
%     % Integration
%     q_ts_diff(:,i) = tmp_ + dt_ * omega_diff;
% 
%     tmp = quat_exp(q_ts_diff(:,i));
%     q_diff_array(:,i) = quat2array(quat_mult(tmp, q_goal));
% end
% plot(q_diff_array','k--')
% hold on;
% plot(q_demo_array{1}','b')

    % Simulation 
    for d=1:length(demos)
        q_curr = array2quat(q_init_all(:,d));
        q_diff_array = q_init_all(:,d);
        quat_err = quat_log(q_curr, array2quat(demos{d}.quat(:, 1)));
        q_ts_des = demos{d}.tsPos(:, sam_ind);
        
        dist_init = norm(q_ts_des(:,1) - q_ts_des(:,end));
        dist_stable_1 = stability_thres * dist_init;
        dist_stable_2 = dist_stable_1 / 2;
        
        for i=1:demo_len-1
            % Project to tangent space 
            q_ts_curr = quat_log(q_curr, q_goal);
    
            % Apply diffeomorphism and Compute Jacobian
            [q_ts_diff, Jac] = fun_jac(q_ts_curr);
    
            % Stable dynamics
            omega_diff = - kq  * Jac * q_ts_curr;
    
            % Integration
            q_ts_diff = q_ts_diff + dt_ * omega_diff;
    
            % Convert to Riemmanian Manifold
            tmp = quat_exp(q_ts_diff);
            q_diff = quat_mult(tmp, q_goal);
    
            % Invert diffeomorphism
%             q_ts_curr = fun_reverse(q_ts_diff);
    
            % if fun_reverse is not working properly, especially at the 
            % beginning of the trajectory, we can simulate the stable
            % log-dynamics
            % Stable dynamics
            omega = - kq * q_ts_curr;
            % Integration
            q_ts_curr = q_ts_curr + dt_ * omega;
    
            tmp = quat_exp(q_ts_curr);
            q_curr = quat_mult(tmp, q_goal);
            
            % Save data 
            q_diff_array = [q_diff_array, quat2array(q_diff)];
            
            % Errors
            quat_err = [quat_err, q_ts_des(:,i+1)-q_ts_diff];
%             quat_err = [quat_err, quat_log(q_diff, array2quat(demos{d}.quat(:, i+1)))];
        end
        q_diff_save{d} = q_diff_array;
        % Cumulative errors
        quat_err2 = dot(quat_err, quat_err);
        e_quat(demo_iter,d) = sqrt(mean(quat_err2));
        
        % Check convergence
        dist_end = norm(quat_err(:,end));
        if(dist_end < dist_stable_1)
            stable_1(demo_iter,d) = 1;
        else
            stable_1(demo_iter,d) = 0;
        end
        if(dist_end < dist_stable_2)
            stable_2(demo_iter,d) = 1;
        else
            stable_2(demo_iter,d) = 0;
        end
        
        PLOT = 1;
        if(PLOT)
            c_green = [0 127 0]/255;
            c_brown = [212 85 0]/255;
            c_red   = [170 0 0]/255;
            c_blue  = [0 113 188]/255;
            c_yello = [236 176 31]/255;
            time_ = dt_*[0:demo_len-1];
            %figure(modIt)
            h1 = figure(1);
            plot(time_, q_diff_save{d}(1,:), 'color', c_yello, 'LineWidth', 3)
            hold on
            plot(time_, q_diff_save{d}(2,:), 'color', c_blue, 'LineWidth', 3)
            plot(time_, q_diff_save{d}(3,:), 'color', c_green, 'LineWidth', 3)
            plot(time_, q_diff_save{d}(4,:), 'color', c_brown, 'LineWidth', 3)

            plot(time_, q_demo_array{d}', 'k--', 'LineWidth', 3)
        end
        
    end   
    axis([0 3 -0.55 1.1])
    axis square;
    set(gca, 'XTick', [0 1.5 3], 'YTick', [-0.5 0 0.5 1], 'FontName','Times','Fontsize',24);
    xlabel('Time [s]', 'interpreter', 'latex', 'FontName','Times','Fontsize',24);
    ylabel('$\mathbf{\Omega}$', 'interpreter', 'latex', 'FontName','Times','Fontsize',24);
    title('UQ','FontName','Times','Fontsize',24);
    
    demo_iter = demo_iter + 1;
end
% Synthetic results
disp(['UQ GMR K=' num2str(nbStates)])

[R, C] = size(e_quat);
e_quat_vec = reshape(e_quat, R*C, 1);
meanE = num2str(mean(e_quat_vec));
minE  = num2str(min(e_quat_vec));
maxE  = num2str(max(e_quat_vec));
stdPE  = num2str(std(e_quat_vec));
disp(['Position error (mean/range/std) ' meanE ' / [' minE '-' maxE ' / +- ' stdPE ']'])

meanT = num2str(mean(gmm_time));
minT =  num2str(min(gmm_time));
maxT =  num2str(max(gmm_time));
stdT =  num2str(std(gmm_time));
disp(['GMM Training Time (mean/range) ' meanT ' / [' minT '-' maxT ' / +- ' stdT ']'])

[R, C] = size(stable_2);
conv_rate_1 = num2str(sum(sum(stable_1))/(R*C));
conv_rate_2 = num2str(sum(sum(stable_2))/(R*C));
disp(['Convergence rate: 90% - ' conv_rate_1 ' , 95% - ' conv_rate_2 ']'])