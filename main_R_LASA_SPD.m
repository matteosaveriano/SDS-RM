clear;
close all;

% Add unit quaternion functions, GMR lib, and LASA dataset
addpath(genpath('Data'))
addpath(genpath('Lib'))


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
for modIt=13 %demos_num_loop
    %% Load training data
    % Load demonstrations
    demos = load_LASA_SPD_models('LASA_HandWriting_SPD/', modIt);
    S_goal = squeeze(demos{1}.spd(:,:,end));
    
    tic;
    
    % Downsample training data     
    N = size(demos{1}.tsPos, 2);
    sam_ind = [1:sampling_rate:N N];
    demo_len = length(sam_ind);
    s_ts_demo = []; 
    s_ts_vel = []; 
    s_ts_base = [];
    s_ts_demo_init_all = [];
    for d=1:length(demos)
        dt_ = sampling_rate*demos{d}.dt;
        
        s_init_all{d} = squeeze(demos{d}.spd(:,:,1));
        s_ts_demo_init_all = [s_ts_demo_init_all, demos{d}.tsPos(:, 1)];
        s_ts_demo = [s_ts_demo, demos{d}.tsPos(:, sam_ind)];
        s_ts_vel = [s_ts_vel [diff(s_ts_demo, [], 2)./dt_, zeros(3,1)]];
            
        % Base stable dynamics
        S_curr = squeeze(demos{d}.spd(:,:,1));

        % Compute control to have zero error after 'demo_len' steps
        kq = demo_len * dt_ / 5;
        for i=1:demo_len   
            s_demo_array{d}(:,i) = symMat2Vec(demos{d}.spd(:, :, sam_ind(i)));
            
            log_s = symMat2Vec(logmap(S_goal, S_curr));
    
            s_dot = -kq * log_s;
    
            % Store data for training
            s_ts_base = [s_ts_base, log_s];
    
            % Integration
            log_s = log_s + s_dot * dt_;
    
            S_curr = expmap(S_goal, vec2symMat(log_s));
        end
    end
    %% Learn a diffeomorphism between tangent spaces
    D = 3;
    Data = [s_ts_base; s_ts_demo];
    scale_ = 100;
    Data = Data / scale_;
    in = 1:D;
    out = D+1:2*D; 

    nbStates = 10;

    [Priors, Mu, Sigma] = EM_init_kmeans(Data, nbStates);
    %Mu(:,1) = [mean(s_ts_demo_init_all, 2); mean(s_ts_demo_init_all, 2)];
    Mu(:,end) = Data(:,end);
    [Priors, Mu, Sigma] = EM(Data, Priors, Mu, Sigma);   
    Mu(:,end) = Data(:,end);

    fun = @(pt) GMR(Priors, Mu, Sigma, pt, in, out);
    fun_reverse = @(pt) GMR(Priors, Mu, Sigma, pt, out, in);
    fun_jac = @(pt) GMR_with_Jacobian(Priors, Mu, Sigma, pt, in, out);
    
    gmm_time(demo_iter) = toc;
%     
%     plot(s_ts_base','r')
%     hold on
%     plot(s_ts_demo','b')
%     plot(fun(s_ts_base)','k--')
%     plot(fun_reverse(s_ts_demo)','m--')
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
        S_curr = s_init_all{d};
        s_diff_array = symMat2Vec(s_init_all{d});
        spd_err = symMat2Vec(logmap(demos{d}.spd(:,:,1),S_curr));
        s_ts_des = demos{d}.tsPos(:, sam_ind);
        
        dist_init = norm(s_ts_des(:,1) - s_ts_des(:,end));
        dist_stable_1 = stability_thres * dist_init;
        dist_stable_2 = dist_stable_1 / 2;
        
        for i=1:demo_len-1
            % Project to tangent space 
            s_ts_curr =  symMat2Vec(logmap(S_goal, S_curr));
      
            % Apply diffeomorphism and Compute Jacobian
            [s_ts_diff, Jac] = fun_jac(s_ts_curr / scale_);

            % Compute Diffeomorphed dynamics
            sd_ts_diff = - kq * Jac * ( s_ts_curr / scale_);
    
            % Intergration and projection to Riemmanian manifold
            s_ts_diff = scale_*(s_ts_diff + dt_ * sd_ts_diff);
    
            % Convert to Riemmanian Manifold
            S_diff = expmap(S_goal, vec2symMat(s_ts_diff));
    
            % Base dynamics
            sd_ts_curr = - kq * s_ts_curr;
            % Intergration and projection to Riemmanian manifold
            s_ts_curr = s_ts_curr + sd_ts_curr * dt_;
    
            % Convert to Riemmanian Manifold
            S_curr = expmap(S_goal, vec2symMat(s_ts_curr));
            
            % Save data 
            s_diff_array = [s_diff_array, symMat2Vec(S_diff)];
            
            % Errors
            spd_err = [spd_err, s_ts_des(:,i+1)-s_ts_diff];
        end
        s_diff_save{d} = s_diff_array;
        % Cumulative errors
        spd_err2 = dot(spd_err, spd_err);
        e_spd(demo_iter,d) = sqrt(mean(spd_err2));
        
        % Check convergence
        dist_end = norm(spd_err(:,end));
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
            
            h1 = figure(modIt);
            hold on
            plot(time_, s_diff_save{d}(1,:), 'color', c_brown, 'LineWidth', 3)
            plot(time_, s_diff_save{d}(2,:), 'color', c_blue, 'LineWidth', 3)
            plot(time_, s_diff_save{d}(3,:), 'color', c_green, 'LineWidth', 3)
            
            plot(time_, s_demo_array{d}', 'k--', 'LineWidth', 3)
        end
    end   
    axis([0 3 -60 180])
    axis square;
    set(gca, 'XTick', [0 1.5 3], 'YTick', [-60 60 180], 'FontName','Times','Fontsize',24);
    xlabel('Time [s]', 'interpreter', 'latex', 'FontName','Times','Fontsize',24);
    ylabel('$\mathbf{\Upsilon}$', 'interpreter', 'latex', 'FontName','Times','Fontsize',24);
    title('SPD','FontName','Times','Fontsize',24);
    
    demo_iter = demo_iter + 1;
end
% Synthetic results
disp(['SPD GMR K=' num2str(nbStates)])

[R, C] = size(e_spd);
e_spd_vec = reshape(e_spd, R*C, 1);
meanE = num2str(mean(e_spd_vec));
minE  = num2str(min(e_spd_vec));
maxE  = num2str(max(e_spd_vec));
stdPE  = num2str(std(e_spd_vec));
disp(['Position error (mean/range/std) ' meanE ' / [' minE '-' maxE ' / +- ' stdPE ']'])

meanT = num2str(mean(gmm_time));
minT =  num2str(min(gmm_time));
maxT =  num2str(max(gmm_time));
stdT =  num2str(std(gmm_time));
disp(['GMM Training Time error (mean/range) ' meanT ' / [' minT '-' maxT ' / +- ' stdT ']'])

[R, C] = size(stable_2);
conv_rate_1 = num2str(sum(sum(stable_1))/(R*C));
conv_rate_2 = num2str(sum(sum(stable_2))/(R*C));
disp(['Convergence rate: 90% - ' conv_rate_1 ' , 95% - ' conv_rate_2 ']'])