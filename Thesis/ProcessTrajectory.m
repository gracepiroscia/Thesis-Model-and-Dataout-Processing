function [A_DATA, V_DATA,T_DATA, D_DATA] = ProcessTrajectory(acc_data, v_data, t_data, d_data)
%PROCESSTRAJECTORY processes raw velocity-time data from the EV buggy,
%ensuring that the trajectory includes braking trajectory only (i.e.
%initial velocity is the highest value and final is the lowest) and smooths
%the curve to prepare for fitting:

    %find initial and final braking velocities
    [~, start_idx] = max(v_data); %find highest and lowest velocities 
    [~, end_idx] = min(v_data);
    
    % Cut trajectory to include braking range only
    v_data = v_data(start_idx:end_idx); 
    t_data = t_data(start_idx:end_idx);
    d_data = d_data(start_idx:end_idx);
    
    % Use moving mean to smooth
    v_data = movmean(v_data,8); %8 data point window for smoothing
    acc_data = movmean(acc_data,8);
    
    A_DATA = acc_data;
    V_DATA = v_data;
    T_DATA = t_data;
    D_DATA = d_data;
    
    
end

