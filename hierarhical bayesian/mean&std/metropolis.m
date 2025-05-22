function [samples, acceptance_rate] = metropolis(log_target_pdf, n_samples, initial_state, proposal_std, n_dim, type)
 
    % Initialize
    samples = zeros(n_dim, n_samples); % Array to store samples
    samples(:,1) = initial_state;    % Set the initial state
    n_accepted = 0;                % Counter for accepted proposals
    
    % Metropolis algorithm loop
    for i = 2:n_samples
        % Propose a new state
        proposal = samples(:, i - 1) + proposal_std * randn(n_dim, 1);
        
        % Compute acceptance probability
        if type == "analytical" %#ok<IFBDUP>
            log_alpha = log_target_pdf(proposal) - log_target_pdf(samples(:, i - 1));
        else
            % figure()
            % plot(log_target_pdf(1,:), log_target_pdf(2,:));
            term1 = interp1(log_target_pdf(1,:), log_target_pdf(2,:), proposal, 'linear', "extrap" );
            term2 = interp1(log_target_pdf(1,:), log_target_pdf(2,:), samples(i - 1), 'linear', "extrap" );
            log_alpha = term1 - term2;
        end
        alpha = min(1, exp(log_alpha));
        
        % Accept or reject the proposal
        if rand < alpha
            samples(:, i) = proposal; % Accept
            n_accepted = n_accepted + 1; % Increment acceptance counter
        else
            samples(:, i) = samples(:, i - 1); % Reject
        end
    end
    
    % Compute the acceptance rate
    acceptance_rate = n_accepted / (n_samples - 1); % Exclude the initial state
end