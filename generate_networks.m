% Generates synthetic networks (MATLAB version)

function generate_networks(new_dir, dir_lfr)
    nb_layers = 2;

    if ~exist(dir_lfr, 'dir')
        mkdir(dir_lfr);
    end
    if ~exist(new_dir, 'dir')
        mkdir(new_dir);
    end
    if ~exist([new_dir '_networks'], 'dir')
        mkdir([new_dir '_networks']);
    end

    % Parameters
    n = 100;
    k = 6;
    maxk = 10;
    mu = 0.05;

    % Config
    list_alpha = [0.2, 0.4, 0.6, 0.8];
    list_mu = mu;
    list_p  = [0.01, 0.05, 0.1, 0.15];
    list_p1 = [0.1, 0.15, 0.2];
    list_p2 = [0.0, 0.1];

    intra_params.n    = n;
    intra_params.k    = k;
    intra_params.maxk = maxk;
    intra_params.mu   = mu;

    % Call create_layers (must be implemented in MATLAB separately)
    create_layers(dir_lfr, nb_layers, intra_params);

    % Iteration
    nb_iteration = length(list_alpha) * length(list_mu) * ...
                   length(list_p) * length(list_p1) * length(list_p2);

    i = 0;
    for alpha = list_alpha
        new_dir_alpha = fullfile(new_dir, ['alpha-' num2str(alpha)]);
        if ~exist(new_dir_alpha, 'dir')
            mkdir(new_dir_alpha);
        end

        for p = list_p
            new_dir_p = fullfile(new_dir_alpha, ['p-' num2str(p)]);
            if ~exist(new_dir_p, 'dir')
                mkdir(new_dir_p);
            end

            for mu_val = list_mu
                new_dir_mu = fullfile(new_dir_p, ['mu-' num2str(mu_val)]);
                if ~exist(new_dir_mu, 'dir')
                    mkdir(new_dir_mu);
                end

                for p1 = list_p1
                    new_dir_p1 = fullfile(new_dir_mu, ['p1-' num2str(p1)]);
                    if ~exist(new_dir_p1, 'dir')
                        mkdir(new_dir_p1);
                    end

                    for p2 = list_p2
                        new_dir_p2 = fullfile(new_dir_p1, ['p2-' num2str(p2)]);
                        if ~exist(new_dir_p2, 'dir')
                            mkdir(new_dir_p2);
                        end

                        i = i + 1;
                        fprintf('%d/%d\n', i, nb_iteration);

                        % Copy layer directories
                        copyfile(fullfile(dir_lfr, 'layer0'), new_dir_p2);
                        copyfile(fullfile(dir_lfr, 'layer1'), new_dir_p2);

                        % Run the Python script
                        cmd = sprintf('python2.7 lfr_multilayer_v3.py %s %s -p %f -a %f -p1 %f -p2 %f', ...
                                      new_dir_p2, new_dir_p2, p, alpha, p1, p2);
                        disp(cmd);
                        system(cmd);

                        % Copy output file if exists
                        srcFile = fullfile(new_dir_p2, 'new_format');
                        if isfile(srcFile)
                            destFile = sprintf('%s_networks/network_%g_%g_%g_%g_%g', ...
                                new_dir, alpha, p, mu_val, p1, p2);
                            copyfile(srcFile, destFile);
                            disp('File copied ----------');
                        end

                        pause(0.15);
                    end
                end
            end
        end
    end
end

