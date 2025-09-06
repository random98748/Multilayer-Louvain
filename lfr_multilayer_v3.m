function lfr_multilayer_v3(outdir, dir_lfr, varargin)
% MATLAB port of lfr_multilayer_v3.py (condensed single file).
% Usage:
%   lfr_multilayer_v3(outdir, dir_lfr, '-p', 0.1, '-a', 0.6, '-p1', 0.15, '-p2', 0.1)
%
% Notes:
% - Expects LFR-generated layers already present in dir_lfr/layer0 and dir_lfr/layer1,
%   each with network.dat (tab-separated u v) and community.dat (tab-separated node communityId).
% - If you need to *create* layers, call create_layers() (provided below) first.
%
% Files written:
%   outdir/couplings
%   outdir/new_format

    if nargin < 2
        error('Usage: lfr_multilayer_v3(outdir, dir_lfr, ["-p", p, "-a", a, "-p1", p1, "-p2", p2, "-d", d])');
    end

    % ---- Parse args (INTER params only, as per the Python) ----
    p  = 0; a  = 0; p1 = 0; p2 = 0; d = 0; %#ok<NASGU>
    i = 1;
    while i <= numel(varargin)
        key = varargin{i};
        if i == numel(varargin)
            error('Missing value for argument %s', key);
        end
        val = varargin{i+1};
        switch key
            case '-p',  p  = double(val);
            case '-a',  a  = double(val);
            case '-p1', p1 = double(val);
            case '-p2', p2 = double(val);
            case '-d',  d  = double(val);  %#ok<NASGU>
            otherwise
                error('Unknown parameter: %s', key);
        end
        i = i + 2;
    end

    nb_layers = 2;

    % ---- Load layers (already generated) ----
    layers = load_layers(dir_lfr); % struct array 1..nb_layers

    % ---- Create cross-layer communities ----
    [cross_communities, intra_coms, nb_intra_coms] = ...
        create_cross_layer_communities(nb_layers, layers, a);
    nb_com = size(cross_communities,1) + nb_intra_coms;

    fprintf('Setting network communities (alpha = %f)\n', a);
    fprintf('%d communities created - %d single layer - %d cross-layer\n', ...
        nb_com, nb_intra_coms, size(cross_communities,1));

    % ---- Build couplings ----
    couplings_links = init_couplings(nb_layers);
    nb_in = apply_p1_parameter(nb_layers, layers, couplings_links, cross_communities, p1);
    nb_out = apply_p2_parameter(nb_layers, layers, couplings_links, intra_coms, cross_communities, p2);
    [nb_in, nb_out] = apply_p_parameter(nb_layers, layers, couplings_links, ...
                                        intra_coms, cross_communities, nb_in, nb_out, p);
    nb_couplings = nb_in + nb_out;

    % ---- Write couplings file ----
    if ~exist(outdir, 'dir'); mkdir(outdir); end
    fid = fopen(fullfile(outdir, 'couplings'), 'w');
    for l = 1:nb_layers-1
        pairs = couplings_links(l).to(l+1).pairs;
        for r = 1:size(pairs,1)
            fprintf(fid, '%s %s\n', pairs{r,1}, pairs{r,2});
        end
    end
    fclose(fid);

    % ---- Write final format ----
    write_soum_format(dir_lfr, outdir, nb_layers, layers, nb_couplings, ...
                      couplings_links, nb_com, cross_communities, intra_coms);

    % ---- Stats ----
    disp('Statistics :');
    for iL = 1:nb_layers
        fprintf('Layer %d\n', iL-1);
        fprintf('  - Nb of links %d\n', layers(iL).nb_links);
    end
    for iL = 1:nb_layers-1
        fprintf('Cross-layer %d-%d\n', iL-1, iL);
        nlinks = size(couplings_links(iL).to(iL+1).pairs,1);
        fprintf('  - Nb of coupling links : %d\n', nlinks);
        fprintf('  - Nb inside coms : %d\n', nb_in);
        fprintf('  - Nb outside coms : %d\n', nb_out);
    end
end

% ======================================================================
% ======================== Helper subfunctions ==========================
% ======================================================================

function layers = load_layers(dir_lfr)
% Loads layers info from dir_lfr/layer0 and dir_lfr/layer1 (community.dat).
% Also computes offsets (from_id_node, from_id_com) and renumbers as in Python.

    sep = sprintf('\t');
    nb_layers = 2;
    layers = repmat(struct( ...
        'com_nodes', containers.Map('KeyType','int32','ValueType','any'), ...
        'node_com',  containers.Map('KeyType','char','ValueType','int32'), ...
        'nodes',     {{}}, ...
        'nb_coms',   0, ...
        'nb_links',  0, ...
        'nb_nodes',  0, ...
        'from_id_node', 0, ...
        'from_id_com',  0), 1, nb_layers);

    % Load layer 0
    [nn0, nl0, nc0, com_nodes0, node_com0, nodes0] = load_com_nodes(0, dir_lfr, sep);
    layers(1).from_id_node = 0;
    layers(1).from_id_com  = 0;
    layers(1).nb_coms = nc0;
    layers(1).nb_links = nl0;
    layers(1).nb_nodes = nn0;
    layers(1).com_nodes = com_nodes0;
    layers(1).node_com  = node_com0;
    layers(1).nodes     = nodes0;

    % Load layer 1
    [nn1, nl1, nc1, com_nodes1, node_com1, nodes1] = load_com_nodes(1, dir_lfr, sep);
    layers(2).from_id_node = layers(1).nb_nodes;
    layers(2).from_id_com  = layers(1).nb_coms;
    layers(2).nb_coms = nc1;
    layers(2).nb_links = nl1;
    layers(2).nb_nodes = nn1;
    layers(2).com_nodes = com_nodes1;
    layers(2).node_com  = node_com1;
    layers(2).nodes     = nodes1;
end

function [nb_nodes, nb_link, nb_coms, com_nodes, node_com, nodes] = load_com_nodes(num, directory, sep)
% Reads community.dat in layer <num>, builds node->com and com->nodes maps.
    dir_layer = fullfile(directory, sprintf('layer%d', num));
    filename_coms = fullfile(dir_layer, 'community.dat');

    com_nodes = containers.Map('KeyType','int32','ValueType','any'); % comId(int32) -> cellstr nodes
    node_com  = containers.Map('KeyType','char','ValueType','int32'); % node(str) -> comId(int32)
    nodes = {};

    nb_link = 0;
    fid = fopen(filename_coms, 'r');
    if fid < 0
        error('Cannot open %s', filename_coms);
    end
    C = textscan(fid, '%s%s', 'Delimiter', sep, 'CollectOutput', true);
    fclose(fid);

    if isempty(C{1})
        nb_nodes = 0; nb_link = 0; nb_coms = 0; return;
    end

    pairs = C{1}; % Nx2 cell array of strings {node, com}
    nb_link = size(pairs,1);
    nodes = pairs(:,1);

    for r = 1:size(pairs,1)
        nodeStr = pairs{r,1};
        comRaw  = pairs{r,2};
        % Python uses: index_com = int(line[1]) - 1
        idx = int32(str2double(comRaw) - 1);
        if ~isKey(com_nodes, idx)
            com_nodes(idx) = {nodeStr};
        else
            com_nodes(idx) = [com_nodes(idx), {nodeStr}];
        end
        node_com(nodeStr) = idx;
    end

    nb_nodes = numel(nodes);
    nb_coms  = numel(com_nodes.keys);

    % Also compute nb_links from network.dat line count
    filename_net = fullfile(dir_layer, 'network.dat');
    nb_link = count_lines(filename_net);
end

function n = count_lines(fname)
    n = 0;
    fid = fopen(fname, 'r');
    if fid < 0
        n = 0; return;
    end
    while true
        tline = fgetl(fid);
        if ~ischar(tline); break; end
        n = n + 1;
    end
    fclose(fid);
end

function layers = create_layers(outdir, nb_layers, dict_intra)
% (Optional) creates layers using external ./benchmark and renumbers IDs.
% dict_intra: struct with fields n, k, maxk, mu
    layers = repmat(struct('nb_coms',0,'nb_links',0,'nb_nodes',0, ...
                           'from_id_node',0,'from_id_com',0), 1, nb_layers);
    for i = 1:nb_layers
        [nb_com_layer, nb_link_layer] = benchmark_lfr(i-1, outdir, dict_intra.n, dict_intra.k, dict_intra.maxk, dict_intra.mu);
        layers(i).nb_coms  = nb_com_layer;
        layers(i).nb_links = nb_link_layer;
        layers(i).nb_nodes = dict_intra.n;
        fprintf('Layer %d created (with %d intra layer coms)\n', i-1, nb_com_layer);
    end

    index_node = 0; index_com = 0;
    for i = 1:nb_layers
        rewrite_nodes_and_coms(i-1, outdir, index_node, index_com);
        layers(i).from_id_node = index_node;
        layers(i).from_id_com  = index_com;
        index_node = index_node + layers(i).nb_nodes;
        index_com  = index_com  + layers(i).nb_coms;
    end
end

function [nb_com_layer, nb_link_layer] = benchmark_lfr(num, directory, n, k, maxk, mu)
% Calls external LFR generator and moves *.dat into layer dir.
    dir_layer = fullfile(directory, sprintf('layer%d', num));
    if ~exist(dir_layer, 'dir'); mkdir(dir_layer); end

    cmd = sprintf('./benchmark -N %d -k %g -maxk %d -mu %g', n, k, maxk, mu);
    disp(cmd);
    system(cmd);

    % move *.dat into dir_layer
    movefile('*.dat', dir_layer, 'f');

    % Count unique communities and edge lines
    comm = fullfile(dir_layer, 'community.dat');
    net  = fullfile(dir_layer, 'network.dat');

    if ~isfile(comm) || ~isfile(net)
        error('benchmark did not produce expected files in %s', dir_layer);
    end

    % unique communities in column 2
    fid = fopen(comm,'r');
    C = textscan(fid, '%s%s', 'Delimiter', sprintf('\t'), 'CollectOutput', true);
    fclose(fid);
    if isempty(C{1}); nb_com_layer = 0; else
        u = unique(C{1}(:,2));
        nb_com_layer = numel(u);
    end

    nb_link_layer = count_lines(net);
end

function rewrite_nodes_and_coms(num, directory, from_node, from_com)
% Adds offsets to node IDs and com IDs (like Python rewrite).
    sep = sprintf('\t');
    dir_layer = fullfile(directory, sprintf('layer%d', num));
    filename_nodes = fullfile(dir_layer, 'network.dat');
    filename_coms  = fullfile(dir_layer, 'community.dat');

    % network.dat
    fin = fopen(filename_nodes, 'r'); fout = fopen([filename_nodes 'bis'], 'w');
    C = textscan(fin, '%s%s', 'Delimiter', sep, 'CollectOutput', true);
    fclose(fin);
    data = C{1};
    for r = 1:size(data,1)
        a = num2str(str2double(data{r,1}) + from_node);
        b = num2str(str2double(data{r,2}) + from_node);
        fprintf(fout, '%s%s%s\n', a, sep, b);
    end
    fclose(fout);
    delete(filename_nodes);
    movefile([filename_nodes 'bis'], filename_nodes, 'f');

    % community.dat
    fin = fopen(filename_coms, 'r'); fout = fopen([filename_coms 'bis'], 'w');
    C = textscan(fin, '%s%s', 'Delimiter', sep, 'CollectOutput', true);
    fclose(fin);
    data = C{1};
    for r = 1:size(data,1)
        a = num2str(str2double(data{r,1}) + from_node);
        b = num2str(str2double(data{r,2}) + from_com);
        fprintf(fout, '%s%s%s\n', a, sep, b);
    end
    fclose(fout);
    delete(filename_coms);
    movefile([filename_coms 'bis'], filename_coms, 'f');
end

function [cross_communities, intra_coms, nb_intra_coms] = create_cross_layer_communities(nb_layers, layers, alpha)
% Returns:
%   cross_communities: Mx2 numeric (for 2 layers) of global community IDs
%   intra_coms: struct with fields 1..nb_layers -> numeric vector of com IDs (global)
%   nb_intra_coms: count
    cross_communities = zeros(0,2);
    intra_coms = struct;
    nb_intra_coms = 0;

    if alpha == 0
        for i = 1:nb_layers
            % global com ids in layer i
            ids = layers(i).from_id_com + (0:layers(i).nb_coms-1);
            intra_coms.(sprintf('L%d', i)) = ids;
            nb_intra_coms = nb_intra_coms + numel(ids);
        end
        return;
    end

    % Pair communities between consecutive layers
    for i = 1:nb_layers-1
        ids1 = layers(i).from_id_com   + (0:layers(i).nb_coms-1);
        ids2 = layers(i+1).from_id_com + (0:layers(i+1).nb_coms-1);

        nb_cross = ceil(min(numel(ids1), numel(ids2)) * alpha);

        [ids1_shuf, ids2_shuf] = deal(ids1(randperm(numel(ids1))), ids2(randperm(numel(ids2))));
        cross_pairs = [ids1_shuf(1:nb_cross).', ids2_shuf(1:nb_cross).'];

        cross_communities = [cross_communities; cross_pairs]; %#ok<AGROW>

        intra1 = setdiff(ids1, cross_pairs(:,1));
        intra2 = setdiff(ids2, cross_pairs(:,2));
        intra_coms.(sprintf('L%d', i))   = intra1;
        intra_coms.(sprintf('L%d', i+1)) = intra2;
        nb_intra_coms = nb_intra_coms + numel(intra1) + numel(intra2);
    end
end

function couplings_links = init_couplings(nb_layers)
% couplings_links(l).to(m).pairs: cell array N x 2 of strings (node ids)
    couplings_links = repmat(struct('to',[]), 1, nb_layers);
    for l = 1:nb_layers
        couplings_links(l).to = repmat(struct('pairs',{{}}), 1, nb_layers);
    end
end

function nb_link_inside = apply_p1_parameter(nb_layers, layers, couplings_links, cross_communities, p1)
% Connects nodes *within* cross-layer communities
    if isempty(cross_communities); nb_link_inside = 0; return; end
    nb_link_inside = 0;

    for i = 1:nb_layers-1
        % ensure list allocated
        if isempty(couplings_links(i).to(i+1).pairs)
            couplings_links(i).to(i+1).pairs = cell(0,2);
        end

        ids1 = layers(i).from_id_com   + (0:layers(i).nb_coms-1);
        ids2 = layers(i+1).from_id_com + (0:layers(i+1).nb_coms-1);

        % cross pairs for these two layers
        mask = ismember(cross_communities(:,1), ids1) & ismember(cross_communities(:,2), ids2);
        cross_l1_l2 = cross_communities(mask,:);

        for r = 1:size(cross_l1_l2,1)
            id1 = int32(cross_l1_l2(r,1));
            id2 = int32(cross_l1_l2(r,2));

            list1 = get_nodes_from_com(layers(i).com_nodes, id1);
            list2 = get_nodes_from_com(layers(i+1).com_nodes, id2);

            n1 = round(numel(list1) * p1);
            n2 = round(numel(list2) * p1);
            s1 = sample_cells(list1, n1);
            s2 = sample_cells(list2, n2);

            S1 = s1; S2 = s2;
            while ~isempty(S1) || ~isempty(S2)
                nmin = min(numel(S1), numel(S2));
                for j = 1:nmin
                    i1 = randi(numel(S1));
                    i2 = randi(numel(S2));
                    node1 = S1{i1}; node2 = S2{i2};

                    couplings_links(i).to(i+1).pairs(end+1,:) = {node1, node2}; %#ok<AGROW>
                    nb_link_inside = nb_link_inside + 1;

                    S1(i1) = []; S2(i2) = [];
                end
                if ~isempty(S1)
                    S2 = sample_cells(s2, min(numel(S1), numel(s2)));
                elseif ~isempty(S2)
                    S1 = sample_cells(s1, min(numel(S2), numel(s1)));
                end
            end
        end
    end
end

function nb_link_outside = apply_p2_parameter(nb_layers, layers, couplings_links, intra_coms, cross_communities, p2)
% Connects nodes *outside* cross-layer communities (from intra-layer communities)
    nb_link_outside = 0;
    if isempty(fieldnames(intra_coms)); return; end

    for i = 1:nb_layers-1
        if isempty(couplings_links(i).to(i+1).pairs)
            couplings_links(i).to(i+1).pairs = cell(0,2);
        end

        ids1 = layers(i).from_id_com   + (0:layers(i).nb_coms-1);
        ids2 = layers(i+1).from_id_com + (0:layers(i+1).nb_coms-1);
        mask = ismember(cross_communities(:,1), ids1) & ismember(cross_communities(:,2), ids2);
        cross_l1_l2 = cross_communities(mask,:);

        list_l1 = {};
        list_l2 = {};

        % intra nodes
        intra1 = intra_coms.(sprintf('L%d', i));
        intra2 = intra_coms.(sprintf('L%d', i+1));
        for id1 = intra1
            list_l1 = [list_l1, get_nodes_from_com(layers(i).com_nodes, int32(id1))]; %#ok<AGROW>
        end
        for id2 = intra2
            list_l2 = [list_l2, get_nodes_from_com(layers(i+1).com_nodes, int32(id2))]; %#ok<AGROW>
        end
        % plus nodes from cross communities (both sides)
        for r = 1:size(cross_l1_l2,1)
            list_l1 = [list_l1, get_nodes_from_com(layers(i).com_nodes, int32(cross_l1_l2(r,1)))]; %#ok<AGROW>
            list_l2 = [list_l2, get_nodes_from_com(layers(i+1).com_nodes, int32(cross_l1_l2(r,2)))]; %#ok<AGROW>
        end

        % sample from intra communities in layer i to any node in layer i+1
        for id1 = intra1
            nodes1 = get_nodes_from_com(layers(i).com_nodes, int32(id1));
            n1 = round(numel(nodes1) * p2);
            s1 = sample_cells(nodes1, n1);
            for j = 1:numel(s1)
                n2 = randi(numel(list_l2));
                couplings_links(i).to(i+1).pairs(end+1,:) = {s1{j}, list_l2{n2}}; %#ok<AGROW>
                nb_link_outside = nb_link_outside + 1;
            end
        end

        % sample from intra communities in layer i+1 to any node in layer i
        for id2 = intra2
            nodes2 = get_nodes_from_com(layers(i+1).com_nodes, int32(id2));
            n2 = round(numel(nodes2) * p2);
            s2 = sample_cells(nodes2, n2);
            for j = 1:numel(s2)
                % avoid duplicates if possible
                attempt = 0; max_attempt = 100;
                while true
                    n1 = randi(numel(list_l1));
                    cand = {list_l1{n1}, s2{j}};
                    if ~pair_exists(couplings_links(i).to(i+1).pairs, cand) || attempt > max_attempt
                        break;
                    end
                    attempt = attempt + 1;
                end
                couplings_links(i).to(i+1).pairs(end+1,:) = cand; %#ok<AGROW>
                nb_link_outside = nb_link_outside + 1;
            end
        end
    end
end

function tf = pair_exists(pairs, cand)
    tf = false;
    for r = 1:size(pairs,1)
        if strcmp(pairs{r,1}, cand{1}) && strcmp(pairs{r,2}, cand{2})
            tf = true; return;
        end
    end
end

function [nb_in, nb_out] = apply_p_parameter(nb_layers, layers, couplings_links, intra_coms, cross_communities, nb_in, nb_out, p)
% Adjusts the fraction p (inside cross-layer communities vs outside)
    if isempty(cross_communities); return; end

    list_com_node = struct; % key as 'l1_l2' -> struct(layerIndex -> set of nodes)
    for i = 1:nb_layers-1
        ids1 = layers(i).from_id_com   + (0:layers(i).nb_coms-1);
        ids2 = layers(i+1).from_id_com + (0:layers(i+1).nb_coms-1);
        mask = ismember(cross_communities(:,1), ids1) & ismember(cross_communities(:,2), ids2);
        cross_l1_l2 = cross_communities(mask,:);

        pairs = couplings_links(i).to(i+1).pairs;
        for r = 1:size(cross_l1_l2,1)
            c1 = cross_l1_l2(r,1);
            c2 = cross_l1_l2(r,2);
            key = sprintf('%d_%d', c1, c2);
            list_com_node.(key) = containers.Map('KeyType','int32','ValueType','any');
            list_com_node.(key)(int32(i-1)) = unique_nodes_in_pairs(pairs, get_nodes_from_com(layers(i).com_nodes, int32(c1)), 1);
            list_com_node.(key)(int32(i))   = unique_nodes_in_pairs(pairs, get_nodes_from_com(layers(i+1).com_nodes, int32(c2)), 2);
        end
    end

    p_temp = (1.0 * nb_in) / max(1, (nb_in + nb_out));
    fprintf('actual p value : %f\n', p_temp);

    nb_new = 0;
    if p > p_temp
        % add inside links between already-used nodes within same cross-pair
        total_possible = 0;
        keys = fieldnames(list_com_node);
        for k = 1:numel(keys)
            M = list_com_node.(keys{k});
            for i = 1:(length(M.keys)-1)
                nA = numel(M(int32(i-1))); nB = numel(M(int32(i)));
                total_possible = total_possible + nA * nB;
            end
        end
        total_possible = total_possible - nb_in;

        while p > p_temp && nb_new < total_possible
            k = randi(numel(keys));
            iLayer = 1; % only 2 layers: i-1 and i in Python terms
            M = list_com_node.(keys{k});
            A = M(int32(iLayer-1));
            B = M(int32(iLayer));

            if isempty(A) || isempty(B), break; end
            node1 = A{randi(numel(A))};
            node2 = B{randi(numel(B))};

            pairs = couplings_links(iLayer).to(iLayer+1).pairs;
            if ~pair_exists(pairs, {node1,node2})
                couplings_links(iLayer).to(iLayer+1).pairs(end+1,:) = {node1,node2};
                nb_new = nb_new + 1;
                p_temp = (1.0 * (nb_in + nb_new)) / max(1, (nb_in + nb_out));
            end
        end
        nb_in = nb_in + nb_new;
        fprintf('Nb of coupling links created (inside cross-layer communities) : %d\n', nb_new);

    elseif p < p_temp
        % add outside links by mixing nodes from different cross-pairs
        keys = fieldnames(list_com_node);
        total_possible = 0;
        for i1 = 1:numel(keys)
            for i2 = 1:numel(keys)
                if i1 == i2, continue; end
                A = list_com_node.(keys{i1})(int32(0));
                B = list_com_node.(keys{i2})(int32(1));
                total_possible = total_possible + numel(A)*numel(B);
            end
        end

        nb_new = 0;
        while p < p_temp && nb_new < total_possible
            if numel(keys) < 2, break; end
            % choose two distinct cross-pairs
            idx = randperm(numel(keys), 2);
            k1 = keys{idx(1)}; k2 = keys{idx(2)};
            A = list_com_node.(k1)(int32(0));
            B = list_com_node.(k2)(int32(1));
            if isempty(A) || isempty(B), break; end

            iLayer = 1; node1 = A{randi(numel(A))}; node2 = B{randi(numel(B))};
            pairs = couplings_links(iLayer).to(iLayer+1).pairs;
            if ~pair_exists(pairs, {node1,node2})
                couplings_links(iLayer).to(iLayer+1).pairs(end+1,:) = {node1,node2};
                nb_new = nb_new + 1;
                p_temp = (1.0 * nb_in) / max(1, (nb_in + nb_out + nb_new));
            end
        end
        nb_out = nb_out + nb_new;
        fprintf('Nb of coupling links created (outside cross-layer communities) : %d\n', nb_new);
    end

    fprintf('p value reached : %f\n', p_temp);
end

function out = unique_nodes_in_pairs(pairs, nodes_set, col)
% nodes from nodes_set that actually appear in pairs(:,col)
    used = {};
    for r = 1:size(pairs,1)
        used{end+1} = pairs{r,col}; %#ok<AGROW>
    end
    out = intersect(nodes_set, used);
end

function nodes = get_nodes_from_com(com_nodes_map, com_id)
    if ~isKey(com_nodes_map, com_id)
        nodes = {};
    else
        nodes = com_nodes_map(com_id);
    end
end

function s = sample_cells(list, k)
% sample k elements from cell array list without replacement
    if isempty(list) || k <= 0
        s = {};
        return;
    end
    k = min(k, numel(list));
    idx = randperm(numel(list), k);
    s = list(idx);
end

function write_soum_format(dir_lfr, outdir, nb_layers, layers, ~, couplings_links, nb_com, cross_communities, intra_coms)
% Writes the "new_format" file as in Python
    sep = sprintf('\t');
    f = fopen(fullfile(outdir, 'new_format'), 'w');

    % #layers
    fprintf(f, '%d\n', nb_layers);

    % For each layer: nodes, edges
    for i = 1:nb_layers
        % nodes line
        fprintf(f, '%s\n', strjoin(layers(i).nodes, ' '));

        % edges of layer i
        filelayer = fullfile(dir_lfr, sprintf('layer%d', i-1), 'network.dat');
        edges = read_edges(filelayer, sep);
        fprintf(f, '%d\n', size(edges,1));
        for r = 1:size(edges,1)
            fprintf(f, '%s %s\n', edges{r,1}, edges{r,2});
        end
    end

    % nb couplings (number of cross-layer pairs (l1,l2) with any edges)
    nb_cpl_blocks = 0;
    for l1 = 1:nb_layers-1
        if ~isempty(couplings_links(l1).to(l1+1).pairs)
            nb_cpl_blocks = nb_cpl_blocks + 1;
        end
    end
    fprintf(f, '%d\n', nb_cpl_blocks);

    % coupling edges blocks
    for l1 = 1:nb_layers-1
        pairs = couplings_links(l1).to(l1+1).pairs;
        if isempty(pairs), continue; end
        fprintf(f, '%d %d\n', l1, l1+1); % Python writes l1+1 l2+1 (1-based). Here: (l1) (l2) already 1-based
        fprintf(f, '%d\n', size(pairs,1));
        for r = 1:size(pairs,1)
            fprintf(f, '%s %s\n', pairs{r,1}, pairs{r,2});
        end
    end

    % nb of communities
    fprintf(f, '%d\n', nb_com);

    % First, cross-layer communities (each line concatenates layer-wise node lists)
    for r = 1:size(cross_communities,1)
        c1 = int32(cross_communities(r,1));
        c2 = int32(cross_communities(r,2));
        wrote = false;
        % layer 1 part
        if isKey(layers(1).com_nodes, c1)
            fprintf(f, '%s ', strjoin(layers(1).com_nodes(c1), ' '));
            wrote = true;
        end
        % layer 2 part
        if isKey(layers(2).com_nodes, c2)
            if ~wrote, fprintf(f, ' '); end
            fprintf(f, '%s ', strjoin(layers(2).com_nodes(c2), ' '));
        end
        fprintf(f, '\n');
    end

    % Then, intra-layer communities
    for i = 1:nb_layers
        ids = intra_coms.(sprintf('L%d', i));
        for id = ids
            cid = int32(id);
            if isKey(layers(i).com_nodes, cid)
                fprintf(f, '%s\n', strjoin(layers(i).com_nodes(cid), ' '));
            else
                fprintf(f, '\n');
            end
        end
    end

    fclose(f);
end

function edges = read_edges(fname, sep)
    fid = fopen(fname, 'r');
    if fid < 0
        error('Cannot open %s', fname);
    end
    C = textscan(fid, '%s%s', 'Delimiter', sep, 'CollectOutput', true);
    fclose(fid);
    edges = C{1};
end

