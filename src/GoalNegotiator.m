classdef GoalNegotiator < matlab.mixin.Copyable
    properties (Access = public)
        agents
        num_players {mustBeInteger}
        init_costs {mustBeNumeric}
        % for all deals 
        all_deals_idxs {mustBeInteger}
        all_costs {mustBeNumeric}
        all_self_utils {mustBeNumeric}
        allocation_array {mustBeNumericOrLogical}
        allocation_cost {mustBeNumeric}
        % upto trust level
        all_utils {mustBeNumeric}
        % optimal bools
        max_sw_bools {mustBeNumericOrLogical}
        max_np_bools {mustBeNumericOrLogical}
        pareto_bools {mustBeNumericOrLogical}

        % stage-wise properties
        stg_trust_hist {mustBeNumeric}
        stg_cncsn_hist
        stg_offer_hist
        stg_deal_bool_hist
        stg_NS_bool_hist
        stg_bar_NS_bool_hist
        stg_rNS_bool_hist
        stg_bar_rNS_bool_hist
        stg_cncsn_limit_hist {mustBeNumeric}
        stg_table_bool_hist
        
        % negotiation
        neg_beta_params {mustBeNumeric}
        neg_cand_deal_idx {mustBeInteger}
    end

    properties (Access = private)
        secret_thres_ {mustBeNumeric}
        trans_thres_ {mustBeNumeric}
        learning_rate_ {mustBeNumeric}
    end
    
    methods
        % constructor
        function self = GoalNegotiator(agents)
            % public
            self.agents = agents;
            self.num_players = numel(agents);

            self.neg_beta_params = repmat([1, 1], self.num_players, self.num_players);
            self.neg_cand_deal_idx = 0;

            % private
            self.secret_thres_ = 1.0;
            self.trans_thres_ = 0.0;
            self.learning_rate_ = 0.9;
        end

        % setter
        function set.agents(self, agents)
            self.agents = agents;
            self.UpdateNumPlayer();
            self.ComputeInitialCosts();
        end
    end

    methods (Access = public)
        function self = RunSequentialNegotiation(self, agents)
            if (nargin > 1)
                self.agents = agents;
            end
            if (self.num_players < 2 || self.num_players > 3)
                error("Invalid # of negotiation players!");
            end

            % initialize
            self.stg_trust_hist = zeros(self.num_players, self.num_players, 1);
            self.stg_deal_bool_hist = cell(0, 1);
            self.stg_NS_bool_hist = cell(0, 1);
            self.stg_bar_NS_bool_hist = cell(0, 1);
            self.stg_rNS_bool_hist = cell(0, 1);
            self.stg_bar_rNS_bool_hist = cell(0, 1);
            self.stg_cncsn_limit_hist = zeros(0, self.num_players);
            self.stg_table_bool_hist = cell(0, self.num_players);
            self.stg_cncsn_hist = cell(0, 1);
            self.stg_offer_hist = cell(0, 1);
            num_goals = zeros(1, self.num_players);
            for id = 1:self.num_players
                self.stg_trust_hist(:,id,1) = self.agents(id).trusts;
                num_goals(id) = size(self.agents(id).goals, 1);
            end
            num_total_goals = sum(num_goals);
            
            % sequential negotiation
            stg_num = 1;
            is_conflict = false;
            while (self.num_players*num_total_goals >= stg_num)
                % select table tasks
                goal_neg_bools = self.SetNegotiationTable(stg_num);
                table_bools = zeros(num_total_goals, 1, "logical");
                num_table_goals = 0;
                table_goals_ids = zeros(0, 1);
                table_goals_idxs = zeros(0, 1);
                pntr = 1;
                for id = 1:self.num_players
                    table_bool = (sum(goal_neg_bools{id}, 2) > 1);
                    table_bools(pntr:pntr+size(self.agents(id).goals, 1)-1) = table_bool;
                    table_goals_ids = cat(1, table_goals_ids, id*ones(sum(table_bool), 1));
                    table_goals_idxs = cat(1, table_goals_idxs, find(table_bool));
                    num_table_goals = num_table_goals + sum(table_bool);
                    pntr = pntr + num_goals(id);
                end
                self.stg_table_bool_hist = cat(1, self.stg_table_bool_hist, goal_neg_bools);
                % % if conflict
                if (is_conflict)
                    break;
                end
                % if there is no table tasks
                if (~num_table_goals)
                    break;
                end
                % if there is no additional table tasks at the non-initial stage and no trust evolution
                is_stalled = true;
                if (stg_num > 1)
                    for id = 1:self.num_players
                        if (~isequal(goal_neg_bools{end, id}, self.stg_table_bool_hist{end-1, id}))
                            is_stalled = false;
                            break;
                        end
                    end
                    if (is_stalled)
                        break;
                    end
                end

                % compute all possible partitions
                stg_perms = GetRepeatedPermutation(1:self.num_players, num_table_goals);
                stg_num_perms = size(stg_perms, 1);
                % check validity of the permutation
                valid_bools = ones(stg_num_perms, 1, "logical");
                for perm_id = 1:stg_num_perms
                    perm = stg_perms(perm_id, :);
                    for g_id = 1:numel(perm)
                        assigned_id = perm(g_id);
                        goal_id = table_goals_ids(g_id);
                        goal_idx = table_goals_idxs(g_id);
                        if ~(ismember(assigned_id, find(goal_neg_bools{goal_id}(goal_idx, :))))
                            valid_bools(perm_id) = false;
                            break;
                        end
                    end
                    if (~valid_bools(perm_id))
                        continue;
                    end
                end
                stg_perms = stg_perms(valid_bools, :);
                stg_num_perms = size(stg_perms, 1);

                stg_deals_idx = zeros(stg_num_perms, self.num_players);
                % compute cost functions

                for perm_id = 1:stg_num_perms
                    perm = stg_perms(perm_id, :);
                    assign_bools = zeros(self.num_players, num_total_goals, "logical");
                    pntr = 1;
                    for id = 1:self.num_players
                        assign_bools(id, pntr:pntr+num_goals(id)-1) = true;
                        assign_bools(id, table_bools) = (perm == id);
                        pntr = pntr + num_goals(id);
                    end
                    stg_deals_idx(perm_id, :) = AlocVec2Index(assign_bools).';
                end

                stg_deals_bools = zeros(size(self.all_deals_idxs, 1), 1, "logical");
                [~,deal_idxs] = ismember(stg_deals_idx, self.all_deals_idxs, "rows");
                stg_deals_bools(deal_idxs) = true;
                self.stg_deal_bool_hist = cat(1, self.stg_deal_bool_hist, stg_deals_bools);

                % Run constrained MCP
                is_conflict = self.RunMCP(stg_deals_bools);
                
                % Draw debug plot
                % self.DrawDebugPlot(100);

                % Trust evloution
                update_trust = self.EvolveTrust();
                self.stg_trust_hist = cat(3, self.stg_trust_hist, update_trust);
                
                stg_num = stg_num + 1;
            end
    
            if (self.neg_cand_deal_idx(end) > 0)
                self.ExecuteOffer();
            end
        end

        function self = ComputeGlobalOptima(self, agents)
            if (nargin > 1)
                self.agents = agents;
            end
            self.GenerateAllDeals();

            % all possible partitions of total goal set
            social_welfare = sum(self.all_self_utils, 2);
            pos_utils = self.all_self_utils.*(self.all_self_utils >= 0);
            nash_product = prod(pos_utils, 2);
            self.max_sw_bools = (social_welfare == max(social_welfare));
            self.max_np_bools = (nash_product == max(nash_product));
        end

        function pareto_bools = ComputeParetoDeals(self, all_utils)
            utils = self.all_self_utils;
            if nargin == 2
                utils = all_utils;
            end
            num_deals = size(utils, 1);
            pareto_bools = 1:num_deals;
            nxt_pt_idx = 1;
            while nxt_pt_idx <= size(utils, 1)
                non_dom_pt_mask = any(utils > utils(nxt_pt_idx, :), 2);
                non_dom_pt_mask(nxt_pt_idx) = true;
                pareto_bools = pareto_bools(non_dom_pt_mask);
                utils = utils(non_dom_pt_mask, :);
                nxt_pt_idx = sum(non_dom_pt_mask(1:nxt_pt_idx)) + 1;
            end
            is_efficient_mask = false(num_deals, 1);
            is_efficient_mask(pareto_bools) = true;
            pareto_bools = is_efficient_mask;
        end
        
        function DrawDebugPlot(self, fid)
            num_stg = size(self.stg_trust_hist, 3);
            max_self_util = max(self.all_self_utils);
            min_self_util = min(self.all_self_utils);
            for i = 1:num_stg
                trust_mtx = self.stg_trust_hist(:, :, i);
                feas_bools = self.stg_deal_bool_hist{i};
                feas_self_utils = self.all_self_utils(feas_bools, :);
                cncsn_limits = self.stg_cncsn_limit_hist(i, :);
                NS_bools = self.stg_NS_bool_hist{i};
                rNS_bools = self.stg_rNS_bool_hist{i};
                self.ComputeAllUtils(trust_mtx);

                % trust guide
                p1_coord = - trust_mtx(2, 1) * max_self_util(2);
                p2_coord = - trust_mtx(1, 2) * max_self_util(1);

                figure(fid)
                fid = fid + 1;
                subplot(1,2,1)
                % draw all possible deals
                plot(self.all_self_utils(:,1), self.all_self_utils(:,2), ".");
                hold on
                plot(feas_self_utils(:, 1), feas_self_utils(:, 2), "g.");
                plot(feas_self_utils(NS_bools, 1), feas_self_utils(NS_bools, 2), "r.");
                plot(feas_self_utils(rNS_bools, 1), feas_self_utils(rNS_bools, 2), "ko");
                plot(self.all_self_utils(self.max_sw_bools, 1), self.all_self_utils(self.max_sw_bools, 2), "r^");
                plot(0, 0, "kx");
                ref_utils = zeros(1, self.num_players);
                for j = 1:i
                    if (numel(self.neg_cand_deal_idx) < 1+j)
                        continue;
                    end
                    nxt_utils = self.all_self_utils(self.neg_cand_deal_idx(1+j), :);
                    util_hist = [ref_utils; nxt_utils];
                    ref_utils = nxt_utils;
                    plot(util_hist(:, 1), util_hist(:, 2), ":m", "linewidth", 1.5);
                end
                plot([0, p1_coord], [0, max_self_util(2)], "--k");
                plot([0, max_self_util(1)], [0, p2_coord], "--k");
                hold off
                axis equal
                grid on
                legend("All deals",...
                    "Feasible deals",...
                    "Negotiation Set (NS)",...
                    "Reduced NS",...
                    "Max Social Welfare",...
                    "Conflict",...
                    "Negotiation",...
                    "Mutual benefit",...
                    "Location",...
                    "best")
                xlabel('$$\tilde{u}^{(1)}$$','fontsize',14,'Interpreter', 'LaTeX')
                ylabel('$$\check{\tilde{u}}^{(2)}$$','fontsize',14,'Interpreter', 'LaTeX')
                xlim([min_self_util(1), max_self_util(1)])
                ylim([min_self_util(2), max_self_util(2)])
                
                subplot(1,2,2)
                feas_utils = self.all_utils(feas_bools, :);
                plot_list = [];
                legend_list = [];
                p = plot(self.all_utils(:, 1), self.all_utils(:, 2), ".");
                plot_list = cat(2, plot_list, p);
                legend_list = cat(2, legend_list, "All deals");
                hold on
                p = plot(feas_utils(:, 1), feas_utils(:, 2), "g.");
                plot_list = cat(2, plot_list, p);
                legend_list = cat(2, legend_list, "Feasible deals");
                p = plot(feas_utils(NS_bools, 1), feas_utils(NS_bools, 2), "r.");
                plot_list = cat(2, plot_list, p);
                legend_list = cat(2, legend_list, "NS deals");
                p = plot(feas_utils(rNS_bools, 1), feas_utils(rNS_bools, 2), "ko");
                plot_list = cat(2, plot_list, p);
                legend_list = cat(2, legend_list, "Reduced NS deals");
                [~, max_idx] = max(prod(feas_utils(NS_bools, :), 2));
                NS_idxs = find(NS_bools);
                max_idx = NS_idxs(max_idx);
                p = plot(feas_utils(max_idx, 1), feas_utils(max_idx, 2), "r^");
                plot_list = cat(2, plot_list, p);
                legend_list = cat(2, legend_list, "Max Nash product");
                p = plot(0, 0, "kx");
                plot_list = cat(2, plot_list, p);
                legend_list = cat(2, legend_list, "Conflict");

                ref_utils = zeros(1, self.num_players);
                for j = 1:i
                    if (numel(self.neg_cand_deal_idx) < 1+j)
                        continue;
                    end
                    nxt_utils = self.all_utils(self.neg_cand_deal_idx(1+j), :);
                    util_hist = [ref_utils; nxt_utils];
                    ref_utils = nxt_utils;
                    p = plot(util_hist(:, 1), util_hist(:, 2), ":m", "linewidth", 1.5);
                end
                plot_list = cat(2, plot_list, p);
                legend_list = cat(2, legend_list, "Negotiation");
                
                p1_min_util = min(self.all_utils(:, 1));
                p2_min_util = min(self.all_utils(:, 2));
                p1_max_util = max(self.all_utils(:, 1));
                p2_max_util = max(self.all_utils(:, 2));
                p = plot([p1_min_util, p1_max_util], cncsn_limits(2)*[1, 1], "-.b");
                plot_list = cat(2, plot_list, p);
                legend_list = cat(2, legend_list, "Concession limit");
                plot(cncsn_limits(1)*[1, 1], [p2_min_util, p2_max_util], "-.b");
                p = plot([0, p1_max_util], 0*[1, 1], "--k");
                plot_list = cat(2, plot_list, p);
                legend_list = cat(2, legend_list, "Mutual benefit");
                plot(0*[1, 1], [0, p2_max_util], "--k");
                hold off
                axis equal
                grid on
                legend(plot_list, legend_list, "Location", "best")
                title(strcat("$\rho^{(1)} = ", num2str(trust_mtx(2,1)),...
                        ", \quad \rho^{(2)} = ", num2str(trust_mtx(1,2)), "$"),...
                        "Interpreter", "latex")
                xlabel('${u}^{(1)}$','fontsize',14,'interpreter','latex')
                % ylabel('$\check{\tilde{u}}^{(2)}$','fontsize',14,'interpreter','latex')
                ylabel('${u}^{(2)}$','fontsize',14,'interpreter','latex')
                xlim([p1_min_util, p1_max_util])
                ylim([p2_min_util, p2_max_util])
            end
        end

        function self = ComputeAllUtils(self, trust_mtx)
            if nargin < 2
                trust_mtx = zeros(self.num_players, self.num_players);
                for id = 1:self.num_players
                    trust_mtx(:, id) = self.agents(id).trusts;
                end
            end
            self.all_utils = self.all_self_utils * trust_mtx;
        end

        function goals = Array2Goals(self, array)
            goals = zeros(0, 2);
            pntr = 0;
            for id = 1:self.num_players
                num_goals = size(self.agents(id).goals, 1);
                add_task_idx = find(array(pntr + (1:num_goals)));
                if ~isempty(add_task_idx)
                    goals = cat(1, goals, self.agents(id).goals(add_task_idx, :));
                end
                pntr = pntr + num_goals;
            end
        end
    end

    methods (Access = private)
        function self = ComputeInitialCosts(self)
            self.init_costs = zeros(self.num_players, 1);
            for i = 1:self.num_players
                self.init_costs(i) = self.SolveTSP_DP(self.agents(i).dist_table);
            end
        end

        function self = GenerateAllDeals(self)
            total_goals = zeros(0, 2);
            total_goals_idxs = zeros(0, 2, "uint8");
            for id = 1:self.num_players
                total_goals = cat(1, total_goals, self.agents(id).goals);
                num_goals = size(self.agents(id).goals, 1);
                goals_idxs = [repmat(id, num_goals, 1), (1:num_goals).'];
                total_goals_idxs = cat(1, total_goals_idxs, goals_idxs);
            end

            num_total_goals = size(total_goals, 1);
            self.allocation_array = Index2AlocVec(1:2^(num_total_goals), num_total_goals);
            self.allocation_cost = zeros(size(self.allocation_array, 1), self.num_players);
            for array_id = 1:size(self.allocation_array, 1)
                goals = self.Array2Goals(self.allocation_array(array_id, :));
                for id = 1:self.num_players
                    copy_agent = copy(self.agents(id));
                    copy_agent.goals = goals;
                    self.allocation_cost(array_id, id) = self.SolveTSP_DP(copy_agent.dist_table);
                end
            end

            % all possible partitions of total goal set
            goals_perms = GetRepeatedPermutation(1:self.num_players, size(total_goals, 1));
            num_deals = size(goals_perms, 1);
            self.all_deals_idxs = zeros(num_deals, self.num_players);
            self.all_costs = zeros(num_deals, self.num_players);
            self.all_self_utils = zeros(num_deals, self.num_players);
            all_deals_idxs = zeros(num_deals, self.num_players);
            num_players = self.num_players;
            % Parallelize all deal computation
            parfor deal_id = 1:num_deals
                perm = goals_perms(deal_id, :);
                perm_aloc_vec = zeros(0, num_total_goals);
                for id = 1:num_players
                    perm_aloc_vec = cat(1, perm_aloc_vec, perm == id);
                end
                allocation_idxs = AlocVec2Index(perm_aloc_vec);
                all_deals_idxs(deal_id, :) = allocation_idxs.';
            end
            self.all_deals_idxs = all_deals_idxs;
            for id = 1:self.num_players
                self.all_costs(:,id) = self.allocation_cost(self.all_deals_idxs(:,id), id);
                self.all_self_utils(:, id) = self.init_costs(id) - self.all_costs(:,id);
            end
            % Compute cost and check pareto optimality
            self.ComputeAllUtils();
            self.pareto_bools = self.ComputeParetoDeals();
        end

        function cost = SolveTSP_DP(self, dist)
            % initialize Process
            cost = 0;
            num_dests = size(dist, 1);
            if (num_dests <= 1)
                return;
            end
            Primes = primes(num_dests * 10);

            % move ego position to the first entry
            dist = [dist(:, end), dist(:, 1:end-1)];
            dist = [dist(end, :); dist(1:end-1, :)];
            % remove diagonal entry
            dist(boolean(eye(num_dests))) = Inf * ones(1, num_dests);
            dist_table = dist;
            % nullify distance to origin
            for i = 2:num_dests
                dist_table(i, 1) = 0;
            end

            num_data_sets = 1;
            for i = 2:num_dests
                num_data_sets = num_data_sets + nchoosek(num_dests, i);
            end
            node_data(num_data_sets).S = [];
            node_data(num_data_sets).l = 0;
            node_data(num_data_sets).cost = Inf;
            node_data(num_data_sets).pre = [];
            node_data(num_data_sets).m = [];
            look_up_table(num_data_sets) = 0;
            node_data(1).S = 1;
            node_data(1).l = 1;
            node_data(1).cost = 0;
            node_data(1).Pre = [];
            node_data(1).m = [];
            for s = 2:num_dests
                node_data(s).S = cat(2, node_data(1).S, s);
                node_data(s).l = s;
                node_data(s).cost = dist_table(1, s);
                node_data(s).Pre = 1;
                node_data(s).m = 1;
                LUT = Compute_LUT(node_data(s).S, s, Primes);
                look_up_table(s) = LUT;
            end
            % index into node_data that marks the beginning of the previous step
            idx_st_prv_step = 2;
            % index into node_data that marks the end of the previous step
            idx_last_step = num_dests;
            % index into node_data that marks the current dataset
            cur_data = idx_last_step;
            
            % Dynamic Programming loop
            candidate_cost = [];
            for s = 3:num_dests
                % generate possible sets with s-1 dest out of the possible N-1 dest
                temp_sets = nchoosek(2:num_dests, s-1);
                num_sets = size(temp_sets);
                % iterate through all the sets
                for j = 1:num_sets(1)
                    % iterate through all the elements in each set
                    for k = 1:num_sets(2)
                        % this is the set S-{k}
                        minus_set = [1, temp_sets(j,1:k-1), temp_sets(j, k+1:num_sets(2))];
                        candidate_cost = zeros(1, length(minus_set));
                        candidate_cost(2:length(minus_set)) = Inf;
                        indices = zeros(1, length(minus_set));
                        % choose a dest in S-{k} that will be last
                        for mm = 2:length(minus_set)
                            look_up_val = Compute_LUT(minus_set, minus_set(mm), Primes);
                            index = find(look_up_val == look_up_table(idx_st_prv_step:idx_last_step));
                            index = index + idx_st_prv_step - 1;
                            if (index == 0)
                                candidate_cost(mm) = Inf;
                            else
                                candidate_cost(mm) = node_data(index).cost + dist_table(temp_sets(j,k), minus_set(mm));
                                indices(mm) = index;
                            end
                        end
                        [mincost, indexcost] = min(candidate_cost(2:end));
                        cur_data = cur_data + 1;
                        node_data(cur_data).S = [1, temp_sets(j, :)];
                        node_data(cur_data).l = temp_sets(j, k);
                        node_data(cur_data).cost = mincost;
                        node_data(cur_data).Pre = indices(indexcost + 1);
                        node_data(cur_data).m = minus_set(indexcost + 1);
                        look_up_table(cur_data) = Compute_LUT(node_data(cur_data).S, temp_sets(j,k), Primes);
                    end
                end
                idx_st_prv_step = idx_last_step + 1;
                idx_last_step = cur_data;
            end
            mm = 0;
            % add the distance back from the last dest to dest 1
            for i = idx_st_prv_step:idx_last_step
                mm = mm + 1;
                candidate_cost(mm) = node_data(i).cost + dist_table(node_data(i).l, 1);
            end
            % find the one that minimizes the total distance
            [cost, ~] = min(candidate_cost);
        end

        function goal_neg_bools = SetNegotiationTable(self, stg_num)
            % select table tasks
            goal_neg_bools = cell(1, self.num_players);
            for id = 1:self.num_players
                trust_bools = self.agents(id).GetTrustTask();
                if (stg_num == 1)
                    % initial negotiation tasks
                    % add all confidential tasks lower than trust
                    goal_neg_bools{id} = trust_bools;
                else
                    trans_bool = self.agents(id).GetTransTask(self.trans_thres_);
                    trans_bools = repmat(trans_bool, [1, self.num_players]);
                    secret_bool = self.agents(id).GetSecretTask(self.secret_thres_);
                    secret_bools = repmat(secret_bool, [1, self.num_players]);
                    non_trans_bool = (~trans_bools) & (~secret_bools) & (trust_bools);
                    % add all transparent goals
                    goal_neg_bools{id} = self.stg_table_bool_hist{end, id};
                    % add one non_transparent goal
                    addable_bool = ~goal_neg_bools{id} & non_trans_bool;
                    task_add_bool = (sum(addable_bool, 2) > 0);
                    if (any(task_add_bool))
                        task_add_idxs = find(task_add_bool);
                        [~,odr] = min(self.agents(id).goal_confs(task_add_bool, :));
                        task_add_idx = task_add_idxs(odr);
                        opponent_idxs = find(addable_bool(task_add_idx, :));
                        [~,odr] = max(self.agents(id).trusts(opponent_idxs));
                        opponent_idx = opponent_idxs(odr);
                        goal_neg_bools{id}(task_add_idx, opponent_idx) = true;
                    end
                end
            end
        end

        function update_trust = EvolveTrust(self)
            % initialize cross trust
            update_trust = self.stg_trust_hist(:, :, end);
            % If the negotiation has single NS deal
            if (sum(self.stg_NS_bool_hist{end}) <= 1)
                return;
            end
            
            feas_utils = self.all_utils(self.stg_deal_bool_hist{end}, :);
            NS_utils = feas_utils(self.stg_NS_bool_hist{end}, :);
            [~, idx] = max(prod(NS_utils, 2));
            agrd_utils = NS_utils(idx, :);
            cncsn_ratio = ones(1, self.num_players);
            for id = 1:self.num_players
                max_util = max(NS_utils(:,id));
                if (max_util < eps)
                    continue;
                end
                cncsn_ratio(id) = (max_util - agrd_utils(id)) / max_util;
            end
            
            num_total_concsn = sum(self.stg_NS_bool_hist{end}) - 1;
            beta_params = self.neg_beta_params(:, :, end);
            update_beta_params = beta_params;
            for id = 1:self.num_players
                if (self.agents(id).deception_bool)
                    continue;
                end
                for op_id = 1:self.num_players
                    % skip self-confidence
                    if (id == op_id)
                        continue;
                    end
                    if (cncsn_ratio(id) == 0 && cncsn_ratio(op_id) == 0)
                        continue;
                    else
                        fairness = cncsn_ratio(id) / (cncsn_ratio(id) + cncsn_ratio(op_id));
                        update_beta_params(id, 2*(op_id-1)+[1, 2]) = ...
                                        beta_params(id, 2*(op_id-1)+[1, 2]) +...
                                        num_total_concsn * [1-fairness, fairness];
                        posterior = update_beta_params(id, 2*(op_id-1)+1) / sum(update_beta_params(id, 2*(op_id-1)+[1, 2])); % Expectation
                        trust_ref = min(1, 2*posterior);
                        update_trust(op_id, id) = (trust_ref - update_trust(op_id, id)) * self.learning_rate_ + update_trust(op_id, id);
                    end
                    self.agents(id).trusts(op_id) = update_trust(op_id, id);
                end
            end
            self.neg_beta_params = cat(3, self.neg_beta_params, update_beta_params);
        end

        function self = ExecuteOffer(self, offer_idx)
            if nargin < 2
                offer_idx = self.neg_cand_deal_idx(end);
            end
            num_goals = zeros(1, self.num_players);
            for id = 1:self.num_players
                num_goals(id) = size(self.agents(id).goals, 1);
            end
            agrd_aloc_vecs = Index2AlocVec(self.all_deals_idxs(offer_idx, :), sum(num_goals));
            agrd_deal = cell(1, self.num_players);
            for id = 1:self.num_players
                agrd_deal{id} = self.Array2Goals(agrd_aloc_vecs(id, :));
            end
            for id = 1:self.num_players
                self.agents(id).goals = agrd_deal{id};
            end
        end

        function is_conflict = RunMCP(self, deal_bools)
            is_conflict = false;
            % compute MCP utils for all possible deals
            deal_idxs = find(deal_bools);
            self.ComputeAllUtils();
            stg_self_utils = self.all_self_utils(deal_bools, :);
            stg_utils = self.all_utils(deal_bools, :);

            % compute concession guide
            MB_bools = prod(stg_utils >= 0, 2);
            PO_bools = self.ComputeParetoDeals(stg_utils);
            NS_bools = PO_bools & MB_bools;
            rNS_bools = NS_bools;
            bar_NS_bools = NS_bools;
            bar_rNS_bools = rNS_bools;
            if self.num_players > 2
                bar_NS_bools = zeros(size(NS_bools), "logical");
                cg_bools = self.ComputeConcessionGuide(deal_bools);
                % Compute negotiation set
                cg_idxs = find(cg_bools);
                cg_utils = stg_utils(cg_bools, :);
                cg_PO_bools = self.ComputeParetoDeals(cg_utils);
                cg_MB_bools = prod(cg_utils >= 0, 2);
                prt_cg_idx = find(cg_PO_bools & cg_MB_bools);
                bar_NS_bools(cg_idxs(prt_cg_idx)) = true;
                bar_rNS_bools = bar_NS_bools;
            end
            cncsn_limits = zeros(1, self.num_players);
            cncsn_hist = zeros(1, self.num_players);
            for id = 1:self.num_players
                cncsn_limits(id) = self.agents(id).ConstructConcessionLimit(stg_self_utils, stg_utils, NS_bools);
                val_bools = (cncsn_limits(id) < stg_utils(:, id));
                rNS_bools = rNS_bools & val_bools;
                bar_rNS_bools = bar_rNS_bools & val_bools;
                cncsn_hist(id) = sum(~(cncsn_limits(id) < stg_utils(NS_bools, id))) + 1;
            end

            self.stg_NS_bool_hist = cat(1, self.stg_NS_bool_hist, NS_bools);
            self.stg_bar_NS_bool_hist = cat(1, self.stg_bar_NS_bool_hist, bar_NS_bools);
            self.stg_rNS_bool_hist = cat(1, self.stg_rNS_bool_hist, rNS_bools);
            self.stg_bar_rNS_bool_hist = cat(1, self.stg_bar_rNS_bool_hist, bar_rNS_bools);
            self.stg_cncsn_limit_hist = cat(1, self.stg_cncsn_limit_hist, cncsn_limits);

            rNS_deal_idxs = deal_idxs(bar_rNS_bools);
            rNS_utils = stg_utils(bar_rNS_bools, :);

            offer_hist = zeros(0, self.num_players);
            rNS_utils_copy = rNS_utils;
            nxt_offer_idxs = zeros(1, self.num_players);
            if ~any(bar_rNS_bools)
                is_conflict = true;
                return;
            end
            while any(rNS_utils_copy ~= -Inf, "all")
                % if no agent concedes, conflict
                if (~any(cncsn_hist(end, :)))
                    self.stg_cncsn_hist = cat(1, self.stg_cncsn_hist, cncsn_hist);
                    self.stg_offer_hist = cat(1, self.stg_offer_hist, offer_hist);
                    break;
                end
                % Offer generation
                offer_idxs = self.GenerateMCPOffers(rNS_utils_copy, nxt_offer_idxs);
                offer = zeros(1, self.num_players);
                for id = 1:self.num_players
                    rNS_utils_copy(offer_idxs(id), id) = -Inf;
                    offer(id) = rNS_deal_idxs(offer_idxs(id));
                end
                offer_hist = cat(1, offer_hist, offer);

                % Offer evaluation
                ex_deal_idx = self.EvaluateMCPOffers(offer_idxs, rNS_utils);

                % if a deal is agreed
                if (ex_deal_idx > 0)
                    self.neg_cand_deal_idx = cat(2, self.neg_cand_deal_idx, rNS_deal_idxs(ex_deal_idx));
                    self.stg_cncsn_hist = cat(1, self.stg_cncsn_hist, cncsn_hist);
                    self.stg_offer_hist = cat(1, self.stg_offer_hist, offer_hist);
                    break;
                end

                % Compute Zeuthen risk criteria
                cncsn_id = self.GetZeuthenCncsn(offer_idxs, rNS_utils);
                for id = 1:self.num_players
                    nxt_offer_idxs(id) = offer_idxs(id);
                    if (id == cncsn_id)
                        [~, nxt_offer_idxs(id)] = max(rNS_utils_copy(:, id));
                    end
                end
                concsn_array = zeros(1, self.num_players, "logical");
                concsn_array(cncsn_id) = true;
                cncsn_hist = cat(1, cncsn_hist, concsn_array);
            end
        end

        function cg_bools = ComputeConcessionGuide(self, deal_bools)
            stg_self_utils = self.all_self_utils(deal_bools, :);
            % compute concession guide
            stg_index = self.all_deals_idxs(deal_bools, :);
            cg_bools = ones(size(stg_self_utils, 1), 1, "logical");
            for id = 1:self.num_players
                if (self.agents(id).deception_bool)
                    continue;
                end
                deal_combs = stg_index;
                checked_bool = zeros(size(stg_index, 1), 1, "logical");
                within_bool = zeros(size(stg_index, 1), 1, "logical");
                op_ids = setdiff(1:self.num_players, id);
                for ii = 1:size(deal_combs, 1)
                    if checked_bool(ii)
                        continue;
                    end
                    sel_bools = (deal_combs(:,id) == deal_combs(ii, id));
                    fam_within_bool = zeros(size(within_bool), "logical");
                    fam_within_bool(sel_bools) = true;
                    checked_bool(sel_bools) = true;
                    sel_idxs = find(sel_bools);
                    sel_utils = stg_self_utils(sel_bools, :);
                    nor_sel_utils = sel_utils - min(sel_utils);
                    prt_PO_bools = self.ComputeParetoDeals(nor_sel_utils(:, op_ids));
                    prt_MB_bools = prod(nor_sel_utils(:, op_ids) >= 0, 2);
                    prt_sel_idx = find(prt_PO_bools & prt_MB_bools);
                    prt_sel_bools = zeros(size(sel_bools), "logical");
                    prt_sel_bools(sel_idxs(prt_sel_idx)) = true;
                    prt_sel_utils = stg_self_utils(prt_sel_bools, :);
                    for op_id = 1:self.num_players
                        if (id == op_id)
                            continue;
                        end
                        max_util = max(stg_self_utils(:, op_id));
                        min_util = min(prt_sel_utils(:, op_id));
                        if (numel(sel_idxs) <= 1 && max_util == max(stg_self_utils(:,op_id)))
                            fam_within_bool(sel_bools) = false;
                            break;
                        end
                        ref_util = (max_util - min_util)/(self.num_players-1) + min_util;
                        lim_util = (max_util - ref_util) * (self.agents(id).trusts(op_id))^(sum(deal_bools)/size(self.all_deals_idxs,1)) + ref_util;
                        lim_bools = (sel_utils(:,op_id) <= lim_util);
                        lim_idxs = (sel_idxs(lim_bools));
                        temp_bool = zeros(size(fam_within_bool), "logical");
                        temp_bool(lim_idxs) = true;
                        fam_within_bool = fam_within_bool & (temp_bool & prt_sel_bools);
                    end
                    within_bool = within_bool | fam_within_bool;
                end
                cg_bools = cg_bools & within_bool;
            end
        end

        function offer_idxs = GenerateMCPOffers(self, utils, nxt_offer_idxs)
            offer_idxs = nxt_offer_idxs;
            for id = 1:self.num_players
                % stick to previous offer
                if (nxt_offer_idxs(id) > 0)
                    continue;
                end
                % no more beneficial deal
                if ~sum(utils(:, id) >= 0)
                    continue;
                end
                [~, offer_idxs(id)] = max(utils(:, id));
            end
        end

        function exct_idx = EvaluateMCPOffers(self, offer_idxs, utils)
            exct_idx = 0;
            accept_bools = zeros(self.num_players, self.num_players, "logical");
            for id = 1:self.num_players
                for op_id = 1:self.num_players
                    if (id == op_id)
                        continue;
                    end
                    if (utils(offer_idxs(id), id) <= utils(offer_idxs(op_id), id))
                        accept_bools(id, op_id) = true;
                    end
                end
            end
            accepted_offeror = (sum(accept_bools) == self.num_players-1);
            if any(accepted_offeror)
                % both accept deals
                offerors_id = find(accepted_offeror);
                if (numel(offerors_id) > 1)
                    exct_idx = offer_idxs(randsample(offerors_id, 1));
                else
                    exct_idx = offer_idxs(offerors_id);
                end
            end
        end

        function self = UpdateNumPlayer(self)
            self.num_players = numel(self.agents);
        end

        function cncsn_id = GetZeuthenCncsn(self, offer_idxs, rNS_utils)
            gen_risks = zeros(self.num_players);
            for id = 1:self.num_players
                nxt_id = rem(id, self.num_players) + 1;
                ref_deal_util = rNS_utils(offer_idxs(id), id);
                com_deal_util = rNS_utils(offer_idxs(nxt_id), id);
                gen_risks(id) = 1;
                if (ref_deal_util > 0)
                    gen_risks(id) = (ref_deal_util - com_deal_util) / ref_deal_util;
                end
            end
            if (gen_risks(1) > gen_risks(2))
                cncsn_id = 2;
            elseif (gen_risks(1) < gen_risks(2))
                cncsn_id = 1;
            else
                cncsn_id = randi([1,2]);
            end
        end
    end
end
