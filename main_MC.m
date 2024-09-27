clear
close all
clc

addpath(genpath(strcat(pwd, "\util")));
addpath(genpath(strcat(pwd, "\src")));

%% Map generation
T = readtable("Init_config.xlsx");
world_size = 10;

% goal map
goal_pos = [T.goalPositionX, T.goalPositionY];
goal_confs = T.goalConfidentiality;
num_goal = size(goal_pos, 1);
init_alloc = T.initialAssignment;

%% Initialize agent position and task
num_players = 2;
player_pos = [T.agentPositionX, T.agentPositionY];
player_pos = player_pos(1:num_players, :);

goal_sets = cell(num_players, 1);
goal_conf_sets = cell(num_players, 1);
for id = 1:num_players
    goal_sets{id} = goal_pos((init_alloc == id), :);
    goal_conf_sets{id} = goal_confs(init_alloc == id);
end

%% Initialize fleet manager
agents = [];
for id = 1:num_players
    agent = Agent(id, num_players, world_size, player_pos(id, :), goal_sets{id});
    agent.goal_confs = goal_conf_sets{id};
    agents = cat(1, agents, agent);
end
fleet_manager = FleetManager(copy(agents), world_size, goal_pos, player_pos);

% Draw initial map
fleet_manager.DrawGlobalMap();
title("Initial Allocation")

%% Global optimum
goal_negotiator = GoalNegotiator(copy(agents));
goal_negotiator.ComputeGlobalOptima();

for i = find(goal_negotiator.max_sw_bools)
    deal_idxs = goal_negotiator.all_deals_idxs(i, :);
    aloc_vecs = Index2AlocVec(deal_idxs, num_goal);
    goal = cell(1, num_players);
    for id = 1:num_players
        goal{id} = goal_negotiator.Array2Goals(aloc_vecs(id, :));
    end
    for id = 1:num_players
        fleet_manager.agents(id).goals = goal{id};
    end
    fleet_manager.DrawGlobalMap();
    title("Max Social Welfare")
end

for i = find(goal_negotiator.max_np_bools)
    deal_idxs = goal_negotiator.all_deals_idxs(i, :);
    aloc_vecs = Index2AlocVec(deal_idxs, num_goal);
    goal = cell(1, num_players);
    for id = 1:num_players
        goal{id} = goal_negotiator.Array2Goals(aloc_vecs(id, :));
    end
    for id = 1:num_players
        fleet_manager.agents(id).goals = goal{id};
    end
    fleet_manager.DrawGlobalMap();
    title("Max Nash Product")
end
ref_utils = goal_negotiator.all_utils;
ref_self_utils = goal_negotiator.all_self_utils;

%% Monte-Carlo Simulation
num_monte = 1e+3;

% sample generator
init_goal_confs = readmatrix("MC_init_goal_confs.xlsx");
P1_init_goal_confs = init_goal_confs(1:size(agents(1).goals, 1), :);
P2_init_goal_confs = init_goal_confs(size(agents(1).goals, 1)+1:end, :);

[~, conflict_offer_idx] = ismember(zeros(1, num_players), goal_negotiator.all_self_utils, "rows");
final_utils = zeros(num_monte, 2);
final_self_utils = zeros(num_monte, 2);
final_trusts = zeros(num_monte, 2);
final_neg_rnds = zeros(num_monte, 1);
final_N_tasks = zeros(num_monte, num_players);
final_deal_idx = zeros(num_monte, 1);

ref_gn = copy(goal_negotiator);
ref_gn.agents = copy(agents);

parfor ii = 1:num_monte
    disp(strcat(num2str(ii), " / ", num2str(num_monte)))
    test_gn = copy(ref_gn);
    test_gn.agents = copy(agents);

    % Choose one of the presets
    % PRESET: P2 is cooperative [Example V.B.1]
    for id = 1:num_players
        test_gn.agents(id).deception_bool = false;
        test_gn.agents(id).trusts = 0.5 * ones(num_players, 1);
        test_gn.agents(id).trusts(id) = 1;
        if id == 1
            test_gn.agents(id).goal_confs = P1_init_goal_confs(:,ii);
        elseif id == 2
            test_gn.agents(id).goal_confs = P2_init_goal_confs(:,ii);
        end
    end

    % PRESET: P2 is selfish, constant inflation on self-util of (all, nothing) [Example V.B.2.+]
    % test_gn.agents(1).deception_bool = false;
    % test_gn.agents(2).deception_bool = true;
    % test_gn.agents(1).trusts(2) = 0.5;
    % test_gn.agents(2).trusts(1) = 0.0;
    % test_gn.agents(1).goal_confs = P1_init_goal_confs(:,ii);
    % test_gn.agents(2).goal_confs = zeros(size(test_gn.agents(2).goal_confs));
    % add_factor = 10;
    % [~,lie_idx] = max(test_gn.all_self_utils(:, 2));
    % test_gn.all_self_utils(lie_idx,2) = test_gn.all_self_utils(lie_idx,2) + add_factor;
    % test_gn.all_costs(lie_idx,2) = test_gn.all_costs(lie_idx,2) - add_factor;

    % PRESET: P2 is selfish, double inflation on positive self-util deals [Example V.B.2.x]
    % test_gn.agents(1).deception_bool = false;
    % test_gn.agents(2).deception_bool = true;
    % test_gn.agents(1).trusts(2) = 0.5;
    % test_gn.agents(2).trusts(1) = 0.0;
    % test_gn.agents(1).goal_confs = P1_init_goal_confs(:,ii);
    % test_gn.agents(2).goal_confs = zeros(size(test_gn.agents(2).goal_confs));
    % mult_factor = 2.0;
    % lie_idxs = test_gn.all_self_utils(:, 2) > 0;
    % test_gn.all_self_utils(lie_idxs, 2) = test_gn.all_self_utils(lie_idxs, 2) * mult_factor;
    % test_gn.all_costs(lie_idxs, 2) = test_gn.init_costs(2) - test_gn.all_self_utils(lie_idxs, 2);

    % PRESET: P2 is selfish, constant inflation, P1 has no anti-deception measures [Example V.C.+]
    % test_gn.agents(1).deception_bool = true;
    % test_gn.agents(2).deception_bool = true;
    % test_gn.agents(2).trusts(1) = 0.0;
    % test_gn.agents(1).goal_confs = P1_init_goal_confs(:,ii);
    % test_gn.agents(2).goal_confs = zeros(size(test_gn.agents(2).goal_confs));
    % add_factor = 10;
    % [~,lie_idx] = max(test_gn.all_self_utils(:, 2));
    % test_gn.all_self_utils(lie_idx,2) = test_gn.all_self_utils(lie_idx,2) + add_factor;
    % test_gn.all_costs(lie_idx,2) = test_gn.all_costs(lie_idx,2) - add_factor;

    % PRESET: P2 is selfish, double inflation, P1 has no anti-deception measures [Example V.C.x]
    % test_gn.agents(1).deception_bool = true;
    % test_gn.agents(2).deception_bool = true;
    % test_gn.agents(2).trusts(1) = 0.0;
    % test_gn.agents(1).goal_confs = P1_init_goal_confs(:,ii);
    % test_gn.agents(2).goal_confs = zeros(size(test_gn.agents(2).goal_confs));
    % mult_factor = 2.0;
    % lie_idxs = test_gn.all_self_utils(:, 2) > 0;
    % test_gn.all_self_utils(lie_idxs, 2) = test_gn.all_self_utils(lie_idxs, 2) * mult_factor;
    % test_gn.all_costs(lie_idxs, 2) = test_gn.init_costs(2) - test_gn.all_self_utils(lie_idxs, 2);

    test_gn.ComputeAllUtils();
    test_gn.RunSequentialNegotiation();
    for id = 1:num_players
        table_num = sum(test_gn.stg_table_bool_hist{end, id}, 2);
        final_N_tasks(ii, id) = sum(table_num - ones(size(table_num)));
    end
    if isempty(test_gn.stg_offer_hist)
        final_neg_rnds(ii) = 0;
        final_trusts(ii, :) = [test_gn.agents(1).trusts(2), test_gn.agents(2).trusts(1)];
        final_utils(ii, :) = [0, 0];
        final_self_utils(ii, :) = [0, 0];
        continue;
    end

    deal_idx = test_gn.stg_offer_hist{end}(end, :);
    trust_hist = test_gn.stg_trust_hist;
    final_deal_idx(ii) = deal_idx(1);
    if (diff(deal_idx))
        [~,final_deal_idx(ii)] = ismember(final_self_utils(ii, :), test_gn.all_utils, "row");
    end
    final_neg_rnds(ii) = size(trust_hist, 3)-1;
    final_utils(ii, :) = test_gn.all_utils(final_deal_idx(ii), :);
    final_self_utils(ii, :) = test_gn.all_self_utils(final_deal_idx(ii), :);
    final_trusts(ii, :) = [trust_hist(2,1,end), trust_hist(1,2,end)];
end
final_num_tasks = sum(final_N_tasks, 2);
final_deal_idx(final_deal_idx == 0) = conflict_offer_idx;

% Monte-Carlo visualizer
figure()
subplot(2,1,1)
histogram(final_self_utils(:,1))
ylabel("$\tilde{u}^{(1)}$", "Interpreter", "latex")
title(strcat("mean: ", num2str(mean(final_self_utils(:,1)))))
grid on
subplot(2,1,2)
histogram(ref_self_utils(final_deal_idx,2))
ylabel("$\tilde{u}^{(2)}$", "Interpreter", "latex")
title(strcat("mean: ", num2str(mean(ref_self_utils(final_deal_idx,2)))))
grid on

figure()
subplot(2,1,1)
histogram(final_trusts(:,1), 0:0.1:1)
grid on
ylabel("$\rho^{(1)}_2$", "Interpreter", "latex")
title(strcat("mean: ", num2str(mean(final_trusts(:,1)))))
subplot(2,1,2)
histogram(final_trusts(:,2), 0:0.1:1)
title(strcat("mean: ", num2str(mean(final_trusts(:,2)))))
ylabel("$\rho^{(2)}_1$", "Interpreter", "latex")
grid on

figure()
histogram(final_neg_rnds)
ylabel("# of rounds")
title(strcat("mean: ", num2str(mean(final_neg_rnds))))
grid on

figure()
histogram(final_num_tasks)
ylabel("# of table tasks")
title(strcat("mean: ", num2str(mean(final_num_tasks))))
grid on

figure()
plot(1:num_monte, final_deal_idx, ".");
ylabel("deal index")
grid on

figure()
histogram(final_deal_idx, 1:size(goal_negotiator.all_deals_idxs, 1))
grid on

figure()
plot(ref_self_utils(:, 1), ref_self_utils(:, 2), ".")
hold on
final_deal_idx(final_deal_idx == 0) = conflict_offer_idx;
copy_final_deal_idx = final_deal_idx;
color_code_list = ["ro", "go", "bo", "co", "mo"];
counting_list = ["1st", "2nd", "3rd", "4th", "5th"];
legend_list = ["All deals", "Conflict"];
peak_modes = [];
plot(0, 0, "xk")
for jj = 1:numel(color_code_list)
    if ~sum(~isnan(copy_final_deal_idx))
        break;
    end
    [peak_modes(jj), freq_peak] = mode(copy_final_deal_idx);
    plot(ref_self_utils(peak_modes(jj), 1), ref_self_utils(peak_modes(jj), 2), color_code_list(jj))
    copy_final_deal_idx(copy_final_deal_idx == peak_modes(jj)) = NaN;
    legend_list = cat(2, legend_list, strcat(counting_list(jj), " mode (", num2str(freq_peak), " / ", num2str(num_monte), ")"));
end
hold off
grid on
legend(legend_list)
xlabel("$\tilde{u}^{(1)}$", "Interpreter", "latex")
ylabel("$\tilde{u}^{(2)}$", "Interpreter", "latex")
