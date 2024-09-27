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

%% Negotiation
% Choose one of the presets
% PRESET: P2 is cooperative [Example V.B.1]
% for id = 1:num_players
%     goal_negotiator.agents(id).deception_bool = false;
%     goal_negotiator.agents(id).trusts = 0.5 * ones(num_players, 1);
%     goal_negotiator.agents(id).trusts(id) = 1;
% end

% PRESET: P2 is selfish, constant inflation on self-util of (all, nothing) [Example V.B.2.+]
% goal_negotiator.agents(1).deception_bool = false;
% goal_negotiator.agents(2).deception_bool = true;
% goal_negotiator.agents(1).trusts(2) = 0.5;
% goal_negotiator.agents(2).trusts(1) = 0.0;
% goal_negotiator.agents(2).goal_confs = zeros(size(goal_negotiator.agents(2).goal_confs));
% add_factor = 10;
% [~,lie_idx] = max(goal_negotiator.all_self_utils(:, 2));
% goal_negotiator.all_self_utils(lie_idx,2) = goal_negotiator.all_self_utils(lie_idx,2) + add_factor;
% goal_negotiator.all_costs(lie_idx,2) = goal_negotiator.all_costs(lie_idx,2) - add_factor;
% goal_negotiator.ComputeAllUtils();

% PRESET: P2 is selfish, double inflation on positive self-util deals [Example V.B.2.x]
goal_negotiator.agents(1).deception_bool = false;
goal_negotiator.agents(2).deception_bool = true;
goal_negotiator.agents(1).trusts(2) = 0.5;
goal_negotiator.agents(2).trusts(1) = 0.0;
goal_negotiator.agents(2).goal_confs = zeros(size(goal_negotiator.agents(2).goal_confs));
mult_factor = 2.0;
lie_idxs = goal_negotiator.all_self_utils(:, 2) > 0;
goal_negotiator.all_self_utils(lie_idxs, 2) = goal_negotiator.all_self_utils(lie_idxs, 2) * mult_factor;
goal_negotiator.all_costs(lie_idxs, 2) = goal_negotiator.init_costs(2) - goal_negotiator.all_self_utils(lie_idxs, 2);
goal_negotiator.ComputeAllUtils();

% PRESET: P2 is selfish, constant inflation, P1 has no anti-deception measures [Example V.C.+]
% goal_negotiator.agents(1).deception_bool = true;
% goal_negotiator.agents(2).deception_bool = true;
% goal_negotiator.agents(1).trusts(2) = 0.5;
% goal_negotiator.agents(2).trusts(1) = 0.0;
% goal_negotiator.agents(2).goal_confs = zeros(size(goal_negotiator.agents(2).goal_confs));
% add_factor = 10;
% [~,lie_idx] = max(goal_negotiator.all_self_utils(:, 2));
% goal_negotiator.all_self_utils(lie_idx,2) = goal_negotiator.all_self_utils(lie_idx,2) + add_factor;
% goal_negotiator.all_costs(lie_idx,2) = goal_negotiator.all_costs(lie_idx,2) - add_factor;
% goal_negotiator.ComputeAllUtils();

% PRESET: P2 is selfish, double inflation, P1 has no anti-deception measures [Example V.C.x]
% goal_negotiator.agents(1).deception_bool = true;
% goal_negotiator.agents(2).deception_bool = true;
% goal_negotiator.agents(1).trusts(2) = 0.5;
% goal_negotiator.agents(2).trusts(1) = 0.0;
% goal_negotiator.agents(2).goal_confs = zeros(size(goal_negotiator.agents(2).goal_confs));
% mult_factor = 2.0;
% lie_idxs = goal_negotiator.all_self_utils(:, 2) > 0;
% goal_negotiator.all_self_utils(lie_idxs, 2) = goal_negotiator.all_self_utils(lie_idxs, 2) * mult_factor;
% goal_negotiator.all_costs(lie_idxs, 2) = goal_negotiator.init_costs(2) - goal_negotiator.all_self_utils(lie_idxs, 2);
% goal_negotiator.ComputeAllUtils();

goal_negotiator.RunSequentialNegotiation();
for id = 1:num_players
    fleet_manager.agents(id).goals = goal_negotiator.agents(id).goals;
end

%% Negotiation results
figure()
color_markers = ["r", "b"];
hold on
for id = 1:num_players
    goal_pos = goal_negotiator.agents(id).goals;
    plot(goal_pos(:,1), goal_pos(:,2), "o", "MarkerEdgeColor", color_markers(id), "MarkerSize", 8)
end
for id = 1:num_players
    goal_pos = goal_sets{id};
    goal_confs = goal_conf_sets{id};
    plot(goal_pos(:,1), goal_pos(:,2), ".", "MarkerEdgeColor", color_markers(id), "MarkerSize", 8)
    text(goal_pos(:,1), goal_pos(:,2)-0.25, num2str(goal_confs, "%.3f"), "HorizontalAlignment", "center");
    for ij = 1:numel(goal_confs)
        text(goal_pos(ij,1), goal_pos(ij,2)+0.4, strcat("$T^{(", num2str(id), ")}_{", num2str(ij), "}$"), "HorizontalAlignment", "center", "Interpreter", "latex");
    end
end
for id = 1:num_players
    text(player_pos(id, 1), player_pos(id, 2), num2str(id), "HorizontalAlignment", "center", "FontWeight","bold");
end
hold off
title("Global Map");
ax = gca;
ax.XTick = 0:1:10;
ax.YTick = 0:1:10;
grid on;
axis equal;
xlim([0, 10]);
ylim([0, 10]);
legend_txt = ["$\mathcal{T}^{(1)}_f$",...
            "$\mathcal{T}^{(2)}_f$",...
            "$\mathcal{T}^{(1)}_0$", "$\mathcal{T}^{(2)}_0$"];
legend(legend_txt, "Location", "best", "Interpreter", "latex")
xlabel("X [m]")
ylabel("Y [m]")

agrd_idx = goal_negotiator.stg_offer_hist{end}(end, 1);
marker_size = 8;
figure()
plot(ref_self_utils(:, 1), ref_self_utils(:, 2), ".")
hold on
plot(ref_self_utils(agrd_idx,1), ref_self_utils(agrd_idx,2), "bo", "MarkerSize", marker_size)
plot(0, 0, ".k", "MarkerSize", marker_size)
[~, max_nash_idx] = max(prod(ref_self_utils, 2));
plot(ref_self_utils(max_nash_idx, 1), ...
    ref_self_utils(max_nash_idx, 2), "r.", "MarkerSize", marker_size)
hold off
grid on
legend("All deals", "Negotiation", "Conflict", "Max Social Welfare")
xlabel("$\tilde{u}^{(1)}$", "Interpreter", "latex", "FontSize", 14)
ylabel("$\tilde{u}^{(2)}$", "Interpreter", "latex", "FontSize", 14)
axis equal
