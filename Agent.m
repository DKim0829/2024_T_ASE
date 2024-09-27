classdef Agent < matlab.mixin.Copyable
    properties (Access = public)
        id {mustBeInteger}
        pos {mustBeNumeric}
        goals {mustBeNumeric}
        goal_confs {mustBeNumeric}
        trusts {mustBeNumeric}
        dist_table {mustBeNumeric}
        deception_bool {mustBeNumericOrLogical}
    end

    properties (Access = private)
        num_players_ {mustBeInteger}
        world_size_ {mustBeNumeric}
    end
    
    methods
        % constructor
        function self = Agent(id, num_players, world_size, pos, goals)
            self.id = id;
            self.num_players_ = num_players;
            self.world_size_ = world_size;
            self.pos = pos;
            self.goals = goals;

            self.goal_confs = zeros(size(goals, 1), 1);
            self.trusts = 1.0 * ones(self.num_players_, 1);
            self.deception_bool = false;

            self.UpdateDistTable();
        end

        % setter
        function set.world_size_(self, world_size)
            self.world_size_ = world_size;
        end
        function set.goals(self, goals)
            self.goals = goals;
            self.UpdateDistTable();
        end
        % getter
        function world_size = get.world_size_(self)
            world_size = self.world_size_;
        end
    end

    methods (Access = public)
        function self = UpdateDistTable(self)
            num_goals = size(self.goals, 1);
            dist_table = zeros(num_goals+1, num_goals+1);
            for i = 1:num_goals
                for j = i+1:num_goals
                    to_goal = self.goals(i, :);
                    from_goal = self.goals(j, :);
                    dist_table(i, j) = norm(from_goal - to_goal, 2);
                end
            end
            for i = 1:num_goals
                dist_table(end, i) = norm(self.pos - self.goals(i, :), 2);
            end
            self.dist_table = dist_table + dist_table.';
        end
        
        function bools = GetTransTask(self, trans_thres)
            bools = (self.goal_confs <= trans_thres);
        end

        function bools = GetSecretTask(self, secret_thres)
            bools = (self.goal_confs >= secret_thres);
        end

        function bools = GetTrustTask(self)
            bools = zeros(numel(self.goal_confs), self.num_players_, "logical");
            for ii = 1:self.num_players_
                bools(:, ii) = (self.goal_confs <= self.trusts(ii));
            end
        end

        function util = ComputeUtil(self, self_utils)
            util = self.trusts.' * self_utils;
        end

        function cncsn_limit = ConstructConcessionLimit(self, stg_self_utils, stg_utils, NS_bools)
            cncsn_limit = 0;
            if (self.deception_bool || ~sum(NS_bools))
                return;
            end
            [~, min_idx] = min(stg_self_utils(:, self.id));
            % min_util = max(stg_utils(min_idx, self.id), 0);
            min_util = stg_utils(min_idx, self.id);
            max_NS_util = max(stg_utils(NS_bools, self.id));
            cncsn_limit = 0.5 * (1 - min(self.trusts)) * (max_NS_util - min_util) + min_util;
        end
    end
end
