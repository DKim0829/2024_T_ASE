classdef FleetManager < matlab.mixin.Copyable
    properties (Access = public)
        num_players {mustBeInteger}
        goal_pos {mustBeNumeric}
        player_pos {mustBeNumeric}
        agents
    end
    properties (Access = private)
        world_size_
    end
    
    methods
        % constructor
        function self = FleetManager(agents, world_size, goal_pos, player_pos)
            self.agents = agents;
            self.num_players = numel(agents);
            self.world_size_ = world_size;
            self.goal_pos = goal_pos;
            self.player_pos = player_pos;
        end

        % setter
        function set.world_size_(self, world_size)
            self.world_size_ = world_size;
        end
        
        % getter
        function world_size = get.world_size_(self)
            world_size = self.world_size_;
        end
    end

    methods (Access = public)
        function self = SetAgent(self, agent, id)
            self.agents(id) = agent;
        end
        
        function agent = GetAgent(self, id)
            agent = self.agents(id);
        end

        function self = DrawGlobalMap(self, fid)
            if (nargin < 2)
                figure()
            else
                figure(fid)
            end
            rem_goals = self.goal_pos;
            legend_txt = [];
            color_code = ["r", "b", "g", "y", "c", "m"];
            for id = 1:self.num_players
                agent = self.GetAgent(id);
                plot(agent.goals(:,1), agent.goals(:,2), strcat("*", color_code(id)));
                hold on
                legend_txt = cat(2, legend_txt, strcat("$\mathcal{T}^{(", num2str(id), ")}_0$"));
                rem_goals = rem_goals(~ismember(rem_goals, self.agents(id).goals, "rows"), :);
            end
            if ~isempty(rem_goals)
                legend_txt = cat(2, legend_txt, "Unassigned");
                plot(self.goal_pos(:,1), self.goal_pos(:,2), "*k")
            end
            for id = 1:self.num_players
                agent = self.GetAgent(id);
                text(agent.pos(1), agent.pos(2), num2str(id), "HorizontalAlignment", "center");
            end
            hold off;
            title("Global Map");
            ax = gca;
            ax.XTick = 0:1:self.world_size_;
            ax.YTick = 0:1:self.world_size_;
            grid on;
            axis equal;
            xlim([0, self.world_size_]);
            ylim([0, self.world_size_]);
            legend(legend_txt, "Location", "best", "Interpreter", "latex")
            xlabel("X [m]")
            ylabel("Y [m]")
        end
    end
end
