classdef myTime < handle
    % Class for time measurement and cycle time 
    
    properties
        name
        time_online
        time_offline
        time_step
        handle_online
        handle_offline
        handle_step
    end
    
    methods
        function obj=myTime(name)
            % Constructor & start of measurement:
            % Inputs:
            % name: (string) tag for time measurement
            obj.handle_offline = tic;
            obj.handle_online = tic;
            obj.handle_step = tic;
            if exist('name','var')
                obj.name = name;
            else
                obj.name =[];
            end
        end
        
        function offline(obj)
            %% Offline Time
            % give time since constructor and this function
            fprintf('%s Offline time: ',obj.name)
            toc(obj.handle_offline) % print
            obj.time_offline = toc(obj.handle_offline);
            obj.handle_online =tic;
            obj.handle_step = tic;
        end
        function online(obj)
            %% Online Time
            % give the online time and the average cycle time
            fprintf('\n')
            fprintf('%s Online time: ',obj.name)
            toc(obj.handle_online)
            obj.time_online = toc(obj.handle_online);
            
            aver = mean(obj.time_step);
            fprintf('Average time per iteration:   %f seconds.\n', aver);
        end
        function step(obj)
            % measure point within the online function
            obj.time_step = [obj.time_step, toc(obj.handle_step)];
            obj.handle_step = tic;
        end
    end
    
end

