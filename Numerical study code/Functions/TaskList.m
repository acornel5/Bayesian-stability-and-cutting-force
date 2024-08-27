classdef TaskList < handle
    % A list of the the current tasks that the program is working on
    
    properties
        Tasks
        Enabled
        ETAFunc
		LogStream
		ETATimer
    end
    
    methods
        function obj = TaskList()
            obj.Tasks = {};
            obj.Enabled = 1;
			obj.LogStream = fopen('Figures//Log.txt', 'w');
        end
        
        function obj = Add(obj, task)
            if obj.Enabled
                obj.Tasks{end+1} = task;
                obj.DisplayTasks;
                obj.Log(task);
            end
        end
        
        function obj = Remove(obj)
            if obj.Enabled
                removedTask = string(obj.Tasks(end));
                obj.Tasks(end) = [];
                obj.DisplayTasks;
                disp(append(removedTask, " complete"))
                obj.Log(append(removedTask, " complete"));
            end
        end
        
        function obj = Update(obj, newTask)
            if obj.Enabled
                if ~isempty(obj.Tasks)
                    obj.Tasks{end} = newTask;
                    obj.DisplayTasks;
                else
                    obj.Add(newTask)
                end
                obj.Log(newTask);
            end
        end
        
        function DisplayTasks(obj)
            if obj.Enabled
                clc
                for i = 1:length(obj.Tasks)
                    disp(obj.Tasks{i})
                end
            end
        end
        
        function InitializeETA(obj, numberToGenerate)
            [obj.ETAFunc, obj.ETATimer] = GenerateETAFunction(numberToGenerate);
            obj.Log(append("Beginning ETA with length ", num2str(numberToGenerate)));
        end
        
        function InitializeCustomETA(obj, customETAFunc)
            obj.ETAFunc = customETAFunc;
            obj.Log("Beginning custom ETA");
        end
        
        function StepETA(obj, varargin)
            if obj.Enabled
                if nargin == 0
                    send(obj.ETAFunc, 1);
					currentTime = obj.ETATimer();
                else
                    send(obj.ETAFunc, varargin);
					currentTime = obj.ETATimer();
                end
            end
        end
        
        function ClearETA(obj)
            obj.ETAFunc = [];
			obj.ETATimer = [];
            obj.Log("ETA completed");
        end
        
        function Clear(obj)
            clc
            obj.Tasks = [];
        end
        
        function Disable(obj)
            clc
            obj.Enabled = false;
        end
        
        function Enable(obj)
            clc
            obj.Enabled = true;
        end

        function Log(obj, str)
            try % Tries to log, but if it can't then it just skips it
                fprintf(obj.LogStream, str);
                fprintf(obj.LogStream, '\n');
            catch
            end
        end
    end
end

