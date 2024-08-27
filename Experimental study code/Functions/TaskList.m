classdef TaskList < handle
    % A list of the the current tasks that the program is working on
    
    properties
        Tasks
        Enabled
        ETAFunc
    end
    
    methods
        function obj = TaskList()
            obj.Tasks = {};
            obj.Enabled = 1;
        end
        
        function obj = Add(obj, task)
            if obj.Enabled
                obj.Tasks{end+1} = task;
                obj.DisplayTasks;
            end
        end
        
        function obj = Remove(obj)
            if obj.Enabled
                removedTask = string(obj.Tasks(end));
                obj.Tasks(end) = [];
                obj.DisplayTasks;
                disp(append(removedTask, " complete"))
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
            obj.ETAFunc = GenerateETAFunction(numberToGenerate);
        end
        
        function InitializeCustomETA(obj, customETAFunc)
            obj.ETAFunc = customETAFunc;
        end
        
        function StepETA(obj, varargin)
            if obj.Enabled
                if nargin == 0
                    send(obj.ETAFunc, 1);
                else
                    send(obj.ETAFunc, varargin);
                end
            end
        end
        
        function ClearETA(obj)
            obj.ETAFunc = [];
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
    end
end

