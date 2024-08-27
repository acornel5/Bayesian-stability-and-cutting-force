% manualTestClassifier.m
% Cornelius A.

% Prompts the user to classify the test as stable/unstable using whatever
% metric they prefer.

function test = manualTestClassifier( ...
    iteration, ...
    test, ...
    param)
    
    while true
        disp("A test is recommended at: ")
        disp(append("Spindle speed: ", num2str(test.RPM), " rpm"))
        disp(append("Axial depth: ", num2str(test.b*1000), " mm"))
        disp(append("Radial width: ", num2str(test.a*1000), " mm"))
        disp(append("Feedrate: ", num2str(test.fz*1000), " mm/flute, ", num2str(test.fz*param.FluteCount*test.RPM*1000), " mm/min"))

        broken = input("Did the endmill break? Enter 1 for yes and 0 for no: ");
        if broken == 1
            disp(append("You have said that the tool broke."))
            if input("Enter 1 if this is correct: ") == 1
                stable = 0;
                chatterFrequency = 0;
                break
            end
        end

        stable = input("Was this test stable? Enter 1 for yes and 0 for no: ");
        chatterFrequency = 0;
        if stable == 0 || stable == 0.5
            chatterFrequency = input("What was the dominant chatter frequency? ");
        end

        power = input("What was the spindle power in watts? ");
        
        if stable == 1
            disp(append("You have said that the test was stable."))
		elseif stable == 0.5 % Marginal
			disp(append("You have said that the test was mariginal, with a dominant frequency of ", num2str(chatterFrequency, 5), " Hz."))
        elseif stable == 0
            disp(append("You have said that the test chattered, with a dominant frequency of ", num2str(chatterFrequency, 5), " Hz."))
        end
        if power ~= 0, disp(strcat("You have said that the cutting power was ", num2str(power), " W.")), end
        if input("Enter 1 if this is correct: ") == 1
            break
        end
    end

    test.Broken = broken;
    test.Stable = stable;
    test.ChatterFrequency = chatterFrequency;
    test.Power = power;
end