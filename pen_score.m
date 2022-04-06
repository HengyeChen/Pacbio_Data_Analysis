function penalty = pen_score(distance,wfactor,linker)
%% this function calculate the penalty score from short linker

if distance > linker% when distance is longer than 20bp, no penalty
    penalty = 0;
else%when distance is negative, high penalty
    penalty = (linker-distance)/10;
end
penalty = penalty*wfactor;
