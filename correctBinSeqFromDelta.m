function fracIncorrect = correctBinSeqFromDelta(params, nets, seqTrials, pL, alpha, beta, eps, nReps, maxsteps, h, thres)
    for i = 1:size(nets, 2)
        N(i) = nets{i}.numnodes;
    end

    for q = 1:nReps
        % initial oscillator conditions
        omegas{1} = rand(1, N(1));
        omegas{2} = rand(1, N(2));
        thetas0{1} = 2 * pi * rand(1, N(1));
        thetas0{2} = 2 * pi * rand(1, N(2));
        lams0 = [0, 0];
        trueSeq{q} = GenRandSeq(seqTrials, pL);
        [~, ~, ~, ~, decision(q, :, :, :)] = pooledInhibBinarySequential(nets, trueSeq{q}, thres, alpha, beta, params(1), params(2), params(3), eps, omegas, thetas0, lams0, maxsteps, h);
        corrects(q, :) = sum(squeeze(decision(q, :, :, 1)), 1) == trueSeq{q};
    end
    fracIncorrect = 1 - sum(sum(corrects)) / (nReps * length(trueSeq{q}));
end