function environment = createEnvironment(numRewards, maxReward, rewardDistFunc, nodeDropout, randomness)
  % Given the reward and randomness parameters, createEnvironment generates a struct 
  % representing the initial conditions for an environment.
  %
  % rewardDistFunc: uniformInteger, uniformCont, gaussianInteger, gaussianCont
  % randomness: none 0, low (0,0.15], medium (0.15-0.3], high (0.3-1]
   
  %% Populate the rewards vector according to numRewards, maxReward and the rewardDistFunc 
  switch rewardDistFunc
    case 'uniformInteger'
      rewards = randi([1 maxReward],numRewards,1);
    case 'uniformCont'
      rewards = (maxReward - 1) * rand([numRewards 1]) + 1;
    case 'gaussianInteger'
      disp('gaussianInteger is not implemented. You can help by coding this yourself.');
      keyboard
    case 'gaussianContinuous'
      disp('gaussianContinuous is not implemented. You can help by coding this yourself.');
      keyboard
    otherwise
      disp(['ERROR: Unrecognized rewardDistFunc term ' rewardDistFunc '. Should be uniformInteger, uniformCont, gaussianInteger, gaussianCont.'])
      keyboard
  end
  
  %% Pick an appropriate value for the environment's stochasticity. This will be 
  % used as noise during value perception.
  switch randomness
    case 'none'
      environmentStochasticity = 0;
    case 'low'
      environmentStochasticity = rand() * 0.15;
    case 'medium'
      environmentStochasticity = 0.15 + (rand() * 0.15);
    case 'high'
      environmentStochasticity = 0.3 + (rand() * 0.7);
    otherwise
      disp(['ERROR: Unrecognized randomness term ' randomness '. Should be none, low, medium or high.'])
      keyboard
  end
  
  %% Ensure that the node dropout rate is in the interval [0,1).
  if nodeDropout >= 1 || nodeDropout < 0
    disp('nodeDropout must be in the interval [0,1).');
    keyboard
  end
  
  %% global urgency starts at 0.
  globalUrgency = 0;
  
  %% At some point, initialize the global urgency update equation.
  %% Or initialize the global urgency equation's coefficients.
  
  %% Create the struct to be returned.
  environment = struct('rewards', rewards, 'stochasticity', environmentStochasticity, 
  'nodeDropRate', nodeDropout, 'urgency', globalUrgency);
  
end