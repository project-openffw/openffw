function JuExact = SaveJuExact(problem,GOAL,pdeSolver);
  if strcmp(problem,'OptimalDesign_Lshape')
  % OD,Lshape
    if strcmp(GOAL,'GOAL_IM')
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P0-RT0, IM
      else
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1, IM
      end
    elseif strcmp(GOAL,'GOAL_DP')
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P0-RT0, DP
      else
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1, DP
      end
    else
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P0-RT0, PR
      else
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1, PR
      end
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(problem,'OptimalDesign_Lshape_exact')
  % OD,Lshape exact
    if strcmp(GOAL,'GOAL_IM')
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P0-RT0, IM
      else
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1, IM
      end
    elseif strcmp(GOAL,'GOAL_DP')
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P0-RT0, DP
      else
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1, DP
      end
    else
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P0-RT0, PR
      else
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1, PR
      end
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(problem,'OptimalDesign_SquareSlit')
  % OD,SquareSlit
    if strcmp(GOAL,'GOAL_IM')
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.0353124120021308; %uniform 3256 AitkenExtrapolation DWR-P0-RT0, IM
      else
        JuExact = 0.1; %uniform 2449 AitkenExtrapolation DWR-P1, IM
      end
    elseif strcmp(GOAL,'GOAL_DP')
      if strcmp(pdeSolver,'ODRT')
        JuExact = -0.895264544837896; %uniform 00 AitkenExtrapolation DWR-P0-RT0, DP
      else
        JuExact = 0.1; %falsch %uniform 2449 AitkenExtrapolation DWR-P1, DP
      end
    else
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.329480291755864; %uniform 00 AitkenExtrapolation DWR-P0-RT0, PR
      else
        JuExact = 0.1; %uniform 2449 AitkenExtrapolation DWR-P1, PR
      end
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(problem,'OptimalDesign_SquareSlit_exact')
  % OD,SquareSlit exact
    if strcmp(GOAL,'GOAL_IM')
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P0-RT0, IM
      else
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1, IM
      end
    elseif strcmp(GOAL,'GOAL_DP')
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P0-RT0, DP
      else
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1, DP
      end
    else
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P0-RT0, PR
      else
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1, PR
      end
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(problem,'TwoWell_Square')
  % TW,Square
    if strcmp(GOAL,'GOAL_IM')
      JuExact = 0.0336997140148101; %uniform 1457 AitkenExtrapolation DWR-P1
    elseif strcmp(GOAL,'GOAL_DP')
      JuExact = -0.60909171428965; %uniform 1457 AitkenExtrapolation DWR-P1
    else
      JuExact = 40.7114749579159;% -0.197445775615336%falsch %uniform 1457 AitkenExtrapolation DWR-P1
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(problem,'TwoWell_SquareExact')
  % TW,SquareExact
    if strcmp(GOAL,'GOAL_IM')
      JuExact = 0.245533802871385; %uniform 1457 AitkenExtrapolation DWR-P1
    elseif strcmp(GOAL,'GOAL_DP')
      JuExact = 3.80526517869181; %uniform 1457 AitkenExtrapolation DWR-P1
    else
      JuExact = 6.08328785969069; %1.81245206523277;%falsch %uniform 1457 AitkenExtrapolation DWR-P1
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(problem,'AbsValue_Square')
  % AV,Square
    if strcmp(GOAL,'GOAL_IM')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    elseif strcmp(GOAL,'GOAL_DP')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    else
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(problem,'AbsValue_Lshape')
  % AV,Lshape
    if strcmp(GOAL,'GOAL_IM')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    elseif strcmp(GOAL,'GOAL_DP')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    else
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(problem,'AbsValue_Square_exact')
  % AV,Square exact
    if strcmp(GOAL,'GOAL_IM')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    elseif strcmp(GOAL,'GOAL_DP')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    else
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  else %AbsValue_Lshape_exact
  % AV,Lshape exact
    if strcmp(GOAL,'GOAL_IM')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    elseif strcmp(GOAL,'GOAL_DP')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    else
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    end
  end