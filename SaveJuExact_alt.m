function JuExact = SaveJuExact(problem,GOAL,pdeSolver);
  if strcmp(problem,'OptimalDesign_Lshape')
  % OD,Lshape
    if strcmp(GOAL,'GOAL_IM')
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.0563419333818107; %uniform 800 AitkenExtrapolation DWR-P0-RT0, IM
      else
        JuExact = 0.057215911128541; %uniform 800 AitkenExtrapolation DWR-P1, IM
      end
    elseif strcmp(GOAL,'GOAL_DP')
      JuExact = 0.0589777082433764; %uniform 2000 AitkenExtrapolation DWR-P1, DP
    else
      JuExact = 0.0589777082433764; %uniform 2000 AitkenExtrapolation DWR-P1, PR
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(problem,'OptimalDesign_Lshape_exact')
  % OD,Lshape exact
    if strcmp(GOAL,'GOAL_IM')
      JuExact = 0.395925305770294; %uniform 2000 AitkenExtrapolation DWR-P1, IM
    elseif strcmp(GOAL,'GOAL_DP')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1, DP
    else
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1, PR
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(problem,'OptimalDesign_SquareSlit')
  % OD,SquareSlit
    if strcmp(GOAL,'GOAL_IM')
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.112963697007926; %uniform 3200 AitkenExtrapolation DWR-P0-RT0, IM
      else
        JuExact = 0.119416556216333; %uniform 2400 AitkenExtrapolation DWR-P1, IM
      end
    elseif strcmp(GOAL,'GOAL_DP')
      JuExact = 0.0749854840150314; %uniform 2400 AitkenExtrapolation DWR-P1, DP
    else
      JuExact = 0.0749854840150314; %uniform 2400 AitkenExtrapolation DWR-P1, PR
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(problem,'OptimalDesign_SquareSlit_exact')
  % OD,SquareSlit exact
    if strcmp(GOAL,'GOAL_IM')
      if strcmp(pdeSolver,'ODRT')
        JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P0-RT0, IM
      else
        JuExact = -0.112061002008049; %uniform 1000 AitkenExtrapolation DWR-P1, IM
      end
    elseif strcmp(GOAL,'GOAL_DP')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    else
      JuExact = 0.427545087815354; %uniform 1000 AitkenExtrapolation DWR-P1, PR
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(problem,'TwoWell_Square')
  % TW,Square
    if strcmp(GOAL,'GOAL_IM')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    elseif strcmp(GOAL,'GOAL_DP')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    else
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  elseif strcmp(problem,'TwoWell_SquareExact')
  % TW,SquareExact
    if strcmp(GOAL,'GOAL_IM')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    elseif strcmp(GOAL,'GOAL_DP')
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
    else
      JuExact = 0.1; %uniform 00 AitkenExtrapolation DWR-P1
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