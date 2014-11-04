# (feas_point) tests whether the given (k,lambda) pair is feasible.
# The returned value is the Mosek solution code
# (type Enum solsta, e.g. MSK_SOL_STA_x)
# 1  = feasible (MSK_SOL_STA_OPTIMAL)
# 8  = nearly feasible (MSK_SOL_STA_NEAR_OPTIMAL)
# 12 = nearly infeasible (MSK_SOL_STA_NEAR_PRIM_INFEAS_CER)
# 5  = infeasible (MSK_SOL_STA_PRIM_INFEAS_CER)

function feas_point( k, lambda )

    # algorithm parameters
    bet = (sqrt(k)-1)/(sqrt(k)+1);
    xt = [ 1 0 0 0
           0 1 0 0
           0 0 1 0 ]
    xtt =  [ 0  -bet   0  bet
             1  1+bet  0  -1-bet
             0    1    0  -1 ]
    q1 = [ 0  1    0 -1
           0 -1/k  0  1 ]

    q2 = [ 0  1   -1 -1
           0 -1/k  0  1 ]

    J = [ 0 1; 1 0]

    tsk = maketask()

    # enable multi-threading
    putintparam(tsk,MSK_IPAR_NUM_THREADS,0) 

    # INTERIOR-POINT TOLERANCES
    # http://docs.mosek.com/7.0/toolbox/The_optimizers_for_continuous_problems.html
    # Primal feasibility tolerance used by the conic interior-point optimizer.
    putdouparam(tsk,MSK_DPAR_INTPNT_CO_TOL_PFEAS,1e-8)    # default: 1e-8   IMPORTANT
    # Dual feasibility tolerance used by the conic interior-point optimizer.
    putdouparam(tsk,MSK_DPAR_INTPNT_CO_TOL_DFEAS,1e-8)    # default: 1e-8   IMPORTANT
    # Relative gap termination tolerance used by the conic interior-point optimizer.
    putdouparam(tsk,MSK_DPAR_INTPNT_CO_TOL_REL_GAP,1e-7)  # default: 1e-7   IMPORTANT
    # Nonlinear solver infeasibility tolerance parameter.
    putdouparam(tsk,MSK_DPAR_INTPNT_TOL_INFEAS,1e-10)     # default: 1e-10
    # Infeasibility tolerance for the conic solver.
    putdouparam(tsk,MSK_DPAR_INTPNT_CO_TOL_INFEAS,1e-10)  # default: 1e-10
    # Optimality tolerance for the conic solver.
    putdouparam(tsk,MSK_DPAR_INTPNT_CO_TOL_MU_RED,1e-8)   # default: 1e-8
    # Optimality tolerance for the conic solver.
    putdouparam(tsk,MSK_DPAR_INTPNT_CO_TOL_NEAR_REL,1e3)  # default: 1e3

    # variables t1, t2
    appendvars(tsk,2)
    putvarname(tsk,1,"t1")
    putvarname(tsk,2,"t2")

    # variables P(3,3) and Y(4,4)
    appendbarvars(tsk,[3 4])
    putbarvarname(tsk,1,"P")
    putbarvarname(tsk,2,"Y")

    # constrain 0 <= ti <= inf
    for j = 1:2
        putvarbound(tsk,j,MSK_BK_LO,0.,Inf)
    end

    # declare constraints (10 semidef, 1 rho-hard, 1 homog-breaking)
    appendcons(tsk,12)

    # linear constraints involving P and Y
    ix =  [1, 2, 3, 2, 3, 3]
    jx =  [1, 1, 1, 2, 2, 3]
    iix = [1, 2, 3, 4, 2, 3, 4, 3, 4, 4]
    jix = [1, 1, 1, 1, 2, 2, 2, 3, 3, 4]
    Pcon = zeros(3,3,10)
    Ycon = zeros(4,4,10)
    for c = 1:10
        Pcon[:,:,c] = xtt[:,iix[c]]*xtt[:,jix[c]]' - lambda*xt[:,iix[c]]*xt[:,jix[c]]'
        Pcon[:,:,c] = (Pcon[:,:,c]+Pcon[:,:,c]')/2
        symaP = appendsparsesymmat(tsk,3,ix,jx,diag(Pcon[ix,jx,c]))
        val = iix[c] == jix[c] ? 1.0 : 0.5
        symaY = appendsparsesymmat(tsk,4,[iix[c]],[jix[c]],[val])
        putbaraij(tsk,c,1, [symaP], [1.0])
        putbaraij(tsk,c,2, [symaY], [1.0])
        putconbound(tsk,c,MSK_BK_FX,0.0,0.0)
        putconname(tsk,c,string("Y(", iix[c], ",", jix[c], ")"))
    end

    # linear constraints involving ti
    Acon = zeros(4,4,2)
    Acon[:,:,1] = q1'*J*q1
    Acon[:,:,2] = q2'*J*q2
    for j = 1:2
        putacol(tsk,j,[1:10],diag(Acon[iix,jix,j]))
    end

    # rho-hardness constraint
    putarow(tsk,11,[1 2],[lambda lambda-1.])
    putconbound(tsk,11,MSK_BK_LO,0.,Inf)
    putconname(tsk,11,"rho-hardness")

    # homogeneity-breaking constraint, t2 == 1
    putarow(tsk,12,[2],[1.0])
    putconbound(tsk,12,MSK_BK_FX,1.,1.)
    putconname(tsk,12,"t2=const")

    optimize(tsk)
    solsta = getsolsta(tsk,MSK_SOL_ITR)

    # if infeasible, return code
    if solsta != MSK_SOL_STA_OPTIMAL
        return int32(solsta)
    else
        t = getxx(tsk,MSK_SOL_ITR)
        p = getbarxj(tsk,MSK_SOL_ITR,1)
        y = getbarxj(tsk,MSK_SOL_ITR,2)
        P = [ p[1] p[2] p[3]
              p[2] p[4] p[5]
              p[3] p[5] p[6] ]
        return P, t[1], t[2]
    end

end
