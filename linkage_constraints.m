%%-----  formulation  --------------------------------------------------
function om = linkage_constraints(om, cost_params)

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

ng= cost_params.ng;
ID_gen = cost_params.ID_gen;
v_IDs = cost_params.v_IDs;
nb = cost_params.nb;
nr_con = cost_params.nc;

npg = size(ID_gen,1);

L_PG = zeros(npg*nr_con,1);
U_PG = zeros(npg*nr_con,1);

A_PG = zeros(npg*nr_con,ng*(nr_con+1));

for g = 1:npg
    for c=1:nr_con
        A_PG(g+(c-1)*npg,ID_gen(g)) = 1;
        A_PG(g+(c-1)*npg,ID_gen(g)+ng*c) = -1;
    end
end

om = add_constraints(om, 'Linkage_Constraints_PG', A_PG, L_PG, U_PG, {'Pg'});

nvg =size(v_IDs,1);

L_VG = zeros(nvg*nr_con,1);
U_VG = zeros(nvg*nr_con,1);

A_VG = zeros(nvg*nr_con,nb*(nr_con+1));

for g = 1:nvg
    for c=1:nr_con
        A_VG(g+(c-1)*nvg,v_IDs(g)) = 1;
        A_VG(g+(c-1)*nvg,v_IDs(g)+nb*c) = -1;
    end
end

om = add_constraints(om, 'Linkage_Constraints_VG', A_VG, L_VG, U_VG, {'Vm'});

end