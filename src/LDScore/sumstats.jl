function estimate_h2(args::Dict)
    # Estimate h² and partitioned h².

    #=
    if haskey(args, "samp_prev") && haskey(args, "pop_prev")
        args["samp_prev"], args["pop_prev"] = map(
            float, [args.samp_prev, args.pop_prev])
    if args.intercept_h2 is not None:
        args.intercept_h2 = float(args.intercept_h2)
    if args.no_intercept:
        args.intercept_h2 = 1
    =#
    M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols = _read_ld_sumstats(
        args, log, args.h2)
    ref_ld = np.array(sumstats[ref_ld_cnames])
    _check_ld_condnum(args, log, ref_ld_cnames)
    _warn_length(log, sumstats)
    n_snp = len(sumstats)
    n_blocks = min(n_snp, args.n_blocks)
    n_annot = len(ref_ld_cnames)
    χ²_max = args.χ²_max
    old_weights = False
    if n_annot == 1:
        if args.two_step is None and args.intercept_h2 is None:
            args.two_step = 30
    else:
        old_weights = True
        if args.χ²_max is None:
            χ²_max = max(0.001*sumstats.N.max(), 80)

    s = lambda x: np.array(x).reshape((n_snp, 1))
    χ² = s(sumstats.Z**2)
    if χ²_max is not None:
        ii = np.ravel(χ² < χ²_max)
        sumstats = sumstats.ix[ii, :]

        n_snp = np.sum(ii)  # lambdas are late-binding, so this works
        ref_ld = np.array(sumstats[ref_ld_cnames])
        χ² = χ²[ii].reshape((n_snp, 1))

    ĥ² = reg.Hsq(χ², ref_ld, s(sumstats[w_ld_cname]), s(sumstats.N),
                     M_annot, n_blocks=n_blocks, intercept=args.intercept_h2,
                     twostep=args.two_step, old_weights=old_weights)


    if args.overlap_annot:
        overlap_matrix, M_tot = _read_annot(args, log)

        # overlap_matrix = overlap_matrix[np.array(~novar_cols), np.array(~novar_cols)]#np.logical_not
        df_results = ĥ²._overlap_output(ref_ld_cnames, overlap_matrix, M_annot, M_tot, args.print_coefficients)
        df_results.to_csv(args.out+'.results', sep="\t", index=False)

    return ĥ²
end
