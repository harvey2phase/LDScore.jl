function get_compression(fh)
    if endswith(fh, "gz") return "gzip"
    elseif endswith(fh, "bz2") return "bz2"
    else return nothing end
end

function sumstats(fh; alleles=false, dropna=true)
    # Parses .sumstats files. See docs/file_formats_sumstats.txt
    #dtype_dict = {"SNP": str,   "Z": float, "N": float, "A1": str, "A2": str}
    compression = get_compression(fh)
    usecols = ["SNP", "Z", "N"]
    if alleles usecols += ["A1", "A2"] end

    try
        x = read_csv(fh, delim_whitespace=True, na_values=".")
    catch e
        println("Improperly formatted sumstats file:", str(e.args))
    end

    if dropna
        x = x.dropna(how="any")
    end

    return x
end

function _read_sumstats(args, fh; alleles=false, dropna=false)
    # Parse summary statistics
    @info "Reading summary statistics from" fh
    sumstats = sumstats(fh, alleles=alleles, dropna=dropna)
    @info "Read summary statistics for" size(sumstats) "SNPs."
    m = len(sumstats)
    sumstats = sumstats.drop_duplicates(subset="SNP")
    if m > size(sumstats)
        log.log(
            "Dropped {M} SNPs with duplicated rs numbers.".format(M=m - len(sumstats)))
    end
    return sumstats
end

function _read_ref_ld(args)
    # Read reference LD Scores
    ref_ld = _read_chr_split_files(args.ref_ld_chr, args.ref_ld,
                                   "reference panel LD Score", ps.ldscore_fromlist)
    log.log(
        "Read reference panel LD Scores for {N} SNPs.".format(N=len(ref_ld)))
    return ref_ld
end

#=
function _read_M(args, n_annot)
    # Read M (--M, --M-file, etc)
    if args.M:
        try:
            M_annot = [float(x) for x in _splitp(args.M)]
        except ValueError as e:
            raise ValueError("Could not cast --M to float: " + str(e.args))
    else:
        if args.ref_ld:
            M_annot = ps.M_fromlist(
                _splitp(args.ref_ld), common=(not args.not_M_5_50))
        elif args.ref_ld_chr:
            M_annot = ps.M_fromlist(
                _splitp(args.ref_ld_chr), _N_CHR, common=(not args.not_M_5_50))

    try:
        M_annot = np.array(M_annot).reshape((1, n_annot))
    except ValueError as e:
        raise ValueError(
            "# terms in --M must match # of LD Scores in --ref-ld.\n" + str(e.args))

    return M_annot
end

function _check_variance(M_annot, ref_ld)
    # Remove zero-variance LD Scores
    ii = ref_ld.ix[:, 1:].var() == 0  # NB there is a SNP column here
    if ii.all()
        raise ValueError("All LD Scores have zero variance.")
    else:
        log.log("Removing partitioned LD Scores with zero variance.")
        ii_snp = np.array([True] + list(~ii))
        ii_m = np.array(~ii)
        ref_ld = ref_ld.ix[:, ii_snp]
        M_annot = M_annot[:, ii_m]

    return M_annot, ref_ld, i
end

function _read_w_ld(args)
    # Read regression SNP LD
    if (args.w_ld and "," in args.w_ld) or (args.w_ld_chr and "," in args.w_ld_chr)
        raise ValueError(
            "--w-ld must point to a single fileset (no commas allowed).")
    w_ld = _read_chr_split_files(args.w_ld_chr, args.w_ld,
                                 "regression weight LD Score", ps.ldscore_fromlist)
    if len(w_ld.columns) != 2:
        raise ValueError("--w-ld may only have one LD Score column.")
    w_ld.columns = ["SNP", "LD_weights"]  # prevent colname conflicts w/ ref ld
    log.log(
        "Read regression weight LD Scores for {N} SNPs.".format(N=len(w_ld)))
    return w_ld
end

function _merge_and_log(ld, sumstats, noun)
    # Wrap smart merge with log messages about # of SNPs
    sumstats = smart_merge(ld, sumstats)
    msg = "After merging with {F}, {N} SNPs remain."
    if len(sumstats) == 0:
        raise ValueError(msg.format(N=len(sumstats), F=noun))
    else:
        log.log(msg.format(N=len(sumstats), F=noun))

    return sumstats
end
=#


# TODO
function _read_ld_sumstats(args, fh; alleles=false)
    sumstats = _read_sumstats(args, fh; alleles=alleles, dropna=dropna)
    ref_ld = _read_ref_ld(args)
    n_annot = size(ref_ld.columns) - 1
    M_annot = _read_M(args, n_annot)
    M_annot, ref_ld, novar_cols = _check_variance(M_annot, ref_ld)
    w_ld = _read_w_ld(args)
    sumstats = _merge_and_log(ref_ld, sumstats, "reference panel LD")
    sumstats = _merge_and_log(sumstats, w_ld, "regression SNP LD")
    w_ld_cname = sumstats.columns[-1]
    ref_ld_cnames = ref_ld.columns[1:len(ref_ld.columns)]

    return (M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols)
end

# TODO
function estimate_h2(args::Dict)
    args = Dict([
        ("h2", "height.sumstats.gz"),
        ("ref-ld-chr", "eur_w_ld_chr/"),
        ("w-ld-chr", "eur_w_ld_chr/"),
        ("out", "height_h2"),
    ])
    # Estimate h² and partitioned h².

    #if haskey(args, "samp_prev") && haskey(args, "pop_prev")
    #    (args["samp_prev"], args["pop_prev"]) = map(
    #        float, [args.samp_prev, args.pop_prev],
    #    )
    #end
    #if !args.intercept_h2 is not None:
    #    args.intercept_h2 = float(args.intercept_h2)
    #if args.no_intercept:
    #    args.intercept_h2 = 1

    M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols = _read_ld_sumstats(
        args, args["h2"])
    ref_ld = np.array(sumstats[ref_ld_cnames])
    _check_ld_condnum(args, ref_ld_cnames)
    _warn_length(sumstats)
    n_snp = len(sumstats)
    n_blocks = min(n_snp, args.n_blocks)
    n_annot = len(ref_ld_cnames)
    χ²_max = args.χ²_max
    old_weights = false
    #=
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
        overlap_matrix, M_tot = _read_annot(args)

        # overlap_matrix = overlap_matrix[np.array(~novar_cols), np.array(~novar_cols)]#np.logical_not
        df_results = ĥ²._overlap_output(ref_ld_cnames, overlap_matrix, M_annot, M_tot, args.print_coefficients)
        df_results.to_csv(args.out+".results", sep="\t", index=False)

    return ĥ²
    =#
end
