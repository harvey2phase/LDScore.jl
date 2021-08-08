function test_print(name, var)
    println("=================================================================")
    println(name)
    println(var)
    println("=================================================================")
end

using CSV
using DataFrames


#using Pkg; Pkg.add("Revise")
#using Revise

include("parse.jl")

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
    # Parser summary statistics
    @info "Reading summary statistics from" fh
    sumstats = parse_sumstats(fh, alleles=alleles, dropna=dropna)
    @info "Read summary statistics for" size(sumstats) "SNPs."
    m = size(sumstats)
    sumstats = sumstats.drop_duplicates(subset="SNP")
    if m > size(sumstats)
        @info "Dropped SNPs with duplicated rs numbers:" m-size(sumstats)
    end
    return sumstats
end

function _read_ref_ld(args)
    # Read reference LD Scores
    ref_ld = _read_chr_split_files(args["ref_ld_chr"], args["ref_ld"],
                                   "reference panel LD Score", ldscore_fromlist)
    @info "Read reference panel LD Scores for {N} SNPs." len(ref_ld)
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

function _load_testset_1(args)
    df_ref_id = DataFrames.DataFrame(CSV.File(args["ref_id"] * ".l2.ldscore"))
    df_h2 = DataFrames.DataFrame(CSV.File(args["h2"]))

    M_annot = [[155881.2526]]
    w_ld_cname = "CHR"
    ref_ld_cnames = ["LD_0"]
    sumstats = DataFrames.innerjoin(df_ref_id, df_h2, on = "SNP")
    novar_cols = Dict([("LD_0", false)])
    return (M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols)
end

# TODO
function _read_ld_sumstats(args; alleles=false, dropna=true)
    return _load_testset_1(args)
end

# TODO
function estimate_h2(args::Dict)
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

    (M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols) = _read_ld_sumstats(
        args)
    ref_ld = (sumstats[!, "LD"])
    n_snp = size(sumstats)[1]
    test_print("n_snp", n_snp)
    n_blocks = min(n_snp, args["n_blocks"])
    n_annot = size(ref_ld_cnames)[1]
    χ²_max = args["χ²_max"]
    old_weights = false

    if n_annot == 1
        if args["two_step"] == nothing && args["intercept_h2"] == nothing
            args["two_step"] = 30
        end
    else
        old_weights = true
        if args["χ²_max"] == nothing
            χ²_max = max(0.001 * max(sumstats[!, "N"]), 80)
        end
    end

    s(x) = reshape(x, (n_snp, 1))
    χ² = s(sumstats[!, "Z"] .^ 2)

    if χ²_max != nothing
        ii = vec(χ² < χ²_max)
        sumstats = sumstats.ix[ii, :]

        n_snp = sum(ii) # TODO double check that this works
        ref_ld = [sumstats[ref_ld_cnames]]
        χ² = reshape(χ²[ii], (n_snp, 1))
    end

    ĥ² = Hsq(
        χ², ref_ld, s(sumstats[!, w_ld_cname]), s(sumstats[!, "N"]), M_annot,
        n_blocks, args["intercept_h2"], false, args["two_step"], old_weights,
    )

    if args["overlap_annot"]
        overlap_matrix, M_tot = _read_annot(args)

        # overlap_matrix = overlap_matrix[np.array(~novar_cols), np.array(~novar_cols)]#np.logical_not
        df_results = ĥ²._overlap_output(ref_ld_cnames, overlap_matrix, M_annot, M_tot, args.print_coefficients)
        df_results.to_csv(args.out+".results", sep="\t", index=False)
    end

    return ĥ²
end
