using CSV
using DataFrames

include("parse.jl")

function test_print(name, value)
    println("=================================================================")
    println(name)
    println(value)
    println("-----------------------------------------------------------------")
end # test_print


# TODO
function _read_M(args, n_annot) end
function _check_variance(M_annot, ref_ld) end
function _read_w_ld(args) end
function _merge_and_log(ld, sumstats, noun) end


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


function _read_ref_ld()
    # Read reference LD Scores
    ref_ld = _read_chr_split_files(
        LDScoreJulia.args["ref-ld-chr"], LDScoreJulia.args["ref-ld"],
        "reference panel LD Score", ldscore_fromlist,
    )
    @info "Num. of SNPs read for reference panel LD Scores:" size(ref_ld)
    return ref_ld
end


function _load_testset_1()
    df_ref_id = DataFrames.DataFrame(CSV.File(LDScoreJulia.args["ref-ld"] * ".l2.ldscore"))
    df_h2 = DataFrames.DataFrame(CSV.File(LDScoreJulia.args["h²"]))

    M_annot = reshape([[155881.2526]], (1, 1))
    w_ld_cname = "CHR"
    ref_ld_cnames = ["LD_0"]
    sumstats = DataFrames.innerjoin(df_ref_id, df_h2, on = "SNP")
    novar_cols = Dict([("LD_0", false)])
    return (M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols)
end


#= TODO: Just loading testing set now,
    need to actually implement to load generic datasets =#
function _read_ld_sumstats(alleles=false, dropna=true)
    return _load_testset_1()
end


""" Estimate h² and partitioned h² """
function estimate_h2()
    args = deepcopy(LDScoreJulia.args)

    # TODO: more checking


    M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols = _read_ld_sumstats()

    ref_ld = sumstats[!, "LD"]
    ref_ld = reshape(ref_ld, (size(ref_ld)[1], 1))
    n_snp = size(sumstats)[1]
    n_blocks = min(n_snp, args["n-blocks"])
    n_annot = size(ref_ld_cnames)[1]
    χ²_max = args["χ²-max"]
    old_weights = false

    if n_annot == 1
        if args["two-step"] == nothing && args["intercept-h²"] == nothing
            args["two-step"] = 30
        end
    else
        old_weights = true
        if args["χ²-max"] == nothing
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
        n_blocks = n_blocks, intercept = args["intercept-h²"],
        slow = false, step1_ii = args["two-step"],
        old_weights = old_weights,
    )

    if args["overlap-annot"]
        overlap_matrix, M_tot = _read_annot(args)

        df_results = ĥ²._overlap_output(
            ref_ld_cnames, overlap_matrix, M_annot, M_tot,
            args["print-coefficients"],
        )
        df_results.to_csv(
            args["out"] * ".results",
            sep = "\t", index = false,
        )
    end

    return ĥ²
end # function estimate_h2
