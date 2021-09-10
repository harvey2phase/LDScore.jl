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


function _read_ref_ld(args)
    # Read reference LD Scores
    ref_ld = _read_chr_split_files(
        args["ref_ld_chr"], args["ref_ld"],
        "reference panel LD Score", ldscore_fromlist,
    )
    @info "Num. of SNPs read for reference panel LD Scores:" size(ref_ld)
    return ref_ld
end


function _load_testset_1(args)
    df_ref_id = DataFrames.DataFrame(CSV.File(args["ref-ld"] * ".l2.ldscore"))
    df_h2 = DataFrames.DataFrame(CSV.File(args["h²"]))

    M_annot = reshape([[155881.2526]], (1, 1))
    w_ld_cname = "CHR"
    ref_ld_cnames = ["LD_0"]
    sumstats = DataFrames.innerjoin(df_ref_id, df_h2, on = "SNP")
    novar_cols = Dict([("LD_0", false)])
    return (M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols)
end


# TODO: Just loading testing set now, need to actually implement to load
#       generic datasets
function _read_ld_sumstats(args; alleles=false, dropna=true)
    return _load_testset_1(args)
end


# TODO: Needs more error-checking
function estimate_h2()
    (
        M_annot, w_ld_cname, ref_ld_cnames, sumstats, novar_cols,
    ) = _read_ld_sumstats(LDScoreJulia.args)

    ref_ld = sumstats[!, "LD"]
    ref_ld = reshape(ref_ld, (size(ref_ld)[1], 1))
    n_snp = size(sumstats)[1]
    n_blocks = min(n_snp, LDScoreJulia.args["n-blocks"])
    n_annot = size(ref_ld_cnames)[1]
    χ²_max = LDScoreJulia.args["χ²-max"]
    old_weights = false

    if n_annot == 1
        if LDScoreJulia.args["two-step"] == nothing && LDScoreJulia.args["intercept-h²"] == nothing
            LDScoreJulia.args["two-step"] = 30
        end
    else
        old_weights = true
        if LDScoreJulia.args["χ²-max"] == nothing
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
        n_blocks, LDScoreJulia.args["intercept-h²"], false, LDScoreJulia.args["two-step"], old_weights,
    )

    if LDScoreJulia.args["overlap_annot"]
        overlap_matrix, M_tot = _read_annot(LDScoreJulia.args)

        df_results = ĥ²._overlap_output(
            ref_ld_cnames, overlap_matrix, M_annot, M_tot,
            LDScoreJulia.args["print_coefficients"],
        )
        df_results.to_csv(LDScoreJulia.args["out"]+".results", sep="\t", index=False)
    end

    return ĥ²
end # estimate_h2
