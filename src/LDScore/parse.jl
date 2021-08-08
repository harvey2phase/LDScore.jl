using CSV


# TODO
#function series_eq(x, y) end
#function read_csv(fh, **kwargs) end
#function sub_chr(s, chrom) end
#function get_present_chrs(fh, num) end
#function which_compression(fh) end
#function read_cts(fh, match_snps) end
#function ldscore_fromlist(flist, num=nothing) end
#function l2_parser(fh, compression) end
#function annot_parser(
#    fh, compression, frqfile_full=nothing, compression_frq=nothing,
#) end
#function frq_parser(fh, compression) end
#function ldscore(fh, num=nothing) end
#function M(fh; num=nothing, N=2, common=False) end
#function M_fromlist(flist; num=nothing, N=2, common=False) end
#function annot(fh_list; num=nothing, frqfile=nothing) end
#function __ID_List_Factory__(
#    colnames, keepcol, fname_end; header=nothing, usecols=nothing,
#) end
#    struct IDContainer end
#        function __init__(container::IDContainer, fname) end
#        function __read__(container::IDContainer, fname) end
#        function loj(container::IDContainer, externalDf) end


function get_compression(fh)
    if endswith(fh, "gz") return "gzip"
    elseif endswith(fh, "bz2") return "bz2"
    else return nothing end
end


function parse_sumstats(fh; alleles=false, dropna=true)
    @info "Reading from:" fh
    x = CSV.File(fh)

    if dropna x = x.dropna(how="any") end # TODO double check Julian equivalent

    @info x
    return x
end
