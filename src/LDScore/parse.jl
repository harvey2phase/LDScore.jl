using Pandas

function series_eq(x, y):
function read_csv(fh, **kwargs):
function sub_chr(s, chrom):
function get_present_chrs(fh, num):
function which_compression(fh):
function read_cts(fh, match_snps):


function get_compression(fh)
    if endswith(fh, "gz") return "gzip"
    elseif endswith(fh, "bz2") return "bz2"
    else return nothing end
end


function parse_sumstats(fh; alleles=false, dropna=true)
    dtype_dict = Dict([
        ("SNP", String),
        #("Z", Float),
      #  ("N", float),
      #  ("A1", str),
      #  ("A2", str),
    ])

    compression = get_compression(fh)
    usecols = ["SNP", "Z", "N"]
    if alleles usecols += ["A1", "A2"] end
    @info "fh" fh

    @info "NEW"
    x = Pandas.read_csv(fh, usecols=usecols, dtype=dtype_dict, compression=compression)
    #=
    try
        x = Pandas.read_csv(fh)
    catch e
        @error "Improperly formatted sumstats file."
    end
    =#

    if dropna
        x = x.dropna(how="any")
    end
    @info x

    return x
end


function ldscore_fromlist(flist, num=None):
function l2_parser(fh, compression):
function annot_parser(fh, compression, frqfile_full=None, compression_frq=None):
function frq_parser(fh, compression):
function ldscore(fh, num=None):
function M(fh, num=None, N=2, common=False):
function M_fromlist(flist, num=None, N=2, common=False):
function annot(fh_list, num=None, frqfile=None):
function __ID_List_Factory__(colnames, keepcol, fname_end, header=None, usecols=None):
    struct IDContainer end
        function __init__(self, fname):
        function __read__(self, fname):
        function loj(self, externalDf):
