using ArgParse

parser = ArgParseSettings()
@add_arg_table! parser begin
    "--two-step"
        arg_type = Float64
        help = """Test statistic bound for use with the two-step estimator.
            Not compatible with --no-intercept and --constrain-intercept."""
    "--h2"
        arg_type = String
        help = """Filename for a .sumstats[.gz] file for one-phenotype LD Score regression.
            --h2 requires at minimum also setting the --ref-ld and --w-ld flags."""
    "--h2-cts"
        arg_type = String
        help = """Filename for a .sumstats[.gz] file for cell-type-specific analysis.
            --h2-cts requires the --ref-ld-chr, --w-ld, and --ref-ld-chr-cts flags."""
    "--rg"
        arg_type = String
        help = "Comma-separated list of prefixes of .chisq filed for genetic correlation estimation."
    "--ref-ld"
        arg_type = String
        help = """Use --ref-ld to tell LDSC which LD Scores to use as the predictors in the LD Score regression.
            LDSC will automatically append .l2.ldscore/.l2.ldscore.gz to the filename prefix."""
    "--n-blocks"
        default = 200
        arg_type = Int64
        help = "Number of block jackknife blocks."
    "--χ²-max"
        arg_type = Float64
        help = "Max χ²."
    "--intercept-h²"
        help = "Intercepts for constrained-intercept single-trait LD Score regression."
    "--overlap-annot"
        default = false
        action = "store_true"
        help = """This flag informs LDSC that the partitioned LD Scores were generates using an
            annot matrix with overlapping categories (i.e., not all row sums equal 1),
            and prevents LDSC from displaying output that is meaningless with overlapping categories."""
end # add_arg_table parser
args = parse_args(ARGS, parser)

#= sanity check
print(args["two-step"])
=#

#= TODO arguments to be implemented
"--out", default = "ldsc", arg_type = String,
    help = "Output filename prefix. If --out is not set, LDSC will use ldsc as the "
    "defualt output filename prefix.")
# Basic LD Score Estimation Flags"
"--bfile", arg_type = String,
    help = "Prefix for Plink .bed/.bim/.fam file")
"--l2", default = False, action = "store_true",
    help = "Estimate l2. Compatible with both jackknife and non-jackknife.")
# Filtering / Data Management for LD Score
"--extract", arg_type = String,
    help = "File with SNPs to include in LD Score estimation. "
    "The file should contain one SNP ID per row.")
"--keep", arg_type = String,
    help = "File with individuals to include in LD Score estimation. "
    "The file should contain one individual ID per row.")
"--ld-wind-snps", arg_type = Int32,
    help = "Specify the window size to be used for estimating LD Scores in units of "
    "# of SNPs. You can only specify one --ld-wind-* option.")
"--ld-wind-kb", arg_type = Float64,
    help = "Specify the window size to be used for estimating LD Scores in units of "
    "kilobase-pairs (kb). You can only specify one --ld-wind-* option.")
"--ld-wind-cm", arg_type = Float64,
    help = "Specify the window size to be used for estimating LD Scores in units of "
    "centiMorgans (cM). You can only specify one --ld-wind-* option.")
"--print-snps", arg_type = String,
    help = "This flag tells LDSC to only print LD Scores for the SNPs listed "
    "(one ID per row) in PRINT_SNPS. The sum r^2 will still include SNPs not in "
    "PRINT_SNPs. This is useful for reducing the number of LD Scores that have to be "
    "read into memory when estimating h2 or rg." )
# Fancy LD Score Estimation Flags
"--annot", arg_type = String,
    help = "Filename prefix for annotation file for partitioned LD Score estimation. "
    "LDSC will automatically append .annot or .annot.gz to the filename prefix. "
    "See docs/file_formats_ld for a definition of the .annot format.")
"--thin-annot", action = "store_true", default = False,
    help = "This flag says your annot files have only annotations, with no SNP, CM, CHR, BP columns.")
"--cts-bin", arg_type = String,
    help = "This flag tells LDSC to compute partitioned LD Scores, where the partition "
    "is defined by cutting one or several continuous variable[s] into bins. "
    "The argument to this flag should be the name of a single file or a comma-separated "
    "list of files. The file format is two columns, with SNP IDs in the first column "
    "and the continuous variable in the second column. ")
"--cts-breaks", arg_type = String,
    help = "Use this flag to specify names for the continuous variables cut into bins "
    "with --cts-bin. For each continuous variable, specify breaks as a comma-separated "
    "list of breakpoints, and separate the breakpoints for each variable with an x. "
    "For example, if binning on MAF and distance to gene (in kb), "
    "you might set --cts-breaks 0.1,0.25,0.4x10,100,1000 ")
"--cts-names", arg_type = String,
    help = "Use this flag to specify names for the continuous variables cut into bins "
    "with --cts-bin. The argument to this flag should be a comma-separated list of "
    "names. For example, if binning on DAF and distance to gene, you might set "
    "--cts-bin DAF,DIST_TO_GENE ")
"--per-allele", default = False, action = "store_true",
    help = "Setting this flag causes LDSC to compute per-allele LD Scores, "
    "i.e., \ell_j : =  \sum_k p_k(1-p_k)r^2_{jk}, where p_k denotes the MAF "
    "of SNP j. ")
"--pq-exp", arg_type = Float64,
    help = "Setting this flag causes LDSC to compute LD Scores with the given scale factor, "
    "i.e., \ell_j : =  \sum_k (p_k(1-p_k))^a r^2_{jk}, where p_k denotes the MAF "
    "of SNP j and a is the argument to --pq-exp. ")
"--no-print-annot", default = False, action = "store_true",
    help = "By defualt, seting --cts-bin or --cts-bin-add causes LDSC to print "
    "the resulting annot matrix. Setting --no-print-annot tells LDSC not "
    "to print the annot matrix. ")
"--maf", arg_type = Float64,
    help = "Minor allele frequency lower bound. Default is MAF > 0.")
# Basic Flags for Working with Variance Components
"--ref-ld-chr", arg_type = String,
    help = "Same as --ref-ld, but will automatically concatenate .l2.ldscore files split "
    "across 22 chromosomes. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz "
    "to the filename prefix. If the filename prefix contains the symbol @, LDSC will "
    "replace the @ symbol with chromosome numbers. Otherwise, LDSC will append chromosome "
    "numbers to the end of the filename prefix."
    "Example 1: --ref-ld-chr ld/ will read ld/1.l2.ldscore.gz ... ld/22.l2.ldscore.gz"
    "Example 2: --ref-ld-chr ld/@_kg will read ld/1_kg.l2.ldscore.gz ... ld/22_kg.l2.ldscore.gz")
"--w-ld", arg_type = String,
    help = "Filename prefix for file with LD Scores with sum r^2 taken over SNPs included "
    "in the regression. LDSC will automatically append .l2.ldscore/.l2.ldscore.gz.")
"--w-ld-chr", arg_type = String,
    help = "Same as --w-ld, but will read files split into 22 chromosomes in the same "
    "manner as --ref-ld-chr.")
"--print-coefficients",default = False,action = "store_true",
    help = "when categories are overlapping, print coefficients as well as heritabilities.")
"--frqfile", arg_type = String,
    help = "For use with --overlap-annot. Provides allele frequencies to prune to common "
    "snps if --not-M-5-50 is not set.")
"--frqfile-chr", arg_type = String,
    help = "Prefix for --frqfile files split over chromosome.")
"--no-intercept", action = "store_true",
    help = "If used with --h2, this constrains the LD Score regression intercept to equal "
    "1. If used with --rg, this constrains the LD Score regression intercepts for the h2 "
    "estimates to be one and the intercept for the genetic covariance estimate to be zero.")
"--intercept-gencov", action = "store", default = None,
    help = "Intercepts for constrained-intercept cross-trait LD Score regression."
    " Must have same length as --rg. The first entry is ignored.")
"--M", arg_type = String,
    help = "# of SNPs (if you don\"t want to use the .l2.M files that came with your .l2.ldscore.gz files)")
"--ref-ld-chr-cts", arg_type = String,
    help = "Name of a file that has a list of file name prefixes for cell-type-specific analysis.")
"--print-all-cts", action = "store_true", default = False)

# Flags for both LD Score estimation and h2/gencor estimation
"--print-cov", default = False, action = "store_true",
    help = "For use with --h2/--rg. This flag tells LDSC to print the "
    "covaraince matrix of the estimates.")
"--print-delete-vals", default = False, action = "store_true",
    help = "If this flag is set, LDSC will print the block jackknife delete-values ("
    "i.e., the regression coefficeints estimated from the data with a block removed). "
    "The delete-values are formatted as a matrix with (# of jackknife blocks) rows and "
    "(# of LD Scores) columns.")
# Flags you should almost never use
"--chunk-size", default = 50, arg_type = Int32,
    help = "Chunk size for LD Score calculation. Use the default.")
"--pickle", default = False, action = "store_true",
    help = "Store .l2.ldscore files as pickles instead of gzipped tab-delimited text.")
"--yes-really", default = False, action = "store_true",
    help = "Yes, I really want to compute whole-chromosome LD Score.")
"--invert-anyway", default = False, action = "store_true",
    help = "Force LDSC to attempt to invert ill-conditioned matrices.")
"--not-M-5-50", default = False, action = "store_true",
    help = "This flag tells LDSC to use the .l2.M file instead of the .l2.M_5_50 file.")
"--return-silly-things", default = False, action = "store_true",
    help = "Force ldsc to return silly genetic correlation estimates.")
"--no-check-alleles", default = False, action = "store_true",
    help = "For rg estimation, skip checking whether the alleles match. This check is "
    "redundant for pairs of chisq files generated using munge_sumstats.py and the "
    "same argument to the --merge-alleles flag.")
# transform to liability scale
"--samp-prev",default = None,
    help = "Sample prevalence of binary phenotype (for conversion to liability scale).")
"--pop-prev",default = None,
    help = "Population prevalence of binary phenotype (for conversion to liability scale).")

if __name__  =  =  "__main__":

    args = parser.parse_args()
    if args.out is None:
        raise ValueError("--out is required.")

    log = Logger(args.out+".log")
    try:
        defaults = vars(parser.parse_args(""))
        opts = vars(args)
        non_defaults = [x for x in opts.keys() if opts[x] ! =  defaults[x]]
        header = MASTHEAD
        header + =  "Call: \n"
        header + =  "./ldsc.py \\\n"
        options = ["--"+x.replace("_","-")+" "+String(opts[x])+" \\" for x in non_defaults]
        header + =  "\n".join(options).replace("True","").replace("False","")
        header = header[0:-1]+"\n"
        log.log(header)
        log.log("Beginning analysis at {T}".format(T = time.ctime()))
        start_time = time.time()
        if args.n_blocks < =  1:
            raise ValueError("--n-blocks must be an integer > 1.")
        if args.bfile is not None:
            if args.l2 is None:
                raise ValueError("Must specify --l2 with --bfile.")
            if args.annot is not None and args.extract is not None:
                raise ValueError("--annot and --extract are currently incompatible.")
            if args.cts_bin is not None and args.extract is not None:
                raise ValueError("--cts-bin and --extract are currently incompatible.")
            if args.annot is not None and args.cts_bin is not None:
                raise ValueError("--annot and --cts-bin are currently incompatible.")
            if (args.cts_bin is not None) ! =  (args.cts_breaks is not None):
                raise ValueError("Must set both or neither of --cts-bin and --cts-breaks.")
            if args.per_allele and args.pq_exp is not None:
                raise ValueError("Cannot set both --per-allele and --pq-exp (--per-allele is equivalent to --pq-exp 1).")
            if args.per_allele:
                args.pq_exp = 1


            ldscore(args, log)
        # summary statistics
        elif (args.h2 or args.rg or args.h2_cts) and (args.ref_ld or args.ref_ld_chr) and (args.w_ld or args.w_ld_chr):
            if args.h2 is not None and args.rg is not None:
                raise ValueError("Cannot set both --h2 and --rg.")
            if args.ref_ld and args.ref_ld_chr:
                raise ValueError("Cannot set both --ref-ld and --ref-ld-chr.")
            if args.w_ld and args.w_ld_chr:
                raise ValueError("Cannot set both --w-ld and --w-ld-chr.")
            if (args.samp_prev is not None) ! =  (args.pop_prev is not None):
                raise ValueError("Must set both or neither of --samp-prev and --pop-prev.")

            if not args.overlap_annot or args.not_M_5_50:
                if args.frqfile is not None or args.frqfile_chr is not None:
                    log.log("The frequency file is unnecessary and is being ignored.")
                    args.frqfile = None
                    args.frqfile_chr = None
            if args.overlap_annot and not args.not_M_5_50:
                if not ((args.frqfile and args.ref_ld) or (args.frqfile_chr and args.ref_ld_chr)):
                    raise ValueError("Must set either --frqfile and --ref-ld or --frqfile-chr and --ref-ld-chr")

            if args.rg:
                sumstats.estimate_rg(args, log)
            elif args.h2:
                sumstats.estimate_h2(args, log)
            elif args.h2_cts:
                sumstats.cell_type_specific(args, log)

            # bad flags
        else:
            print header
            print "Error: no analysis selected."
            print "ldsc.py -h describes options."
    except Exception:
        ex_type, ex, tb = sys.exc_info()
        log.log( traceback.format_exc(ex) )
        raise
    finally:
        log.log("Analysis finished at {T}".format(T = time.ctime()) )
        time_elapsed = round(time.time()-start_time,2)
        log.log("Total time elapsed: {T}".format(T = sec_to_str(time_elapsed)))
=#
