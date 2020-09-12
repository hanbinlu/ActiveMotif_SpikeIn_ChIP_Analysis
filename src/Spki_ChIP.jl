#' # Spike-in ChIP Analysis
#' ## Modify the following parameters to control the pipeline
#' ### Program Control Parameters
ncores, picard = 35,  "/home/software/picard.jar";
mapq = 30; #currently only filter unique mapped reads by mapq cutoff
cleanup, mktrack = true, true
#' ### Experiment Parameters
refgenome = "/home/software/bowtie2-2.2.9/genome/mm9/mm9";
spkgenome = "/home/software/bowtie2-2.2.9/genome/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome";
fastq_folder = "/Extension_HDD2/Hanbin/ProjectW/200511_HL02_ChIP/HL02_raw_Fastq/";
out_dir = ""
fastq = Dict{String,String}(
    "M2_3C10_0hr" => "HL02_Idx17_S102_L004_R1_001.fastq.gz",
    "M2_3C10_16hr" => "HL02_Idx24_S103_L004_R1_001.fastq.gz",
)
n, samples = length(fastq), ["M2_3C10_0hr", "M2_3C10_16hr"]
samples = [joinpath(out_dir, s) for s in samples]
#+
#' ## Mapping to reference genome
#+
function extract_bowtie_result(mapstat)
    lines = split(mapstat, '\n')
    tot = parse(Int, split(lines[1])[1])
    uniq = parse(Int, split(lines[4])[1])
    tot, uniq
end
#data to capture
total, bowtie_unique, rmdup_unique, norm = zeros(Int, n), [zeros(Int, 2) for i=1:n], 
    [zeros(Int, 2) for i=1:n], zeros(Int, n);
#+
for (i, sample) in enumerate(samples)
    fq = fastq[sample]
    input = joinpath(fastq_folder, fq)
    sam = "$sample.ref.sam"
    err = Pipe()
    println("Mapping $sample to reference genome")
    run(pipeline(`bowtie2 -p $ncores -x $refgenome $input`, stdout=sam, stderr=err))
    close(err.in)
    stats = read(err.out)
    tot, uniq = extract_bowtie_result(String(stats))
    total[i] = tot
    bowtie_unique[i][1] = uniq
    print(String(stats))
end
#' ## Mapping to spikin genome
#+
for (i, sample) in enumerate(samples)
    fq = fastq[sample]
    input = joinpath(fastq_folder, fq)
    sam = "$sample.spki.sam"
    err = Pipe()
    println("Mapping $sample to spike-in genome")
    run(pipeline(`bowtie2 -p $ncores -x $spkgenome $input`, stdout=sam, stderr=err))
    close(err.in)
    stats = String(read(err.out))
    tot, uniq = extract_bowtie_result(stats)
    total[i] = tot
    bowtie_unique[i][2] = uniq
    println(stats)
end
#' ## Sort and remove duplication
#+
function rmdup(sample, flag)
    genome = flag == :ref ? "ref" : "spki"
    input, sorted_bam, rmdup_bam = "$sample.$genome.sam", "$sample.$genome.sorted.bam", "$sample.$genome.rmdup.bam"
    cmd_sort = Cmd(`samtools sort -@ $ncores -o $sorted_bam $input`)
    cmd_rmdup = Cmd(`java -jar $picard MarkDuplicates REMOVE_DUPLICATES=true I=$sorted_bam O=$rmdup_bam M=/dev/null`)
    run(cmd_sort)
    run(cmd_rmdup)
    nothing
end
#+
for sample in samples
    rmdup(sample, :ref)
    rmdup(sample, :spki)
end  
#+
#' ## Normalization factors
using XAM
function count_bam_uniq(sample, flag)
    genome = flag == :ref ? "ref" : "spki"
    bam = "$sample.$genome.rmdup.bam"
    count = 0
    reader = open(BAM.Reader, bam)
    for record in reader
        if BAM.ismapped(record) && Int(BAM.mappingquality(record)) >= mapq
            count += 1
        end
    end
    close(reader)
    count
end
for (i,sample) in enumerate(samples)
    counts = [count_bam_uniq(sample, :ref), count_bam_uniq(sample, :spki)]
    rmdup_unique[i] = counts
end
#' ## Report
using DataFrames, PrettyTables
norm_to = minimum([v[2] for v in rmdup_unique])
df = DataFrame([
    :Sample => samples, 
    :Sequence_Depth => total,
    :Unique_Ref => [v[1] for v in bowtie_unique],
    :Unique_Spki => [v[2] for v in bowtie_unique],
    :Rmdup_Unique_Ref => [v[1] for v in rmdup_unique],
    :Rmdup_Unique_Spki => [v[2] for v in rmdup_unique],
    :Pulldown_Ratio=>[Int(ceil(v[1]*norm_to/v[2])) for v in rmdup_unique],
])
pretty_table(df)
#' ## Create Normalized Track
using Random
using BGZFStreams
if mktrack
    for (i,sample) in enumerate(samples)
        ratio = norm_to/rmdup_unique[i][2]
        bam = "$sample.ref.rmdup.bam"
        reader = open(BAM.Reader, bam)
        header = BAM.header(reader)
        out_bam = "$sample.spki_normed.bam"
        out = BAM.Writer(BGZFStream(open(out_bam, "w"), "w"), header)
        for record in reader
            if BAM.ismapped(record) && Int(BAM.mappingquality(record)) >= mapq
                rand() <= ratio && write(out, record)
            end
        end
        close(reader)
        close(out)
        #create track
        out_sam = "$sample.spki_normed.sam"
        run(`samtools view -h -o  $out_sam $out_bam`)
        run(`makeTagDirectory $(sample).spki_normed/ $out_sam -keepAll`)
        run(`makeUCSCfile $(sample).spki_normed/ -o auto -raw`)
    end
end
#' ## Clean intermediate file
if cleanup
    for f in readdir()
        if endswith(f, "sam") || endswith(f, "bam")
            rm(f)
        end
    end
end