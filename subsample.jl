using BioSequences
using FASTX
using GZip
using CodecZlib

rawfastq_path = "/lovelace/echo/sequencing/mgx/rawfastq/"
outpath = "output/"
isdir(outpath) || mkdir(outpath)

samples = readlines("samples.txt")
fastqs = readdir(rawfastq_path)

function subsample(records, size)
    l = length(records)
    if size < l
        idx = rand(1:l, size)
        return view(records, idx)
    else
        return records
    end
end

function write_subsample(path, ss)
    io = FASTQ.Writer(GzipCompressorStream(open(path, "w")))
    for s in ss
        write(io, s)
    end
    close(io)
end

for sample in samples
    records = FASTQ.Record[]
    for fastq in filter(f-> occursin(replace(sample, "_"=>"-"), f) || occursin(sample, f), fastqs)
        file = joinpath(rawfastq_path, fastq)
        append!(records, FASTQ.Reader(GzipDecompressorStream(open(file))))
    end
    
    for i in 1:4
        write_subsample(joinpath(outpath, "$(sample)_10k_$i.fastq.gz"), subsample(records, 10_000))
        write_subsample(joinpath(outpath, "$(sample)_100k_$i.fastq.gz"), subsample(records, 100_000))
        write_subsample(joinpath(outpath, "$(sample)_1000k_$i.fastq.gz"), subsample(records, 1_000_000))
    end
end
