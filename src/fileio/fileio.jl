module Fileio

using DelimitedFiles, MatrixMarket, UnicodePlots, MAT, CSV, DataFrames
export readmm,
       readtx,
       readgl,
       showx


function readmm(filename)
    # read MatrixMarket file
    X=mmread(filename);
end

function readcsv(filename)
    # read CSV file with header and row name
    df = CSV.File(filename; datarow=2) |> DataFrame!
    X=convert(Matrix, df[:,2:end])
    genelist=df[:,1]
    return X,genelist    
end

function readmt(filename)
    # read Matlab Mat file 
    file=matopen(filename)
    X=read(file,"X")
    genelist=read(file,"genelist")
    close(file)
    X = convert(Array{Float64,2}, X)
    return X,genelist
end

function readtx(filename)
    # read DLM text file
    X=readdlm(filename,',',Int16)
end

function readgl(filename::String,colidx::Integer=1)
    # read genelist
    genelist=readdlm(filename,'\t',String)
    genelist=vec(genelist[:,colidx])
end

function showx(X)
    spy(X)
end

end