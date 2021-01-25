module FileIO

using DelimitedFiles, MatrixMarket, MAT, CSV, DataFrames, HDF5
export readtxt,   # Text file   
	   readmtx,   # MatrixMarket file
	   readmat,   # MATLAB mat file
	   readhdf,   # HDF5 file
       readgenelist       


function readmtx(filename::AbstractString)
    if isfile(filename)
	# read MatrixMarket file
	X=MatrixMarket.mmread(filename);
    else
        error("File $(filename) does not exist.")
    end
    return X
end

function readcsv(filename::AbstractString)
    # read CSV file with header and row name
    df = CSV.File(filename; datarow=2) |> DataFrame!
    X=convert(Matrix, df[:,2:end])
    genelist=df[:,1]
    return X,genelist    
end

function readmat(filename::AbstractString)
    # read Matlab Mat file 
    file=matopen(filename)
    X=read(file,"X")
    genelist=read(file,"genelist")
    close(file)
    X = convert(Array{Float64,2}, X)
    return X,genelist
end

function readtxt(filename::AbstractString)
    # read DLM text file
    X=readdlm(filename,',',Int16)
end

function readhdf(filename::AbstractString)
    # https://anndata.readthedocs.io/en/latest/anndata.AnnData.html
end

function readgenelist(filename::AbstractString,colidx::Integer=1)
    # read genelist
    genelist=readdlm(filename,'\t',String)
    genelist=vec(genelist[:,colidx])
end

# using UnicodePlots
# function showx(X)
#    spy(X)
# end

end
