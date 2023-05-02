module FileIO

using DelimitedFiles, MatrixMarket, MAT, CSV, DataFrames, HDF5
export readtxt,   # Text file   
	   readmtx,   # MatrixMarket file
	   readmat,   # MATLAB mat file
	   readhdf,   # HDF5 file
       read10xh5, # filtered_feature_bc_matrix.h5
       readgenelist       

"""
    readmtx(filename::AbstractString) -> AbstractMatrix

Load and return the content of the MatrixMarket file.

# Arguments

- `filename::AbstractString`: the filename of the file to load,

# Keywords

# Returns
- `AbstractMatrix`: the UMI matrix of the file content,

# Throws
- `Error`: in the case of file does not exist.
"""

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
    g=df[:,1]
    return X,g    
end

function readmat(filename::AbstractString)
    # read Matlab Mat file X is sparse matrix, g is a cellstr with genenames
    file=matopen(filename)
    X=read(file,"X")
    g=read(file,"g")
    close(file)
    #X = convert(Array{Float64,2}, X)
    #X = convert(Matrix{Int16}, X)
    return X,g
end

function readtxt(filename::AbstractString)
    # read DLM text file
    X=readdlm(filename,',',Int16)
end

function readhdf(filename::AbstractString)
    # https://anndata.readthedocs.io/en/latest/anndata.AnnData.html
end

function read10xh5(filename::AbstractString)
    # filtered_feature_bc_matrix.h5
    # https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices
    f=h5open(filename,"r")
    d=f["matrix/data"]
    data=read(d)
    indices=read(f["matrix/indices"])
    indptr=read(f["matrix/indptr"])
    shape=read(f["matrix/shape"])    
    X=zeros(eltype(data),shape[1],shape[2])
    for k=1:(length(indptr)-1)
        idx=(indptr[k]+1):indptr[k+1]
        y=indices[idx].+1
        X[y,k].=data[idx]
    end
    g=read(f["matrix/features/name"])
    close(f)
    return X,g
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