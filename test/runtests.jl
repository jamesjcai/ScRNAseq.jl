# push!(LOAD_PATH,"E:/GitHub/ScRNAseq.jl/src/");
using ScRNAseq                          # obviously the tests use the ScTenifoldNet module... 
using Test                              # and the Base.Test module... 
tests = ["code_test1"]         # the test file names are stored as strings... 
for t in tests
  include("$(t).jl")                         # ... so that they can be evaluated in a loop 
end
