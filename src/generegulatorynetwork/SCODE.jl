using LinearAlgebra


using DelimitedFiles
fdata = "C:\\Users\\jcai.AUTH\\Downloads\\SCODE-master\\SCODE-master\\data2\\data.txt" # ARGS[1]
ftime = "C:\\Users\\jcai.AUTH\\Downloads\\SCODE-master\\SCODE-master\\data2\\time.txt" # ARGS[2]
X=readdlm(fdata,'\t')
pseudotime=readdlm(ftime,'\t')[:,2]
pseudotime = pseudotime/maximum(pseudotime)


cnum=length(pseudotime)
tfnum = 100   #parse(Int32,ARGS[4])
pnum = 4      #parse(Int32,ARGS[5])
maxite = 10

maxB = 2.0
minB = -10.0

W = zeros(tfnum, pnum)
Z = zeros(pnum, cnum)
WZ = zeros(tfnum, cnum)

new_B = zeros(pnum)
old_B = zeros(pnum)

#initialization
RSS = Inf
for i in 1:pnum
    new_B[i] = rand()*(maxB-minB)+minB
    old_B[i] = new_B[i]
end

function sample_Z()
    for i in 1:pnum
        for j in 1:cnum
            Z[i,j] = exp(new_B[i]*pseudotime[j]) + (rand()*0.002-0.001)
        end
    end
end

#optimization
for ite in 1:maxite
    #sampling B
    target = rand(1:pnum)
    global new_B[target] = rand()*(maxB-minB)+minB

    #for last calc
    if ite==maxite
        global new_B = copy(old_B)
    end

    #sample Z from new B
    sample_Z()

    #regression
    for i in 1:tfnum
        W[i,:] = (Z * Z')^-1 * (Z * X[i,:])
    end
    #prediction
    WZ = W * Z

    #RSS
    tmp_RSS = 0
    for i in 1:tfnum
        tmp_RSS += sum((X[i,:]-WZ[i,:]).*(X[i,:]-WZ[i,:]))
    end

    if tmp_RSS < RSS
        global RSS = tmp_RSS
        global old_B[target] = new_B[target]
    else
        global new_B[target] = old_B[target]
    end
end

#open("RSS.txt", "w") do f
#    println(f, RSS)
#end

#infer A
B = zeros(pnum, pnum)
for i in 1:pnum
    B[i,i] = new_B[i]
end
invW = pinv(W)
A = W * B * invW

#writedlm("A.txt",A)
#writedlm("B.txt",B)

