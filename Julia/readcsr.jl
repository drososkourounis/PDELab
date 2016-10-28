function readcsr(filename::String)
    
    print("- Opening file ", filename, " for reading... \n");
    f = open(filename);

    nrows = parse(Int, readline(f));
    ncols = parse(Int, readline(f));
    nnz   = parse(Int, readline(f));

    print("  CSR matrix with ", nrows, " rows, ", ncols, " columns, and ", nnz, " nonzeros. \n  ... \n");

    prows = Array{Int,1}(nrows+1);
    pcols = Array{Int,1}(nnz);
    pdata = Array{Float64,1}(nnz);

    for i in 1:nrows+1
        prows[i] = parse(Int, readline(f));
    end

    for i in 1:nnz
        pcols[i] = parse(Int, readline(f));
    end

    for i in 1:nnz
        pdata[i] = parse(Float64, readline(f));
    end


    #print(prows, "\n", pcols, "\n", pdata, "\n");

    if prows[end] == nnz
        prows = prows+1;
        pcols = pcols+1;
    end

    println("-  Done!");
    close(f);

    #return (prows, pcols, pdata);
    
    return (SparseMatrixCSC(nrows, ncols, prows, pcols, pdata))'

end
