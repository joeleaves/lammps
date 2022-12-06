function [rdf,err] = get_gr(g,length_gr)
%Compute the RDF and the standard error from a two-column file based on a
%LAMMPS output for g(r). That file is g and the number of lines in each
%snapshot of the RDF is length_gr.

    n = length(g);
    rdf = zeros(length_gr,1);
    rdf2 = zeros(length_gr,1);
    k = 1;
    start = 1;
    while(start+length_gr < n)
        start = (k-1)*length_gr + k;
        stop = start + length_gr;
        g_tmp = g(start:stop);
        g_tmp(1) = [];
        k = k + 1;
     [   rdf = rdf + g_tmp;
        rdf2 = rdf2 + g_tmp.^2;
    end
    rdf = rdf/k;
    rdf2 = rdf2/k;
    err = (rdf2 - rdf.^2).^(1/2);


end