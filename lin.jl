using Oscar;

# map systems to their newton polytope
# identify common edges
# greed: prefer edges that appear in more vs appear in less
# identify translations (iso to monomials) required to attain these
# greed, prefer monomials which appear more vs less
# identify circuits of associated matrix by columns, try to maximise


# form matrix out of systems
# form n-minors out of this matrix

S,(x,y)=laurent_polynomial_ring(QQ,["x","y"])
# coeff does not work for laurent_polynomial_ring, cheap workaround
R,(a,b)=polynomial_ring(QQ,["a","b"])

# note: rand(::QQMPolyRing
system=[x^2+y^2+2x+3y+7x^2*y^3,x^2+y^7+2x^2+3y^3+7x+2y+x^7*y^7,x^2+y^2+1]
system=[x^5*(1+x+y+x*y),x+y+x*y,x^3+y]
#system=[58*x^4*y^6 + 97*x^3*y^9, 12*x^8*y^9 + 92*x^5*y^9 + 29*x^5*y^2 + 99*x^4*y^5 + 92*x^4*y + 76*x^2*y^3 + 44*x*y^2, 27*x^10*y^5 + 87*x^9*y^8 + 70*x^8*y^5 + 100*x^3*y^8 + 36*x^2*y^8, 74*x^9*y^9 + x^8*y^2 + 22*x^5*y^4 + 7*x^4*y^8 + 2*x^4*y^3 + 71*x^4*y^2, 13*x^10*y + 95*x^7*y^2 + 18*x^7*y + 35*x^6*y^9 + 3*x^2*y^7, 79*x^9*y^9 + 98*x^8*y + 87*x^6*y^8 + 15*x^5*y^5, 99*x^10*y^7 + 49*x^8*y^7 + 100*x^7*y^9 + 7*x^6*y^5 + 83*x^6*y + 15*x^5*y^5 + 61*x^4*y^9 + 64*x^4*y^3, 73*x^9*y^2 + 79*x^7*y^10 + 12*x^7*y^9 + 67*x^7*y^8 + 15*x^4*y^6, 58*x^5*y^6, 61*x^10*y^8 + 18*x^8*y^9 + 85*x^8*y^4 + 34*x^5*y^9 + 48*x^3*y^7 + 9*x^3*y^3 + 80*x*y^4]
#monomialSupport = unique(sort(vcat(collect.(monomials.(system))...)))
#conv_hulls=[convex_hull(collect(Nemo.exponent_vectors(f))) for f in system]


# TODO: Deal with polynomial ring in a nonglobal way?
function greedy_vertex_alignment(system)
    # sorting matters for greedy algo, consider the best one TODO
    vertices=[collect(Nemo.exponent_vectors(f)) for f in system]
    # sort vertices and then the system by ascending number of monomials
    perm=sortperm(vertices)
    vertices=vertices[perm]
    system=system[perm]
    # iterate through the monomials of system i compared to i+1
    # find optimal path between these, greedy choice is to combine these all
    # if m is max monomial length, and n the number of polynomials, max paths computed is (n-1)*m^2
    # going through all of them would yield m^n paths

    # leave last fixed, we transform the rest
    monomialTranslations=[]
    for i in 1:(length(system)-1)
        # score is ratio of zero minors compared to all minors
        max_score=-1
        max_monomial=0
        for v1 in vertices[i]
            for v2 in vertices[i+1]
                monomial=prod(gens(S).^(v2-v1))
                # doesn't work properly TODO
                lowest_term=sort(v2-v1)[1]+100
                mod_system=system[:]
                mod_system=push!(deleteat!(mod_system,i),monomial*system[i])
                global_translation=prod(gens(S).^[lowest_term for j in gens(S)])
                mod_system=[f*global_translation for f in mod_system]
                # computing the score this way is far too slow
                # order of computing minors is high
                # what should score be? currently just number of nontrivial zero minors TODO
                score=number_of_zero_minors(matrix_from_system([evaluate(f,gens(R)) for f in mod_system]))
                if score>max_score
                    print(score)
                    max_score=score
                    max_monomial=monomial
                end
                println("iteration done")
            end
        end
        push!(monomialTranslations,max_monomial)
    end
    push!(monomialTranslations,1)
    println(system[1])
    println(monomialTranslations)
    return [prod(monomialTranslations[i:length(monomialTranslations)])*system[i] for i in 1:length(system)]
end


# return a list of vectors paired with a point
function convex_hull_edge_vectors(convs)
    i=0
    listofedges=[]
    for conv in convs
        i+=1
        # structure of an edge element: vector describing the edge, point of origin, and the i-th system its associated to (needed for applying the monomial map later)
        push!(listofedges,[[m[1]-m[2],m[2],i] for m in vertices.(faces(conv,1))])
    end
    return listofedges
end

function aligned_edges(listofedges)
    #maxdim=length(listofedges)
    # we choose k elements from the list
    #nchoosek=[collect(AbstractAlgebra.combinations(listofedges,k)) for k in 2:maxdim]
    # above, nchoosek is in the form [ [subsets of 2], [subsets of 3], ... [subsets of n]]. concatenate these
    #nchoosek=vcat(nchoosek...)
    #possible_combinations=[collect(Iterators.product(j...)) for j in nchoosek]
    # Iterators.product returns a matrix, convert this into a vector
    #possible_combinations=[vec(j) for j in possible_combinations]
    #possible_combinations=vcat(possible_combinations...)
    # if we have perfect matches, greedily take those first
    #perfect_alignment=filter!(m->prod([c[1]==m[1][1]||c[1]==-m[1][1] for c in m]),possible_combinations)
    #perfect_alignment=filter!(m->m!=[],perfect_alignment)
    # filter this down further by seeing which edge alignments can be attained by a monomial map
    # map to gens(S) TODO
    edges=vcat(listofedges...)
    edges=unique!([e[1] for e in edges])
    # eliminate edges that are equal up to sign TODO
    delta=[[length(filter(m->m[1]==e,sys)) for e in edges] for sys in listofedges]
   # perfect_maps=[[gens(S).^Int.(x[2]-alignment[1][2]) for x in alignment] for alignment in perfect_alignment]
    # define maps TODO
    return delta
    # otherwise, settle for parallel edges
    # TODO
end


function matrix_from_system(pol_system)
    mons=unique(collect(Iterators.flatten([collect(monomials(f)) for f in pol_system])))
    S = matrix_space(QQ, length(pol_system), length(mons))
    M_list=collect(Iterators.flatten(([[coeff(f,m) for m in mons] for f in pol_system])))
    M=matrix(QQ,length(pol_system),length(mons),M_list)
    return M
end

function is_minor_trivial(cols)
    T = tropical_semiring()
    # go through all permutations of columns, if there exists one arrangement where column i has i-th element nonzero, return true.
    incidenceMatrix = matrix(T,[ [-Int(!iszero(c)) for c in col] for col in cols])
    return det(incidenceMatrix)==T(-nrows(incidenceMatrix))
     #for m in collect(Combinatorics.permutations(cols))
         #truthlist=[iszero(l[get_position(l,m)]) for l in m]
         #if ~(true in truthlist)
             #return true
         #end
     #end
     #return false
end

# rename to ratio_of_zero_minors or something
# returns number of zero minors that we consider
function number_of_zero_minors(mat)
    cols=[mat[:,i] for i in 1:number_of_columns(mat)]
 square_dim=min(number_of_columns(mat),number_of_rows(mat))
    minor=AbstractAlgebra.combinations(cols,square_dim)
    minor=filter!(x->is_minor_trivial(x),minor)
    minor=[hcat(m...) for m in minor]
    niceminors=filter!(x->det(x)==0,minor)
    # return the ratio of (nontrivial) minors to total minors
    #return length(niceminors)/length(minors(mat,square_dim))
    return length(niceminors)
    #println(mat)
    #mnrs=minors(mat,square_dim)
    #zero_mnrs=filter!(m->m==0,mnrs)
    #println(length(zero_mnrs))
    #println(length(mnrs))
    #return length(zero_mnrs)/length(mnrs)
end
