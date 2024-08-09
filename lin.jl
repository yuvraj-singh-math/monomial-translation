using Oscar;
using Combinatorics;

# map systems to their newton polytope
# identify common edges
# greed: prefer edges that appear in more vs appear in less
# identify translations (iso to monomials) required to attain these
# greed, prefer monomials which appear more vs less
# identify circuits of associated matrix by columns, try to maximise


# form matrix out of systems
# form n-minors out of this matrix

S,(x,y)=polynomial_ring(QQ,["x","y"])
system=[x^2+y^2+2x+3y+7x^2*y^3,x^2+y^7+2x^2+3y^3+7x+2y+x^7*y^7,x^2+y^2+1]
conv_hulls=[convex_hull(collect(Nemo.exponent_vectors(f))) for f in system]

# return a list of vectors paired with a point
function convex_hull_edge_vectors(conv)
    # vertices.faces(conv,1) selects all edges, returned as the two points it lies between
    # we calculate the difference
    edges=[[m[1]-m[2],m[2]] for m in vertices.(faces(conv,1))]
    return edges
end

function aligned_edges(listofedges)
    maxdim=length(listofedges)
    # we choose k elements from the list
    nchoosek=[collect(Combinatorics.combinations(listofedges,k)) for k in 2:maxdim]
    possible_combinations=[collect(Iterators.product(c...)) for j in nchoosek for c in j]
    println(possible_combinations)
    # if we have perfect matches, greedily take those first
    possible_combinations=[vec(j) for j in possible_combinations]
    compared=[filter!(m->prod([c[1]==m[1][1]||c[1]==-m[1][1] for c in m]),j) for j in possible_combinations]
    compared=filter!(m->m!=[],compared)
    return compared
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

function get_position(element,list)
    n=findfirst(x->x==element,list)
    return n
end

function is_minor_trivial(cols)
    # go through all permutations of columns, if there exists one arrangement where column i has i-th element nonzero, return true.
    for m in collect(Combinatorics.permutations(cols))
        println(m)
        truthlist=[iszero(l[get_position(l,m)]) for l in m]
        if ~(true in truthlist)
            return true
        end
    end
    return false
end

# returns number of zero minors that we consider
function number_of_zero_minors(mat)
    cols=[mat[:,i] for i in 1:number_of_columns(mat)]
    square_dim=min(number_of_columns(mat),number_of_rows(mat))
    minor=collect(Combinatorics.combinations(cols,square_dim))
    minor=filter!(x->is_minor_trivial(x),minor)
    minor=[hcat(m...) for m in minor]
    number=length(filter!(x->det(x)==0,minor))
    return number
end
