using Oscar;

# map systems to their newton polytope
# identify common edges
# greed: prefer edges that appear in more vs appear in less
# identify translations (iso to monomials) required to attain these
# greed, prefer monomials which appear more vs less
# identify circuits of associated matrix by columns, try to maximise


# form matrix out of systems
# form n-minors out of this matrix

include("BIOMD0000001004.jl")

# paramsRing
# polRing

phi=hom(polRing,polRing,c->evaluate(c,rand(Int8,length(gens(paramsRing)))),gens(polRing))

# map parameters to random values
system=phi.(chemSystem)
# some systems have repeated equations, or zero equations, as they arise from ODEs
system=unique(system)
system=filter(l->!iszero(l),system)
system=[prod(gens(polRing).^abs.(rand(Int8,length(gens(polRing)))))*f for f in system]

#system=[x^2+y^2+2x+3y+7x^2*y^3,x^2+y^7+2x^2+3y^3+7x+2y+x^7*y^7,x^2+y^2+1]
#system=[x^5*(1+x+y+x*y),x+y+x*y,x^3+y]
#system=[58*x^4*y^6 + 97*x^3*y^9, 12*x^8*y^9 + 92*x^5*y^9 + 29*x^5*y^2 + 99*x^4*y^5 + 92*x^4*y + 76*x^2*y^3 + 44*x*y^2, 27*x^10*y^5 + 87*x^9*y^8 + 70*x^8*y^5 + 100*x^3*y^8 + 36*x^2*y^8, 74*x^9*y^9 + x^8*y^2 + 22*x^5*y^4 + 7*x^4*y^8 + 2*x^4*y^3 + 71*x^4*y^2, 13*x^10*y + 95*x^7*y^2 + 18*x^7*y + 35*x^6*y^9 + 3*x^2*y^7, 79*x^9*y^9 + 98*x^8*y + 87*x^6*y^8 + 15*x^5*y^5, 99*x^10*y^7 + 49*x^8*y^7 + 100*x^7*y^9 + 7*x^6*y^5 + 83*x^6*y + 15*x^5*y^5 + 61*x^4*y^9 + 64*x^4*y^3, 73*x^9*y^2 + 79*x^7*y^10 + 12*x^7*y^9 + 67*x^7*y^8 + 15*x^4*y^6, 58*x^5*y^6, 61*x^10*y^8 + 18*x^8*y^9 + 85*x^8*y^4 + 34*x^5*y^9 + 48*x^3*y^7 + 9*x^3*y^3 + 80*x*y^4]

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
    # our score is the norm of a vector of number of relevant zero minors
    # and the reciprocal of the number of columns
    base_score=norm([number_of_zero_minors(matrix_from_system(system)), 1/number_of_columns(matrix_from_system(system))])
    base_score=1
    for i in 1:(length(system)-1)
        max_score=base_score
        max_monomial=1
        for v1 in vertices[i]
            for v2 in vertices[i+1]
                lowest_term=sort(v2-v1)[1]
                monomial=prod(gens(polRing).^(v2-v1+[abs(lowest_term) for j in v2]))
                mod_system=system[:]
                mod_system=push!(deleteat!(mod_system,i),monomial*system[i])
                global_translation=prod(gens(polRing).^[abs(lowest_term) for j in v2])
                mod_system=[f*global_translation for f in mod_system]
                mat=matrix_from_system(mod_system)
                num=number_of_zero_minors(mat)
                len=number_of_columns(mat)
                zerominors=number_of_zero_minors(mat)
                #score=norm([zerominors,1/len])
                score=1/len
                if score>max_score
                    max_score=score
                    max_monomial=monomial
                end
            end
        end
        push!(monomialTranslations,max_monomial)
    end
    push!(monomialTranslations,1)
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

function matrix_from_system(pol_system)
    mons=unique(collect(Iterators.flatten([collect(monomials(f)) for f in pol_system])))
    S = matrix_space(QQ, length(pol_system), length(mons))
    M_list=collect(Iterators.flatten(([[coeff(f,m) for m in mons] for f in pol_system])))
    M=matrix(QQ,length(pol_system),length(mons),M_list)
    return M
end

# returns number of zero minors that we consider
function number_of_zero_minors(mat)
    # we tag each column so that we can consider columns that appear more than once
    cols=[[i,mat[:,i]] for i in 1:number_of_columns(mat)]
    cols_per_sys=[filter(m->!iszero(m[2][i]),cols) for i in 1:number_of_rows(mat)]
    square_dim=min(number_of_columns(mat),number_of_rows(mat))
    # We assume the number of rows is larger than the number of columns for now, TODO
    minors=collect(Iterators.product(cols_per_sys...))
    # delete those with repeat columns since they will always be det 0
    minors=filter(m->length(m)==length(unique(m)),minors)
    minors=[[c[2] for c in m] for m in minors]
    #minor=[hcat(m...) for m in minor]
    niceminors=filter(m->det(hcat(m...))==0,minors)
    # return the number of relevant zero minors
    return length(niceminors)
end
