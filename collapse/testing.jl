#################################################################
#                   Includes and Modules                        #
#################################################################

using LinearAlgebra, Arpack, Random, SparseArrays
using BenchmarkTools, Printf, TimerOutputs
Random.seed!(12345)

#include("../m1Accel.jl")

include("../HOLaGraf_lsmr.jl")
using .HOLaGraf_lsmr
include("../generateDelauney.jl")

# Plotting thingies 
using Plots, ColorSchemes, Plots.PlotMeasures,  LaTeXStrings
pgfplotsx()
theme(:mute)
Plots.scalefontsizes(1.25)
cols=ColorSchemes.Spectral_11;
 
# Sampling
using StatsBase
using GR: delaunay
using Distributions
import Base.Iterators: flatten
#################################################################
#                End of Includes and Modules                    #
#################################################################


#################################################################
#                   Functions of Sets                          #
#################################################################

"""
    getEdges2Trig( edges2, trigs2 )

Create a vector of sets of triangles adjacent to edges
"""
function getEdges2Trig( edges2, trigs2 )
      edg2Trig = [];
      for i in axes(edges2, 1)
            tmp = [];
            for j in axes(trigs2, 1)
                  if ((trigs2[j, 1] == edges2[i, 1]) && (trigs2[j, 2] == edges2[i, 2])) || ((trigs2[j, 1] == edges2[i, 1]) && (trigs2[j, 3] == edges2[i, 2])) || ((trigs2[j, 2] == edges2[i, 1]) && (trigs2[j, 3] == edges2[i, 2]))
                        tmp=[tmp; j];
                  end
            end
            edg2Trig = [edg2Trig; Set(tmp)];
      end
      return edg2Trig
end
      
"""
    getTrig2Edg( edges2, trigs2, edg2Trig )

create a vector of sets of edges adjacent to triangles
"""
function getTrig2Edg( edges2, trigs2, edg2Trig )
      trig2Edg = [ Set{Int64}([]) for i in axes(trigs2, 1) ];
      for i in axes(edges2, 1)
            for j in edg2Trig[i]
                  push!(trig2Edg[j], i)
            end
      end
      return trig2Edg
end
      
#################################################################
#             END of Functions of Sets                          #
#################################################################



#################################################################
#                Functions of Delaunay                          #
#################################################################

function sparseDelaunay( ; N = 10, ν = 0.4 )
      n  = N + 4;

      points, edges, trigs = generateDelauney( N )
      ν_init = size(edges, 1 ) / ( ( N + 4) * ( N + 3 ) / 2 );

      backlash = Int( - size( edges, 1 ) + round( ν * n * (n-1) / 2 ) ) 

      if backlash < 0
            for i in 1 : -backlash
                  ind = getIndx2Kill( edges ) ;
                  edges, trigs = killEdge(ind, n, edges, trigs)
            end
      else
            allEdges = Array{Integer}(undef, 0, 2)
            for i in 1:(n-1) 
                  for j in (i+1):n
                        allEdges = [ allEdges; i j  ]
                  end
            end
            for i in 1:size(edges, 1)
                  indx = findall(all(allEdges .== edges[i, :]', dims=2))[1][1]
                  allEdges = allEdges[ 1:size(allEdges, 1) .!= indx, : ]
            end
            for i in 1 : backlash      
                  ind, allEdges = getNewEdge2(n, edges, allEdges);
                  edges, trigs = addEdge(ind, n, edges, trigs)
            end
      end

      return n, points, ν_init, edges, trigs

end

#################################################################
#            END of Functions of Delaunay                       #
#################################################################




function greedyCollapse( edges, trigs )
      edge2Trig = getEdges2Trig(edges, trigs)
      trig2Edg = getTrig2Edg(edges, trigs, edge2Trig)

      Ls = [ length(edge2Trig[i]) for i in 1 : size( edges, 1 ) ]
      Free = Set([ i for i in 1 : size(edges, 1 ) if Ls[i] == 1 ])

      Σ = [];
      Τ = [ ]; 
      while !isempty(Free)
            σ = pop!(Free) 
            Σ = [ Σ; σ ]
            τ = first( edge2Trig[ σ ] ) 
            Τ = [ Τ ; τ ]
            for ξ in trig2Edg[ τ ]
                  setdiff!( edge2Trig[ ξ ], Set([ τ ]) )
                  Ls[ ξ ] = Ls[ ξ ] - 1
                  ( Ls[ ξ ] == 1 ) && ( union!(Free, Set( [ ξ ] ) ) )
                  ( Ls[ ξ ] == 0 ) && ( setdiff!(Free, Set( [ ξ ] ) ) )
            end
      end
      return sum( Ls ) == 0, Σ, Τ, edge2Trig, trig2Edg, Ls, Free  
end

function greedyCollapseShort( edge2Trig, trig2Edg, perm, edges, iperm)
      
      Ls = [ length(edge2Trig[i]) for i in 1 : size( edges, 1 ) ]
      Free = Set([ i for i in 1 : size(edges, 1 ) if Ls[i] == 1 ])

      Σ = [];
      Τ = [ ]; 
      while !isempty(Free)
            σ = pop!(Free) 
            Σ = [ Σ; σ ]
            τ = first( edge2Trig[ σ ] ) 
            Τ = [ Τ ; τ ]
            #in = findall(perm .== τ)[1][1] 
            #in = iperm[τ]
            in = τ
            for ξ in trig2Edg[ in ]
                  setdiff!( edge2Trig[ ξ ], Set([ τ ]) )
                  Ls[ ξ ] = Ls[ ξ ] - 1
                  ( Ls[ ξ ] == 1 ) && ( union!(Free, Set( [ ξ ] ) ) )
                  ( Ls[ ξ ] == 0 ) && ( setdiff!(Free, Set( [ ξ ] ) ) )
            end
      end
      return sum( Ls ) <= 0, Σ, Τ, edge2Trig, trig2Edg, Ls, Free  
end

function generateComplete( n )
      edges = Array{Integer}(undef, 0, 2)
      for i in 1 : n-1
            for j in i+1 : n
                  edges = [edges; i j ]
            end
      end
      trigs = Array{Integer}(undef, 0, 3)
      for i in 1 : n-2
            for j in i+1 : n-1
                  for k in j+1 : n
                        trigs = [ trigs; i j k ]
                  end
            end
      end

      return n, edges, trigs
end

function indicator( ind, m )
      res = zeros(m)
      res[ind] .= 1
      return diagm(res)
end



function B2fromTrig(edges, trigs)
    
      m = size(edges, 1);
      if length(trigs)==1
          return spzeros(m, 1)
      end
      del = size(trigs, 1);
      B2 = spzeros(m, del);
      
      for i in 1:del
          B2[findfirst( all( [trigs[i, 1], trigs[i, 2]]' .== edges, dims=2 )[:, 1] ), i]=1;
          B2[findfirst( all( [trigs[i, 1], trigs[i, 3]]' .== edges, dims=2 )[:, 1] ), i]=-1;
          B2[findfirst( all( [trigs[i, 2], trigs[i, 3]]' .== edges, dims=2 )[:, 1] ), i]=1;
      end 
      return B2
end
 
function Perm(v)
      m = size(v, 1)
      res = zeros( m, m )
      for i in 1 : m 
            res[ v[i], i ] = 1
      end
      return res
end

function condPlus( A; thr = 1e-4 )
      m = size(A, 1);
      #σ = svds(L1up, nsv = m - 1)[1].S;
      σ = svd(Matrix(A)).S;
      return maximum( σ ) / minimum( σ[ abs.(σ) .> thr ])
end



function cgls( A, b; tol = 1e-3, maxit = 500 )
      x0 = zeros( size(A, 2), 1 )
      
      r0 = b - A * x0
      p0 = A' * r0
      s0 = p0

      γ = norm(s0)^2
      β = 2.0
      it = 0 
      for i in 1 : maxit
            it = i
            qi = A * p0
            ai = γ / norm(qi)^2
            x0 = x0 + ai * p0
            r0 = r0 - ai * qi

            if norm(r0, Inf) < tol
                  break
            end

            s0 = A' * r0
            β = norm( s0 )^2 / γ
            γ = norm( s0 )^2

            p0 = s0 + β * p0
      end

      return x0, it 
end

function getAllEdges( n )
      allEdges = Array{Integer}(undef, 0, 2)
      for i in 1:(n-1) 
            for j in (i+1):n
                  allEdges = [ allEdges; i j  ]
            end
      end
      return allEdges
end

function repeatTries(N, add, rep)
      κ_original = zeros( rep )
      κ_precon = zeros( rep )
      κ_ldl = zeros( rep )
      it_ldl = zeros( rep )
      it_original = zeros( rep )
      it_precon = zeros( rep )
      m_sizes = zeros( rep )
      Δ_sizes = zeros( rep )

      for repIt in 1 : rep
            points, edges, trigs = generateDelauney( N )
            edges2, trigs2 = deepcopy(edges), deepcopy(trigs)
            n = N + 4

            ν_Δ = size(edges, 1) / binomial( n, 2 )

            #add = 60 # how many edges do we need to add
            
            allEdges = getAllEdges(n)
            for i in axes(edges2, 1)
                  indx = findall(all(allEdges .== edges2[i, :]', dims=2))[1][1]
                  allEdges = allEdges[ 1:size(allEdges, 1) .!= indx, : ]
            end
            for i in 1 : add      
                  ind, allEdges = getNewEdge2(n, edges2, allEdges);
                  edges2, trigs2 = addEdge(ind, n, edges2, trigs2)
            end


            w = zeros( size(trigs2, 1), 1 )
            m = size( edges2, 1 )
            Δ = size( trigs2, 1 )

            m_sizes[ repIt ] = m
            Δ_sizes[ repIt ] = Δ

            if size(Δ, 1) > 0.95*m/4*log(4*m)
                  break
            end

            edg2Trig = getEdges2Trig( edges2, trigs2 )
            trig2Edge = getTrig2Edg( edges2, trigs2, edg2Trig )

            #w_e = abs.( rand( Cauchy(0, 1) , size(edges2, 1) ) )
            #w_e[ w_e .> 3.5 ] .= 3.5
            w_e = abs.( rand( size(edges2, 1) ) )
            for i in 1 : Δ
                  w[i] = minimum( w_e[ collect( trig2Edge[i] ) ]  )
            end

            #μ1 = 1.0; s1 = 1/3; μ2 = 0.2;  s2 = 0.01;
            #=μ1 = 1.0; s1 = 1/3; μ2 = 0.5;  s2 = 0.1;
            for i in 1 : Δ
                  #global w, trigs2, trigs
                  if sum( all( trigs .== trigs2[i, :]'; dims = 2)) == 1   
                        if rand() > 3/Δ
                              w[ i ] = abs( randn()*s1 + μ1 )
                        else
                              w[ i ] = abs( randn()*s2 + μ2/2 )
                        end 
                  else
                        w[ i ] = abs( randn()*s2 + μ2 )
                  end
            end=#

            W = diagm(vec(sqrt.(w)))
            B2 = B2fromTrig( edges2, trigs2 )
            Lu =  B2 * W * W * B2'

            original = condPlus( Lu )

            perm = sortperm( w, dims =1, rev = true )

            fl, Σ1, Τ1, edge2Trig, trig2Edg, Ls, Free = greedyCollapse( edg2Trig, trig2Edge) 

            sub = [ ] 
            ind = 1

            subEdg2Trigs = [ Set([ ]) for i in 1 : m ]
            subTrig2Edg = [ ]

            while ind <= Δ
                  #global subEdg2Trigs, subTrig2Edg, trig2Edg, sub, ind, edges2
                  tmpSubEdg2Trigs = deepcopy( subEdg2Trigs )
                  tmpSubTrig2Edg = deepcopy( subTrig2Edg )
                  tmpSub = deepcopy( sub )

                  tmpSub = [ tmpSub; perm[ind] ]
                  tmpSubTrig2Edg = [ tmpSubTrig2Edg; trig2Edge[ perm[ind] ] ]
                  for e in trig2Edge[ perm[ind] ]
                        tmpSubEdg2Trigs[e] = union( tmpSubEdg2Trigs[e], perm[ind] )
                  end

                  fl, Σ, Τ, edge2Trig, trig2Edg, Ls, Free = greedyCollapseShort( tmpSubEdg2Trigs, tmpSubTrig2Edg, tmpSub, edges2 )
                  
                  if fl
                        sub = [ sub; perm[ind] ]
                        subTrig2Edg = [ subTrig2Edg; trig2Edge[ perm[ind] ] ]
                        for e in trig2Edge[ perm[ind] ]
                              subEdg2Trigs[e] = union( subEdg2Trigs[e], perm[ind] )
                        end
                  end

                  ind = ind + 1
            end

            filter = [Τ1; sub]
            Π = indicator( sort(filter), Δ )
            trigs3 = trigs2[ sort(filter), :]
            edge2Trig2 = getEdges2Trig( edges2, trigs3 )
            trig2Edge2 = getTrig2Edg( edges2, trigs3, edge2Trig2 )

            fl, Σ, Τ, edge2Trig2, trig2Edg2, Ls, Free = greedyCollapse( edges2, trigs3 )

            Σ_full = [ Σ; sort(collect(setdiff( Set(1:m), Set(Σ)))) ]
            Τ_full = [ filter[Τ]; sort(collect(setdiff( Set(1:Δ), Set(filter[Τ]))))  ]

            P1 = Perm( Σ_full )
            P2 = Perm( Τ_full )

            C = P1 * B2 * W * Π * P2
            
            Lu2 = pinv(C) * P1 * Lu * P1' * pinv(C')
            Lu2 = Lu2[1:size(sub,1), 1:size(sub,1)]

            κ_original[ repIt ] = condPlus(Lu)
            κ_precon[ repIt ] = condPlus(Lu2)
            _, it_original[ repIt ] = cgls(Lu, Lu*ones( size(Lu, 1) ))
            _, it_precon[ repIt ] = cgls(Lu2, Lu2*ones( size(Lu2, 1) ))
      end

      return κ_original, κ_precon, κ_ldl, it_original, it_precon, it_ldl, m_sizes, Δ_sizes
end




N = 28
points, edges, trigs = generateDelauney( N )
edges2, trigs2 = deepcopy(edges), deepcopy(trigs)
n = N + 4

ν_Δ = size(edges, 1) / binomial( n, 2 )


add = 102 # how many edges do we need to add

allEdges = Array{Integer}(undef, 0, 2)
for i in 1:(n-1) 
      global allEdges, n
      for j in (i+1):n
            allEdges = [ allEdges; i j  ]
      end
end
for i in axes(edges, 1)
      global allEdges, edges
      indx = findall(all(allEdges .== edges[i, :]', dims=2))[1][1]
      allEdges = allEdges[ 1:size(allEdges, 1) .!= indx, : ]
end
for i in 1 : add      
      global n, edges, edges2, trigs2, allEdges
      ind, allEdges = getNewEdge2(n, edges, allEdges);
      edges2, trigs2 = addEdge(ind, n, edges2, trigs2)
end



w = zeros( size(trigs2, 1), 1 )
m = size( edges2, 1 )
Δ = size( trigs2, 1 )

μ1 = 1.0; s1 = 1/3; μ2 = 0.2;  s2 = 0.01;

for i in 1 : Δ
      global w, trigs2, trigs
      if sum( all( trigs .== trigs2[i, :]'; dims = 2)) == 1   
            w[ i ] = abs( randn()*s1 + μ1 )
      else
            w[ i ] = abs( randn()*s2 + μ2 )
      end
end

W = diagm(vec(sqrt.(w)))
B2 = B2fromTrig( edges2, trigs2 )
Lu =  B2 * W * W * B2'

original = condPlus( Lu )


perm = sortperm( w, dims =1, rev = true )
iperm = invperm(vec(perm))

edge2Trig = getEdges2Trig( edges2, trigs2 )
trig2Edge = getTrig2Edg( edges2, trigs2, edge2Trig )

fl, Σ1, Τ1, edge2Trig2, trig2Edg2, Ls, Free = greedyCollapseShort( edge2Trig, trig2Edge, vec(1:Δ), edges2, vec(1:Δ) )

leftovers = sort(unique(collect(flatten(collect.(edge2Trig)))))
leftoversOP = sort(iperm[leftovers])
sub = [ ] 
ind = 1

subEdg2Trigs = [ Set([ ]) for i in 1 : m ]
#subTrig2Edg = [ ]
iperm2 = zeros(Int, Δ)
for i in eachindex( leftoversOP)
#while ind <= Δ
      global subEdg2Trigs, subTrig2Edg, trig2Edg, sub,  edges2, edge2Trig, perm, leftovers, iperm2
            ind = leftoversOP[i]
      #subEdg2Trigs = [ Set([ ]) for i in 1 : m ]
      #if perm[ind] in leftovers
            sub = [ sub; perm[ind] ]
            iperm2[perm[ind]] = i
            #println(sub)
            #println(iperm2)
            #sub = [ sub; ind ]
            #subTrig2Edg = [ subTrig2Edg; trig2Edge[ perm[ind] ] ]
            #subTrig2Edg = [  trig2Edge[ sub_ind ] for sub_ind in sub ]
            #for e in trig2Edge[ perm[ind] ]
            #      union!( subEdg2Trigs[e], perm[ind] )
            #end
            for sub_ind in sub
                  for e in trig2Edge[ sub_ind ]
                        union!( subEdg2Trigs[e], sub_ind )
                  end
            end
            #println(subEdg2Trigs)

            fl, _, _, _, _, Ls, Free = greedyCollapseShort( subEdg2Trigs[:], trig2Edge[:], sub, edges2, iperm2 )
            #println(fl)
            #println()
            #@printf "length: %d   |   is it succes? %d \n " size(sub, 1) Int(fl)
            if !fl
                  pop!(sub)
                  #pop!(subTrig2Edg)
                  for e in trig2Edge[ perm[ind] ]
                        setdiff!( subEdg2Trigs[e], perm[ind] )
                  end
            end
      #end
      #ind = ind + 1
end




filter = [ Τ1; sub]
Π = indicator( sort(filter), Δ )
trigs3 = trigs2[ sort(filter), :]
edge2Trig2 = getEdges2Trig( edges2, trigs3 )
trig2Edge2 = getTrig2Edg( edges2, trigs3, edge2Trig2 )

fl, Σ, Τ, edge2Trig2, trig2Edg2, Ls, Free = greedyCollapse( edges2, trigs3 )

Σ_full = [ Σ; sort(collect(setdiff( Set(1:m), Set(Σ)))) ]
Τ_full = [ filter[Τ]; sort(collect(setdiff( Set(1:Δ), Set(filter[Τ]))))  ]

P1 = Perm( Σ_full )
P2 = Perm( Τ_full )

C = P1 * B2 * W * Π * P2
D = B2 * W * Π
E = B2 * W

fl, condPlus( Lu ), condPlus( pinv(C) * P1 * Lu * P1' * pinv(C') ), norm(sqrt.(W)*Π)/norm(sqrt.(W)) #size(nullspace( Π ), 2), Δ- size(nullspace( D ), 2), minimum(w) 

U = nullspace(Lu)

C2 = ichol(Lu + 0.001 * U * U')

condPlus( pinv(C2) * Lu * pinv(C2') )


function B1fromEdges(n, edges)
      m = size(edges, 1);
      B1 = spzeros(n, m);
      
      for i in 1:m
          B1[edges[i, 1], i] = -1;
          B1[edges[i, 2], i] = 1;
      end
      return B1
  end

B1 = B1fromEdges(n, edges2)

x = rand(size(Lu, 1))
t = @elapsed Lu * x
t = @elapsed C2 \ (Lu * (C2' \ x ))
sP1, sP1T = sparse(P1), sparse(P1') 

x=rand(160)
t = @elapsed Krylov.lsmr( C1, sP1 * ( Lu  * (sP1T * Krylov.lsmr(C1', x; rtol=1e-1)[1]) ); rtol = 1e-1 )  